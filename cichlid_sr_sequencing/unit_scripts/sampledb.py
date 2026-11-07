"""Checks and the Blumer update of the SampleDatabase. Python 3.7 compatible.

    python3 unit_scripts/sampledb.py check [--no-blumer]
        READ-ONLY. Structure, duplicates, orphans between SampleDatabase / DNAReadsDatabase / DownloadDatabase, codings, blanks, spelling variants of species names (and, with --s1, how the sheet compares with Blumer).
    python3 unit_scripts/sampledb.py blumer [--out-dir D] [--write]
        Blumer's tables are read from data/blumer/ of the repository (--s1 / --s7 take other files). Adds Blumer et al.'s information to the samples whose SampleID is a Blumer BioSample, and records the SRA species of every sample that has runs in the DownloadDatabase.
        WITHOUT --write it changes nothing and writes the reports (what would change, cell by cell) into --out-dir.

What the update does, per sample that is in Table S1 (matched on S1.Biosample = SampleID):
  Species            set to Blumer's `genus species` (overwrites; the SRA name is kept in SRASpecies, so nothing is lost)
  Sex, Ecogroup      filled if blank; a value that differs is KEPT and listed as a conflict (the sheet carries manual curation: Error, Exclude, ...). Ecogroup from S1's clade:
                     Mbuna, Shallow -> Shallow_Benthic, Deep -> Deep_Benthic, Utaka, Diplo -> Diplotaxodon, Rhampho -> Rhamphochromis, AstCal -> AC
  LabReared          0 if blank (every Blumer sample is wild-caught)
  BlumerSequenceID   S1.sequence_id (= S7.sample_id);  BlumerSpecimenID   S1.specimen_id        (filled if blank, conflicts kept)
  Inversion2/9/10/11/13   Table S7 (chr2, chr9, chr10, chr11, chr13), filled if blank, conflicts kept
For EVERY sample with runs in the DownloadDatabase (tabs Blumer, Malinsky): SRASpecies = the Organism of its runs (several different ones are joined with ' | ').
New columns are inserted before the empty 'Unnamed' columns at the end; existing columns keep their positions. Only the SampleDatabase tab is written, once, after a CSV backup, and read back cell by cell."""
import argparse, collections, csv, os, re, sys, time
import pandas as pd

REPO_DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'data', 'blumer')       # Blumer's tables live in the repository
DEFAULT_S1, DEFAULT_S7 = os.path.join(REPO_DATA, 'Table_S1.tsv'), os.path.join(REPO_DATA, 'Table_S7.tsv')
ECO_MAP = {'Mbuna': 'Mbuna', 'Shallow': 'Shallow_Benthic', 'Deep': 'Deep_Benthic', 'Utaka': 'Utaka', 'Diplo': 'Diplotaxodon', 'Rhampho': 'Rhamphochromis', 'AstCal': 'AC'}
INVERSIONS = collections.OrderedDict([('chr2', 'Inversion2'), ('chr9', 'Inversion9'), ('chr10', 'Inversion10'), ('chr11', 'Inversion11'), ('chr13', 'Inversion13')])
TEXT_NEW = ['SRASpecies', 'BlumerSequenceID', 'BlumerSpecimenID']
NUMERIC_NEW = ['Inversion2', 'Inversion9', 'Inversion11', 'Inversion13']
TOUCHED = ['Sex', 'Species', 'Ecogroup', 'LabReared', 'Inversion10'] + TEXT_NEW + NUMERIC_NEW
DEFAULT_INVERSION_SOURCE = 'Blumer et al. 2025 Table S7'
SEX_OK = {'', 'M', 'F', 'U', 'M?', 'F?'}


def norm(v):
    """'' for None, NaN and whitespace, else the stripped string."""
    if v is None:
        return ''
    try:
        if pd.isna(v):
            return ''
    except (TypeError, ValueError):
        pass
    s = str(v).strip()
    return '' if s.lower() == 'nan' else s


def same_value(a, b):
    a, b = norm(a), norm(b)
    if a == b:
        return True
    try:
        return float(a) == float(b)
    except ValueError:
        return False


# ------------------------------------------------------------------------------------------------------------------ Blumer's tables and the SRA species
def read_blumer(s1_path, s7_path):
    """S1 joined with S7 (S7.sample_id = S1.sequence_id). Returns (table, problems)."""
    s1 = pd.read_csv(s1_path, sep='\t', dtype=str, keep_default_na=False); s7 = pd.read_csv(s7_path, sep='\t', dtype=str, keep_default_na=False)
    need1, need7 = {'sequence_id', 'specimen_id', 'genus', 'species', 'sex', 'clade', 'Biosample'}, {'sample_id'} | set(INVERSIONS)
    if need1 - set(s1.columns) or need7 - set(s7.columns):
        raise SystemExit('Table S1 / S7 lack columns: %s %s' % (sorted(need1 - set(s1.columns)), sorted(need7 - set(s7.columns))))
    problems = []
    for what, ids in (('S7 sample_id not in S1 sequence_id', set(s7.sample_id) - set(s1.sequence_id)), ('S1 sequence_id not in S7 sample_id', set(s1.sequence_id) - set(s7.sample_id))):
        if ids:
            problems.append('%s: %d (e.g. %s)' % (what, len(ids), sorted(ids)[:3]))
    dup = s1.Biosample[s1.Biosample.duplicated() & ~s1.Biosample.isin(['', 'NA'])].unique()
    if len(dup):
        problems.append('Biosample listed twice in S1: %s' % list(dup)[:5])
    b = s1.merge(s7, left_on='sequence_id', right_on='sample_id', how='left')
    b['species_name'] = [' '.join((g.strip() + ' ' + s.strip()).split()) if g.strip() not in ('', 'NA') and s.strip() not in ('', 'NA') else '' for g, s in zip(b.genus, b.species)]
    b['ecogroup'] = [ECO_MAP.get(c.strip(), '') for c in b.clade]
    return b, problems


def sra_species(reads_dt, tabs):
    """{SampleID: 'Organism' or 'A | B'} from the DownloadDatabase tabs, through RunID -> SampleID of DNAReadsDatabase."""
    run2sample = dict(zip(reads_dt.RunID.astype(str), reads_dt.SampleID.astype(str)))
    seen = collections.OrderedDict()
    for tab, df in tabs.items():
        if not len(df) or 'Organism' not in df.columns:
            continue
        for run, org in zip(df['Run'].astype(str), df['Organism']):
            s, org = run2sample.get(run), norm(org)
            if s and org and org not in seen.setdefault(s, []):
                seen[s].append(org)
    return {s: ' | '.join(v) for s, v in seen.items()}, sorted(s for s, v in seen.items() if len(v) > 1)


# ------------------------------------------------------------------------------------------------------------------ the plan (pure: no sheet access)
def inversion_source_text(df):
    """The InversionSource that the existing Blumer rows carry (the most frequent non-blank value among rows with a BlumerSequenceID), so that new rows read the same; a default if there is none."""
    if 'InversionSource' not in df.columns or 'BlumerSequenceID' not in df.columns:
        return DEFAULT_INVERSION_SOURCE
    vals = [norm(v) for v, k in zip(df['InversionSource'], df['BlumerSequenceID']) if norm(v) != '' and norm(k) != '']
    return collections.Counter(vals).most_common(1)[0][0] if vals else DEFAULT_INVERSION_SOURCE


def plan(sample_dt, b, sra):
    """-> (new frame, changes, summary). changes: [{SampleID, column, old, new, action}] with action filled | overwritten | kept (differs) | skipped (...)."""
    df = sample_dt.copy().reset_index(drop=True)
    unnamed = [c for c in df.columns if str(c).startswith('Unnamed:')]
    base = [c for c in df.columns if c not in unnamed]
    for c in TEXT_NEW + NUMERIC_NEW + ['Inversion10']:
        if c not in df.columns:
            df[c] = pd.Series([''] * len(df), dtype=object) if c in TEXT_NEW else pd.Series([float('nan')] * len(df), dtype=float)
    new_cols = [c for c in TEXT_NEW + NUMERIC_NEW if c not in base]
    df = df[base + [c for c in ['Inversion10'] if c not in base] + new_cols + unnamed]
    source_text = inversion_source_text(df) if 'InversionSource' in df.columns else None
    for c in TOUCHED + (['InversionSource'] if source_text else []):
        df[c] = df[c].astype(object)
    idx = {str(s): i for i, s in enumerate(df.SampleID)}
    changes, matched = [], set()

    def put(sid, col, new, mode):
        i = idx[sid]; old = df.at[i, col]
        if same_value(old, new):
            return
        if mode == 'fill' and norm(old) != '':
            changes.append(dict(SampleID=sid, column=col, old=norm(old), new=norm(new), action='kept (differs)')); return
        df.at[i, col] = new
        changes.append(dict(SampleID=sid, column=col, old=norm(old), new=norm(new), action='filled' if norm(old) == '' else 'overwritten'))

    for r in b.itertuples():
        sid = str(r.Biosample)
        if sid not in idx:
            continue
        matched.add(sid)
        put(sid, 'Sex', r.sex.strip(), 'fill') if r.sex.strip() else None
        if r.species_name:
            put(sid, 'Species', r.species_name, 'overwrite')
        else:
            changes.append(dict(SampleID=sid, column='Species', old=norm(df.at[idx[sid], 'Species']), new='', action='skipped (Blumer species missing)'))
        if r.ecogroup:
            put(sid, 'Ecogroup', r.ecogroup, 'fill')
        else:
            changes.append(dict(SampleID=sid, column='Ecogroup', old=norm(df.at[idx[sid], 'Ecogroup']), new='', action='skipped (clade %r not mapped)' % r.clade))
        put(sid, 'LabReared', 0, 'fill')
        put(sid, 'BlumerSequenceID', r.sequence_id.strip(), 'fill'); put(sid, 'BlumerSpecimenID', r.specimen_id.strip(), 'fill') if r.specimen_id.strip() not in ('', 'NA') else None
        filled_before = len(changes)
        for chr_, col in INVERSIONS.items():
            v = norm(getattr(r, chr_, ''))
            if v != '':
                put(sid, col, int(float(v)), 'fill')
        if source_text and any(ch['action'] == 'filled' and ch['column'] in INVERSIONS.values() for ch in changes[filled_before:]):
            put(sid, 'InversionSource', source_text, 'fill')                                   # only where this run filled an inversion value
    for sid, org in sra.items():
        if sid in idx:
            put(sid, 'SRASpecies', org, 'overwrite')
    for c in NUMERIC_NEW:                                                                      # numbers where there are numbers, empty cells elsewhere
        df[c] = [float(x) if norm(x) != '' else float('nan') for x in df[c]]
    summary = collections.OrderedDict()
    for ch in changes:
        k = (ch['column'], ch['action'].split(' (')[0]); summary[k] = summary.get(k, 0) + 1
    return df, changes, dict(matched=matched, counts=summary, new_columns=new_cols, inversion_source=source_text)


def species_pairs(changes):
    c = collections.Counter((ch['old'], ch['new']) for ch in changes if ch['column'] == 'Species' and ch['action'] == 'overwritten')
    return sorted(c.items(), key=lambda kv: -kv[1])


# ------------------------------------------------------------------------------------------------------------------ writing and read-back
def _ensure_grid(fm, tab, rows, cols, log):
    try:
        ws = fm.gc.open_by_key(fm.s_ID).worksheet(tab)
        if ws.row_count < rows or ws.col_count < cols:
            log('  %s: grid enlarged from %dx%d to %dx%d' % (tab, ws.row_count, ws.col_count, max(rows, ws.row_count), max(cols, ws.col_count)))
            ws.resize(rows=max(rows, ws.row_count), cols=max(cols, ws.col_count))
    except Exception as e:
        log('  %s: could not check the grid size (%s: %s); continuing' % (tab, type(e).__name__, str(e)[:100]))


def default_writer(fm, df, log):
    from gspread_dataframe import set_with_dataframe
    for k in range(4):
        try:
            set_with_dataframe(fm.gc.open_by_key(fm.s_ID).worksheet('SampleDatabase'), df); return
        except Exception as e:
            if k == 3:
                raise
            log('  write failed (%s: %s); retrying in %d s' % (type(e).__name__, str(e)[:80], 20 * (k + 1))); time.sleep(20 * (k + 1))


def readback_problems(expected, found, limit=8):
    """Cell-by-cell comparison of what was written with what the sheet holds now (empty and NaN are the same, 2 and 2.0 are the same)."""
    out = []
    if len(expected) != len(found):
        return ['row count: wrote %d, found %d' % (len(expected), len(found))]
    if list(expected.SampleID.astype(str)) != list(found.SampleID.astype(str)):
        return ['the SampleID order changed']
    cols = [c for c in expected.columns if not str(c).startswith('Unnamed:')]
    missing = [c for c in cols if c not in found.columns]
    if missing:
        return ['columns missing after the write: %s' % missing]
    for c in cols:
        bad = [i for i, (a, b) in enumerate(zip(expected[c], found[c])) if not same_value(a, b)]
        if bad:
            out.append('%s: %d cells differ (first: row %d, wrote %r, found %r)' % (c, len(bad), bad[0] + 2, norm(expected[c].iloc[bad[0]]), norm(found[c].iloc[bad[0]])))
    return out[:limit]


# ------------------------------------------------------------------------------------------------------------------ checks (read-only)
def spelling_key(s):
    return re.sub(r"[^a-z0-9]", '', re.sub(r"\bsp\.?\b|['\"-]", '', norm(s).lower()))


def orphan_rows(sample_dt, reads_dt, sra, blumer=None):
    """Samples with reads and no SampleDatabase row: [(SampleID, runs, projects, SRA species, in Blumer's S1)]."""
    known = set(sample_dt.SampleID.astype(str)); in_s1 = set(blumer[0].Biosample) if blumer is not None else set()
    by = collections.OrderedDict()
    for sid, run, proj in zip(reads_dt.SampleID.astype(str), reads_dt.RunID.astype(str), reads_dt['ProjectID'].astype(str) if 'ProjectID' in reads_dt.columns else [''] * len(reads_dt)):
        if sid not in known:
            e = by.setdefault(sid, ([], []))
            e[0].append(run); e[1].append(proj) if proj not in e[1] else None
    return [(sid, ','.join(r), ','.join(p), sra.get(sid, ''), 'TRUE' if sid in in_s1 else 'FALSE') for sid, (r, p) in by.items()]


def run_checks(sample_dt, reads_dt, tabs, blumer=None, out=print):
    n_att = 0
    def say(ok, title, detail=''):
        nonlocal n_att
        n_att += 0 if ok else 1
        out('%-9s %s%s' % ('ok' if ok else 'ATTENTION', title, ('\n          ' + detail.replace('\n', '\n          ')) if detail else ''))
    ids = sample_dt.SampleID.astype(str)
    out('SampleDatabase: %d rows, %d columns; DNAReadsDatabase: %d rows' % (len(sample_dt), len(sample_dt.columns), len(reads_dt)))
    blank_cols = [c for c in sample_dt.columns if all(norm(x) == '' for x in sample_dt[c])]
    say(True, 'columns: %s' % [c for c in sample_dt.columns if c not in blank_cols], '%d entirely empty columns (%s...): trailing ones are harmless; headers reading "Unnamed: N" can be deleted by hand' % (len(blank_cols), blank_cols[:3]) if blank_cols else '')
    d = ids[ids.duplicated()].unique()
    say(len(d) == 0, 'duplicate SampleIDs: %d' % len(d), str(list(d)[:10]) if len(d) else '')
    rr = reads_dt.RunID.astype(str); d = rr[rr.duplicated()].unique()
    say(len(d) == 0, 'duplicate RunIDs in DNAReadsDatabase: %d' % len(d), str(list(d)[:10]) if len(d) else '')
    rs = set(reads_dt.SampleID.astype(str)); orphans = sorted(rs - set(ids))
    say(len(orphans) == 0, 'samples with reads but no SampleDatabase row: %d' % len(orphans), 'e.g. %s' % orphans[:12] if orphans else '')
    noreads = sample_dt[~ids.isin(rs)]; eco = collections.Counter(norm(x) or '(blank)' for x in noreads.Ecogroup) if 'Ecogroup' in noreads else {}
    say(True, 'SampleDatabase rows without reads: %d' % len(noreads), 'by Ecogroup: %s' % dict(eco) if len(noreads) else '')
    allruns = set(rr)
    for tab, df in tabs.items():
        if not len(df) or 'DownloadStatus' not in df.columns:
            continue
        done = df[df.DownloadStatus == 'in_database']; miss = [r for r in done.Run if str(r) not in allruns]
        say(len(miss) == 0, 'DownloadDatabase %s: %d runs say in_database, %d of them are not in DNAReadsDatabase' % (tab, len(done), len(miss)), str(miss[:10]) if miss else '')
        owners = {str(r): None for r in done.Run}; wrong = []
        sample_of = dict(zip(reads_dt.RunID.astype(str), reads_dt.SampleID.astype(str)))
        for r in done.itertuples():
            s = sample_of.get(str(r.Run)); ok = {norm(getattr(r, 'BioSample', '')), norm(df.loc[r.Index, 'Library Name'] if 'Library Name' in df.columns else '')}
            if s is not None and s not in ok:
                wrong.append((str(r.Run), s))
        say(len(wrong) == 0, 'DownloadDatabase %s: runs whose SampleID is neither the BioSample nor the library name: %d' % (tab, len(wrong)), str(wrong[:6]) if wrong else '')
    for col, allowed in (('Sex', SEX_OK), ('LabReared', {'', '0', '1'})):
        if col in sample_dt.columns:
            vals = collections.Counter(norm(x) if col == 'Sex' else (str(int(float(norm(x)))) if norm(x) not in ('',) and re.match(r'^[0-9.]+$', norm(x)) else norm(x)) for x in sample_dt[col])
            odd = {k: v for k, v in vals.items() if k not in allowed}
            say(not odd, '%s values: %s' % (col, dict(vals)), 'unexpected: %s' % odd if odd else '')
    if 'Ecogroup' in sample_dt.columns:
        say(True, 'Ecogroup values: %s' % dict(collections.Counter(norm(x) or '(blank)' for x in sample_dt.Ecogroup)))
    have = sample_dt[ids.isin(rs)]
    blanks = {c: sum(1 for x in have[c] if norm(x) == '') for c in ('Sex', 'Species', 'Ecogroup', 'LabReared') if c in have.columns}
    say(not any(blanks.values()), 'blank key fields among the %d samples that have reads: %s' % (len(have), blanks))
    if 'Species' in sample_dt.columns:
        groups = collections.defaultdict(collections.Counter)
        for s in sample_dt.Species:
            if norm(s):
                groups[spelling_key(s)][norm(s)] += 1
        multi = {k: dict(v) for k, v in groups.items() if len(v) > 1}
        say(not multi, 'species names: %d distinct; %d taxa written in more than one way' % (sum(len(v) for v in groups.values()), len(multi)), '\n'.join('%s' % v for v in list(multi.values())[:8]) if multi else '')
    if blumer is not None:
        b, problems = blumer
        say(not problems, "Blumer's tables join", '; '.join(problems) if problems else '')
        known = {str(s) for s in ids}; um = b[~b.Biosample.isin(known)]
        say(len(um) == 0, 'Blumer samples (Table S1) that are not in the SampleDatabase: %d of %d' % (len(um), len(b)), 'e.g. %s' % list(um.Biosample)[:8] if len(um) else '')
        m = b[b.Biosample.isin(known)].merge(sample_dt.drop_duplicates('SampleID'), left_on='Biosample', right_on='SampleID', suffixes=('', '_sheet'))
        say(True, 'matched samples: %d; Species differs from Blumer for %d, Sex for %d, Ecogroup for %d' % (len(m), sum(not same_value(a, c) for a, c in zip(m.species_name, m.Species)), sum(not same_value(a, c) for a, c in zip(m.sex, m.Sex)), sum(not same_value(a, c) for a, c in zip(m.ecogroup, m.Ecogroup))))
    out('%d check(s) need attention' % n_att)
    return n_att


# ------------------------------------------------------------------------------------------------------------------ command line
def write_reports(out_dir, changes, summary, pairs, unmatched, not_in_s1, multi, log):
    os.makedirs(out_dir, exist_ok=True)
    def tsv(name, rows, cols):
        with open(os.path.join(out_dir, name), 'w', newline='') as f:
            w = csv.writer(f, delimiter='\t'); w.writerow(cols); [w.writerow([r.get(c, '') if isinstance(r, dict) else r[i] for i, c in enumerate(cols)]) for r in rows]
    tsv('sampledb_changes.tsv', changes, ['SampleID', 'column', 'old', 'new', 'action'])
    tsv('sampledb_species_pairs.tsv', [(n, o, nw) for (o, nw), n in pairs], ['samples', 'SRA_or_old_species', 'Blumer_species'])
    tsv('sampledb_blumer_unmatched.tsv', [(x,) for x in unmatched], ['S1_Biosample_not_in_sheet'])
    tsv('sampledb_blumer_runs_not_in_S1.tsv', [(x,) for x in not_in_s1], ['SampleID_with_Blumer_runs_not_in_S1'])
    tsv('sampledb_multiple_sra_species.tsv', [(x,) for x in multi], ['SampleID'])


def main(argv=None, fm=None, writer=default_writer, out=print):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest='cmd', required=True)
    c = sub.add_parser('check'); c.add_argument('--s1', default=DEFAULT_S1); c.add_argument('--s7', default=DEFAULT_S7); c.add_argument('--no-blumer', action='store_true', help='do not compare with Blumer et al.'); c.add_argument('--out-dir', default=None); c.add_argument('--tabs', nargs='+', default=['Blumer', 'Malinsky'])
    b = sub.add_parser('blumer'); b.add_argument('--s1', default=DEFAULT_S1); b.add_argument('--s7', default=DEFAULT_S7); b.add_argument('--tabs', nargs='+', default=['Blumer', 'Malinsky'])
    b.add_argument('--out-dir', default=None); b.add_argument('--write', action='store_true')
    a = ap.parse_args(argv)
    if fm is None:
        sys.path.insert(0, os.getcwd())
        from helper_modules.file_manager import FileManager
        fm = FileManager()
    fm.readSampleDatabase()
    tabs = {}
    for t in a.tabs:
        try:
            tabs[t] = fm.readDownloadTab(t)
        except AttributeError:
            raise SystemExit('this file_manager.py has no readDownloadTab (the download database methods): use the checkout that has them (mzebra)')
    if a.cmd == 'check':
        blumer = read_blumer(a.s1, a.s7) if (not a.no_blumer and os.path.exists(a.s1) and os.path.exists(a.s7)) else None
        n = run_checks(fm.sample_dt, fm.reads_dt, tabs, blumer=blumer, out=out)
        orph = orphan_rows(fm.sample_dt, fm.reads_dt, sra_species(fm.reads_dt, tabs)[0], blumer=blumer)
        if orph:
            od = a.out_dir or '%ssampledb_reports_%s/' % (getattr(fm, 'localTempDir', './'), time.strftime('%Y%m%d-%H%M%S')); os.makedirs(od, exist_ok=True)
            with open(od + 'sampledb_orphans.tsv', 'w', newline='') as f:
                w = csv.writer(f, delimiter='\t'); w.writerow(['SampleID', 'RunIDs', 'ProjectIDs', 'SRA_species', 'in_Blumer_S1']); w.writerows(orph)
            out('%d samples with reads but no SampleDatabase row, listed with their runs, projects and SRA species in %ssampledb_orphans.tsv' % (len(orph), od))
        return n
    for f in (a.s1, a.s7):
        if not os.path.exists(f):
            raise SystemExit('%s not found: Blumer\'s tables are in data/blumer/ of the repository (git pull?), or pass --s1 / --s7' % f)
    blumer, problems = read_blumer(a.s1, a.s7)
    for p in problems:
        out('NOTE: ' + p)
    sra, multi = sra_species(fm.reads_dt, tabs)
    before = fm.sample_dt.copy()
    new, changes, summary = plan(before, blumer, sra)
    pairs = species_pairs(changes); known = set(before.SampleID.astype(str))
    unmatched = [s for s in blumer.Biosample if s not in known]
    runs_in_s1 = set(blumer.Biosample); blumer_tab = tabs.get('Blumer')
    run2sample = dict(zip(fm.reads_dt.RunID.astype(str), fm.reads_dt.SampleID.astype(str)))
    not_in_s1 = sorted({run2sample[str(r)] for r in (blumer_tab.Run if blumer_tab is not None and len(blumer_tab) else []) if str(r) in run2sample and run2sample[str(r)] not in runs_in_s1})
    out_dir = a.out_dir or '%ssampledb_reports_%s/' % (getattr(fm, 'localTempDir', './'), time.strftime('%Y%m%d-%H%M%S'))
    write_reports(out_dir, changes, summary, pairs, unmatched, not_in_s1, multi, out)
    out('samples in Table S1: %d, matched to the SampleDatabase: %d, not in the sheet: %d; %d sheet samples have Blumer-project runs but are not in S1' % (len(blumer), len(summary['matched']), len(unmatched), len(not_in_s1)))
    out('new columns: %s' % (summary['new_columns'] or 'none'))
    if summary['inversion_source']:
        out('InversionSource of the rows filled now: %r (the value of the existing Blumer rows)' % summary['inversion_source'])
    for (col, act), n in sorted(summary['counts'].items()):
        out('  %-18s %-12s %5d' % (col, act, n))
    out('Species: %d distinct (SRA or old name -> Blumer name) pairs, most frequent: %s' % (len(pairs), '; '.join('%s -> %s (%d)' % (o, n, k) for (o, n), k in pairs[:4])))
    out('SRASpecies for %d samples (%d have more than one organism among their runs)' % (len(sra), len(multi)))
    out('reports in %s' % out_dir)
    if not a.write:
        out('report only: nothing written (add --write to apply)'); return new
    stamp = time.strftime('%Y%m%d-%H%M%S'); bdir = '%sdatabase_backup_%s/' % (getattr(fm, 'localTempDir', './'), stamp); os.makedirs(bdir, exist_ok=True)
    before.to_csv(bdir + 'SampleDatabase.csv', index=False); out('  the SampleDatabase tab was saved to %s (Google Sheets also keeps its version history)' % bdir)
    _ensure_grid(fm, 'SampleDatabase', len(new) + 1, len(new.columns), out)
    writer(fm, new, out)
    fm.readSampleDatabase(); bad = readback_problems(new, fm.sample_dt)
    if bad:
        raise SystemExit('the sheet does not hold what was written:\n  ' + '\n  '.join(bad) + '\nrestore from %s (or the version history)' % bdir)
    out('  verified cell by cell: %d rows, %d columns' % (len(new), len([c for c in new.columns if not str(c).startswith('Unnamed:')])))
    return new


if __name__ == '__main__':
    main()
