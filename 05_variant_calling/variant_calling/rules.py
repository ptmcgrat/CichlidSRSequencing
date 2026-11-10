"""The cohort check (module 05): the rules a cohort must pass before any calling, by kind of cohort (SPEC decision 61).

Reads only the SampleDatabase and the AlignmentDatabase (the ledger), so it needs no BAM and costs seconds, where a call costs hours.
`family` and `sex_set` follow sexloci 0.1.26's preflight (module 06 applies it again to the VCF), with these deliberate differences:
  * a sample with no alignment record BLOCKS here (sexloci only notes it): 05 calls from the BAM;
  * a `family` cohort BLOCKS when a listed parent is not in it, and when no member lists parents (sexloci warns: it can still use the offspring);
`bc_parents` is the cohort of a CROSS (two species, F1s), which sexloci refuses; its VCF is for module 07, and it has its own rules.
Levels: BLOCK (no call is started), WARN (starts, the manifest records it), INFO."""
import os
from collections import Counter

import pandas as pd
import yaml

from cichlid_common.common_light import config_path
from variant_calling.cohort import blank, empty_selectors, norm_key, norm_species, parse_parents, parse_sex, resolve_members, validate

DEFAULT_RULES = dict(min_per_sex_block=3, min_per_sex_warn=6, family_min_offspring_warn=12, coverage_offspring_warn=2.0, coverage_other_warn=5.0, coverage_unusable=3.0,
                     coverage_parent_warn=10.0, bc_parents_min_coverage=5.0, confound_columns=['BroodID', 'Project', 'LabReared'])


def load_rules(path=None):
    rules = dict(DEFAULT_RULES)
    p = path or config_path('cohort_rules.yaml')
    if os.path.exists(p):
        rules.update(yaml.safe_load(open(p)) or {})
    return rules


def short(ids, n=8):
    ids = list(ids)
    return ', '.join(map(str, ids[:n])) + (f' ... (+{len(ids) - n})' if len(ids) > n else '')


class Finding:
    def __init__(self, level, code, message, samples=None):
        self.level, self.code, self.message, self.samples = level, code, message, list(samples or [])

    def to_dict(self):
        return dict(level=self.level, code=self.code, message=self.message, samples=self.samples)


class Report:
    def __init__(self, name, kind, genome):
        self.name, self.kind, self.genome, self.members, self.findings = name, kind, genome, [], []

    def add(self, level, code, message, samples=None):
        self.findings.append(Finding(level, code, message, samples))

    @property
    def blocked(self):
        return any(f.level == 'BLOCK' for f in self.findings)

    def count(self, level):
        return sum(1 for f in self.findings if f.level == level)

    def codes(self, level=None):
        return [f.code for f in self.findings if level in (None, f.level)]

    def to_dict(self):
        return dict(cohort=self.name, kind=self.kind, genome=self.genome, n_members=len(self.members), members=self.members, blocked=self.blocked,
                    n_block=self.count('BLOCK'), n_warn=self.count('WARN'), n_info=self.count('INFO'), findings=[f.to_dict() for f in self.findings])

    def text(self):
        out = [f'cohort {self.name} ({self.kind}, {self.genome}): {len(self.members)} samples; {self.count("BLOCK")} BLOCK, {self.count("WARN")} WARN, {self.count("INFO")} INFO']
        for lvl in ('BLOCK', 'WARN', 'INFO'):
            for f in self.findings:
                if f.level == lvl:
                    out.append(f'[{lvl}] {f.code}: {f.message}')
        out.append('RESULT: ' + ('BLOCKED: nothing is called until the BLOCK lines are fixed' if self.blocked else 'ok to build'))
        return '\n'.join(out)


# ------------------------------------------------------------------ the pieces
def _pedigree(members, rows):
    """Families of the cohort from the Parents column: {(sire, dam): [offspring]}, the listed parents present in the cohort, and the problems."""
    known = set(rows.index)
    fams, malformed, partial, missing = {}, [], {}, {}
    mset = set(members)
    for s in members:
        st, sire, dam = parse_parents(rows.at[s, 'Parents'] if 'Parents' in rows.columns else None, known)
        if st == 'partial_single':                                    # one ID: its sex says whether it is the sire or the dam
            sx = parse_sex(rows.at[sire, 'Sex']) if sire in rows.index else None
            st, sire, dam = ('partial', sire, None) if sx == 1 else ('partial', None, sire) if sx == 0 else ('malformed', None, None)
        if st == 'malformed':
            malformed.append(s)
        elif st in ('ok', 'partial'):
            fams.setdefault((sire, dam), []).append(s)
            if st == 'partial':
                partial.setdefault((sire, dam), []).append(s)
            for p in (sire, dam):
                if p and p not in mset:
                    missing.setdefault(p, []).append(s)
    parents = {p for (a, b) in fams for p in (a, b) if p and p in mset}
    return dict(families=fams, parents=parents, malformed=malformed, partial=partial, missing_parents=missing)


def _common(R, spec, members, sample_dt, ledger, rules):
    table = sample_dt.dropna(how='all')
    miss = [s for s in members if s not in set(table.SampleID.astype(str))]
    if miss:
        keys = {}
        for t in table.SampleID.astype(str):
            keys.setdefault(norm_key(t), []).append(t)
        near = {m: keys[norm_key(m)][0] for m in miss if len(keys.get(norm_key(m), [])) == 1}
        hint = (' Likely matches (same letters and digits): ' + '; '.join(f'{m} -> {t}' for m, t in list(near.items())[:5]) + '.') if near else ''
        R.add('BLOCK', 'sample_not_in_table', f'{len(miss)} of {len(members)} members have no row in the SampleDatabase: {short(miss)}.{hint}', miss)
    dup = [s for s, n in Counter(table.SampleID.astype(str)).items() if n > 1 and s in set(members)]
    if dup:
        R.add('BLOCK', 'duplicate_rows', f'{len(dup)} members have more than one row in the SampleDatabase: {short(dup)}.', dup)
    rows = table.drop_duplicates('SampleID').set_index('SampleID')
    present = [s for s in members if s in rows.index]
    flagged = [s for s in present if str(rows.at[s, 'Ecogroup'] if 'Ecogroup' in rows.columns else '') in ('Error', 'Exclude')]
    if flagged:
        R.add('WARN', 'flagged_ecogroup', f"{len(flagged)} members are marked Ecogroup Error or Exclude in the sheet (included by name): {short(flagged)}.", flagged)
    # ---- the ledger: 05 calls from the BAM, so a member needs an alignment to this genome
    led = ledger.dropna(how='all').drop_duplicates('SampleID', keep='last').set_index('SampleID') if ledger is not None and len(ledger) else pd.DataFrame()
    noal = [s for s in members if s not in led.index]
    if noal:
        R.add('BLOCK', 'no_alignment', f"{len(noal)} of {len(members)} members have no alignment row in the AlignmentDatabase of {R.genome} (05 calls from the BAM): {short(noal)}.", noal)
    if 'GenomeVersion' in led.columns:
        bad = [s for s in members if s in led.index and not blank(led.at[s, 'GenomeVersion']) and str(led.at[s, 'GenomeVersion']) != R.genome]
        if bad:
            found = Counter(str(led.at[s, 'GenomeVersion']) for s in bad)
            R.add('BLOCK', 'genome_version_mismatch', f"{len(bad)} members were aligned to another genome than {R.genome} (found {dict(found)}): {short(bad)}.", bad)
    cov = pd.to_numeric(led['Coverage'], errors='coerce') if 'Coverage' in led.columns else pd.Series(dtype=float)
    return rows, present, led, cov


def _sex_and_species(R, spec, members, rows, present):
    r = rows.loc[present]
    bad = [s for s in present if parse_sex(r.at[s, 'Sex']) is None]
    if bad:
        vals = Counter('blank' if blank(r.at[s, 'Sex']) else repr(str(r.at[s, 'Sex'])) for s in bad)
        R.add('BLOCK', 'sex_not_mf', f"{len(bad)} of {len(present)} members have no valid sex ({', '.join(f'{k}: {v}' for k, v in vals.items())}); every member, parents included, must be exactly M or F: {short(bad)}.", bad)
    sp = {s: norm_species(r.at[s, 'Species']) for s in present}
    nosp = [s for s, v in sp.items() if not v]
    if nosp:
        R.add('BLOCK', 'species_missing', f'{len(nosp)} members have no Species: {short(nosp)}.', nosp)
    distinct = Counter(v for v in sp.values() if v)
    if len(distinct) > 1:
        R.add('BLOCK', 'species_mismatch', 'the cohort holds ' + str(len(distinct)) + ' species (' + '; '.join(f'{k}: {v}' for k, v in distinct.items()) + '): one VCF is one species for module 06; split it, or use kind bc_parents for a cross.', present)
    elif spec.get('species') and distinct and next(iter(distinct)) != norm_species(spec['species']):
        R.add('BLOCK', 'species_not_as_declared', f"the cohort declares {spec['species']!r} but its members are {next(iter(distinct))!r}.", present)


def _per_sex(R, rules, members, rows, present, parents):
    analysed = [s for s in present if s not in parents]
    sx = [parse_sex(rows.at[s, 'Sex']) for s in analysed]
    nm, nf = sum(1 for x in sx if x == 1), sum(1 for x in sx if x == 0)
    if min(nm, nf) < rules['min_per_sex_block']:
        R.add('BLOCK', 'too_few_per_sex', f"only {nm} males and {nf} females are analysed (at least {rules['min_per_sex_block']} of each; parents are not counted).")
    elif min(nm, nf) < rules['min_per_sex_warn']:
        R.add('WARN', 'low_per_sex', f"{nm} males and {nf} females: power is low (at least {rules['min_per_sex_warn']} of each recommended).")
    return analysed


def _coverage(R, rules, members, cov, fams, parents, analysed):
    nocov = [s for s in members if s in set(cov.index) and cov.get(s) != cov.get(s)]
    offspring = {s for v in fams.values() for s in v}
    c = lambda s: float(cov[s]) if s in cov.index and cov[s] == cov[s] else None
    lo = [s for s in offspring if c(s) is not None and c(s) < rules['coverage_offspring_warn']]
    if lo:
        R.add('WARN', 'low_coverage_offspring', f"{len(lo)} family offspring are below {rules['coverage_offspring_warn']:g}x: " + short(f'{s} ({c(s):.1f}x)' for s in lo), lo)
    lp = [s for s in parents if c(s) is not None and c(s) < rules['coverage_parent_warn']]
    if lp:
        R.add('WARN', 'low_coverage_parent', f"{len(lp)} parents are below {rules['coverage_parent_warn']:g}x: " + short(f'{s} ({c(s):.1f}x)' for s in lp), lp)
    other = [s for s in analysed if s not in offspring]
    lw = [s for s in other if c(s) is not None and c(s) < rules['coverage_other_warn']]
    if lw:
        unusable = [s for s in lw if c(s) < rules['coverage_unusable']]
        R.add('WARN', 'low_coverage', f"{len(lw)} fish are below {rules['coverage_other_warn']:g}x" + (f' ({len(unusable)} below {rules["coverage_unusable"]:g}x, unusable for the sex contrasts)' if unusable else '') + ': ' + short(f'{s} ({c(s):.1f}x)' for s in lw), lw)
    nv = [s for s in members if s in cov.index and cov[s] != cov[s]]
    if nv:
        R.add('INFO', 'no_coverage_value', f'{len(nv)} members have an alignment row without a Coverage value: {short(nv)}.', nv)


# ------------------------------------------------------------------ the kinds
def _brood_consistency(R, members, rows):
    """A BroodID is a family of two parents: its children must all name the same sire and dam (SPEC module 05; the sheet's BroodID and Parents columns)."""
    by_brood = {}
    for s in members:
        b = rows.at[s, 'BroodID'] if 'BroodID' in rows.columns else ''
        if not blank(b):
            by_brood.setdefault(str(b), []).append(s)
    for b, kids in sorted(by_brood.items()):
        pairs, none = Counter(), []
        for s in kids:
            st, sire, dam = parse_parents(rows.at[s, 'Parents'] if 'Parents' in rows.columns else None)
            if st == 'none':
                none.append(s)
            else:
                pairs[(sire, dam)] += 1
        if len(pairs) > 1:
            R.add('WARN', 'brood_parents_inconsistent', f"brood {b}: its children name {len(pairs)} different parent pairs ({'; '.join(f'{x} x {y}: {n}' for (x, y), n in pairs.items())}), but a BroodID is one family of two parents.", kids)
        if none:
            R.add('WARN', 'brood_member_without_parents', f"brood {b}: {len(none)} of {len(kids)} children list no parents in the Parents column: {short(none)}.", none)


def _family(R, spec, rules, members, rows, present, cov):
    ped = _pedigree(members, rows)
    _brood_consistency(R, members, rows)
    fams, parents = ped['families'], ped['parents']
    if not fams:
        R.add('BLOCK', 'no_pedigree', "no member lists parents ('sire,dam' in the Parents column): a family cohort needs offspring with their parents (sexloci would only warn).")
    for p, offs in sorted(ped['missing_parents'].items()):
        R.add('BLOCK', 'parent_not_in_cohort', f'parent {p} of {len(offs)} offspring ({short(offs, 4)}) is not in the cohort (sexloci would leave the brood out of the transmission test; a family cohort exists to give it the parents).', [p] + offs)
    if ped['partial']:
        R.add('WARN', 'parent_unknown', f"{sum(len(v) for v in ped['partial'].values())} offspring list only one parent: " + short(s for v in ped['partial'].values() for s in v), [s for v in ped['partial'].values() for s in v])
    if ped['malformed']:
        R.add('WARN', 'parents_malformed', f"the Parents field is not 'sire,dam' for {len(ped['malformed'])} members (treated as no pedigree): " + short(ped['malformed']), ped['malformed'])
    for (sire, dam), offs in sorted(fams.items(), key=lambda kv: str(kv[0])):
        if sire in rows.index and parse_sex(rows.at[sire, 'Sex']) == 0 or dam in rows.index and parse_sex(rows.at[dam, 'Sex']) == 1:
            R.add('WARN', 'parent_sex_mismatch', f'brood {sire} x {dam}: the first parent should be male (sire) and the second female (dam).', [p for p in (sire, dam) if p])
        if len(offs) < rules['family_min_offspring_warn']:
            R.add('WARN', 'small_family', f"brood {sire} x {dam} has {len(offs)} offspring (at least {rules['family_min_offspring_warn']} recommended).", offs)
    analysed = _per_sex(R, rules, members, rows, present, parents)
    _coverage(R, rules, members, cov, fams, parents, analysed)


def _related_pairs(members, rows):
    """Pairs of members that the sheet says are relatives: siblings (same non-blank BroodID, or the same sire x dam) and parent with offspring."""
    pairs = set()
    mset = set(members)
    brood = {}
    for s in members:
        if s in rows.index:
            b = rows.at[s, 'BroodID'] if 'BroodID' in rows.columns else None
            if not blank(b):
                brood.setdefault(str(b), []).append(s)
            st, sire, dam = parse_parents(rows.at[s, 'Parents'] if 'Parents' in rows.columns else None)
            if st in ('ok', 'partial'):
                for p in (sire, dam):
                    if p in mset:
                        pairs.add(tuple(sorted((s, p))))
                if st == 'ok':
                    brood.setdefault(f'{sire}x{dam}', []).append(s)
    for v in brood.values():
        v = sorted(set(v))
        pairs.update((a, b) for i, a in enumerate(v) for b in v[i + 1:])
    return sorted(pairs)


def _sex_set(R, spec, rules, members, rows, present, cov):
    ped = _pedigree(members, rows)
    analysed = _per_sex(R, rules, members, rows, present, ped['parents'])
    _coverage(R, rules, members, cov, ped['families'], ped['parents'], analysed)
    lab = [s for s in analysed if str(rows.at[s, 'LabReared'] if 'LabReared' in rows.columns else '').strip().lower() in ('1', '1.0', 'true', 'yes')
           and blank(rows.at[s, 'BroodID'] if 'BroodID' in rows.columns else '') and blank(rows.at[s, 'Parents'] if 'Parents' in rows.columns else '')]
    if lab:
        R.add('INFO', 'relatedness_unknown', f"{len(lab)} of {len(analysed)} analysed fish are lab-reared with no BroodID or Parents in the sheet: the sheet cannot rule out relatives, and inbred stock makes them likely "
              '(kinship from the genotypes is planned for after calling).', lab)
    rel = _related_pairs(members, rows)
    if rel:
        inv = Counter(x for p in rel for x in p)
        R.add('WARN', 'related_members', f"the sheet lists {len(rel)} pairs of relatives (siblings or parent and offspring) in a pool treated as unrelated animals; module 06 counts the pool as one independent unit, and inbred stock makes close kin more likely: "
               + short([f'{a}/{b}' for a, b in rel], 5), sorted(inv))
    for col in rules.get('confound_columns') or []:
        if col not in rows.columns:
            continue
        sex = {s: parse_sex(rows.at[s, 'Sex']) for s in analysed}
        lev = {}
        for s in analysed:
            v = rows.at[s, col]
            if not blank(v) and sex[s] is not None:
                lev.setdefault(str(v), []).append(sex[s])
        multi = {k: v for k, v in lev.items() if len(v) >= 2}
        if len(multi) >= 2 and all(len(set(v)) == 1 for v in multi.values()) and len({x for v in multi.values() for x in v}) == 2:
            R.add('WARN', 'sex_confounded', f'sex is perfectly confounded with {col}: every {col} level that has two or more fish holds one sex only ({short(multi)}).')


def _bc_parents(R, spec, rules, members, rows, present, cov):
    ped = _pedigree(members, rows)
    for p, offs in sorted(ped['missing_parents'].items()):
        R.add('BLOCK', 'pedigree_incomplete', f'{short(offs, 4)} list {p} as a parent, but {p} is not in the cohort (the BC cross needs the pedigree inside it).', [p] + offs)
    founders = [s for s in present if parse_parents(rows.at[s, 'Parents'] if 'Parents' in rows.columns else None)[0] == 'none']
    species = Counter(norm_species(rows.at[s, 'Species']) or '(blank)' for s in present)
    R.add('INFO', 'cohort_summary', f"{len(founders)} founders (no parents listed) and {len(present) - len(founders)} with parents; species: " + '; '.join(f'{k}: {v}' for k, v in species.items()) + '.')
    nosex = [s for s in present if parse_sex(rows.at[s, 'Sex']) is None]
    if nosex:
        R.add('WARN', 'sex_unknown', f'{len(nosex)} members have no valid sex (not needed to call, needed later): {short(nosex)}.', nosex)
    c = lambda s: float(cov[s]) if s in cov.index and cov[s] == cov[s] else None
    low = [s for s in members if c(s) is not None and c(s) < rules['bc_parents_min_coverage']]
    if low:
        R.add('BLOCK', 'coverage_too_low', f"{len(low)} members are below {rules['bc_parents_min_coverage']:g}x, too shallow to call genotypes that define informative sites (a provisional threshold, an open decision): " + short(f'{s} ({c(s):.1f}x)' for s in low), low)
    warn = [s for s in members if c(s) is not None and rules['bc_parents_min_coverage'] <= c(s) < rules['coverage_parent_warn']]
    if warn:
        R.add('WARN', 'low_coverage', f"{len(warn)} members are below {rules['coverage_parent_warn']:g}x: " + short(f'{s} ({c(s):.1f}x)' for s in warn), warn)


# ------------------------------------------------------------------ the entry point
def check_cohort(name, spec, sample_dt, ledger_dt, rules=None):
    """The Report of cohort `name` (a definition from load_cohorts) against the SampleDatabase and AlignmentDatabase DataFrames."""
    R = Report(name, spec.get('kind'), spec.get('genome'))
    problems = validate(spec)
    if problems:
        for p in problems:
            R.add('BLOCK', 'bad_definition', p)
        return R
    rules = rules or load_rules()
    members, notes = resolve_members(sample_dt, spec)
    R.members = members
    for n in notes:
        R.add('INFO', 'selection', n)
    for msg in empty_selectors(sample_dt, spec):
        R.add('WARN', 'selector_empty', msg)
    if not members:
        R.add('BLOCK', 'no_members', 'the selectors match no sample.')
        return R
    rows, present, led, cov = _common(R, spec, members, sample_dt, ledger_dt, rules)
    members = present          # the rules below read each member's sheet row; members without a row are already reported (sample_not_in_table)
    kind = spec['kind']
    if kind in ('family', 'sex_set'):
        _sex_and_species(R, spec, members, rows, present)
        (_family if kind == 'family' else _sex_set)(R, spec, rules, members, rows, present, cov)
    else:
        _bc_parents(R, spec, rules, members, rows, present, cov)
    return R
