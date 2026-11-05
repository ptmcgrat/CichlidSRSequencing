"""Genome registry: one folder (<CichlidSequencingData>/YH_SD_Paper/genomes/<key>/) with every genome used in the project, each with a manifest
(accession, source URL, md5 of the download and of the FASTA, contig stats). Entries are never overwritten once their manifest exists.
Download from NCBI uses E-utilities (assembly name -> accession -> FTP path) and md5checksums.txt; the heavy lifting is curl/gunzip/samtools."""
import gzip, hashlib, json, os, re, shutil, subprocess, urllib.parse, urllib.request, datetime
import yaml
from .common_light import sh, config_path

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'


def load_registry(path=None):
    return yaml.safe_load(open(path or config_path('genomes.yaml')))['genomes']


# ------------------------------------------------------------------ pure helpers (unit-tested)
def guess_accession(assembly_name):
    """GenBank assemblies named ASM<digits>v<k> usually have accession GCA_0<digits>1.<k>. A GUESS: always verified against NCBI before use."""
    m = re.fullmatch(r'ASM(\d+)v(\d+)', assembly_name or '')
    return f'GCA_0{m.group(1)}1.{m.group(2)}' if m else None


def parse_esearch(payload):
    return list(payload.get('esearchresult', {}).get('idlist', []))


def parse_esummary(payload):
    """Returns a list of dicts (accession, assembly_name, organism, ftp_path) from an esummary JSON payload (db=assembly)."""
    out = []
    res = payload.get('result', {})
    for uid in res.get('uids', []):
        r = res[uid]
        ftp = r.get('ftppath_genbank') or r.get('ftppath_refseq') or ''
        out.append(dict(uid=uid, accession=r.get('assemblyaccession'), assembly_name=r.get('assemblyname'),
                        organism=r.get('speciesname') or r.get('organism'), ftp_path=ftp, status=r.get('assemblystatus')))
    return out


def ftp_to_https(ftp_path):
    return ftp_path.replace('ftp://', 'https://', 1)


def download_urls(ftp_path):
    """ftp_path like .../GCA_047663971.1_ASM4766397v1 -> (genomic.fna.gz URL, md5checksums.txt URL, file name)."""
    base = ftp_to_https(ftp_path).rstrip('/')
    stem = base.rsplit('/', 1)[-1]
    fn = f'{stem}_genomic.fna.gz'
    return f'{base}/{fn}', f'{base}/md5checksums.txt', fn


def parse_md5checksums(text):
    d = {}
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 2:
            d[parts[1].lstrip('./')] = parts[0]
    return d


def md5_of(path, chunk=1 << 22):
    h = hashlib.md5()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(chunk), b''):
            h.update(b)
    return h.hexdigest()


def fasta_stats(path):
    """(n_contigs, total_bp, longest) of a plain FASTA, streaming."""
    n = tot = cur = longest = 0
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                longest = max(longest, cur); cur = 0; n += 1
            else:
                cur += len(line.strip()); tot += len(line.strip())
    longest = max(longest, cur)
    return n, tot, longest


# ------------------------------------------------------------------ NCBI
def _get_json(url):
    with urllib.request.urlopen(url, timeout=60) as r:
        return json.loads(r.read().decode())


def resolve_ncbi(entry):
    """Find the assembly record by accession (if given) or assembly name; cross-check name and accession. Returns the esummary dict."""
    terms = []
    if entry.get('accession'):
        terms.append(f"{entry['accession']}[Assembly Accession]")
    if entry.get('assembly_name'):
        terms.append(f"{entry['assembly_name']}[Assembly Name]")
    g = guess_accession(entry.get('assembly_name'))
    if g and not entry.get('accession'):
        terms.append(f'{g}[Assembly Accession]')
    last_err = None
    for term in terms:
        try:
            ids = parse_esearch(_get_json(f"{EUTILS}esearch.fcgi?db=assembly&retmode=json&term={urllib.parse.quote(term)}"))
            if not ids:
                continue
            recs = parse_esummary(_get_json(f"{EUTILS}esummary.fcgi?db=assembly&retmode=json&id={','.join(ids)}"))
            for r in recs:
                if not entry.get('assembly_name') or r['assembly_name'] == entry['assembly_name']:
                    return r
        except Exception as e:                      # noqa
            last_err = e
    raise RuntimeError(f"could not resolve {entry.get('assembly_name') or entry.get('accession')} at NCBI ({last_err}). "
                       'If it is not public yet, supply the FASTA with  source: local, path: <fasta>  in config/genomes.yaml')


def fetch_ncbi(entry, outdir):
    rec = resolve_ncbi(entry)
    gz_url, md5_url, fn = download_urls(rec['ftp_path'])
    os.makedirs(outdir, exist_ok=True)
    gz = f'{outdir}/{fn}'
    sh(['curl', '-L', '--fail', '-C', '-', '--retry', '5', '-o', gz, gz_url])
    md5_remote = parse_md5checksums(sh(['curl', '-L', '--fail', '-s', md5_url]).stdout).get(fn)
    md5_gz = md5_of(gz)
    if md5_remote and md5_gz != md5_remote:
        raise RuntimeError(f'md5 mismatch for {fn}: {md5_gz} != {md5_remote}')
    return rec, gz, gz_url, md5_gz, bool(md5_remote)


def gunzip(gz, out):
    with gzip.open(gz, 'rb') as fi, open(out, 'wb') as fo:
        shutil.copyfileobj(fi, fo, 1 << 22)


# ------------------------------------------------------------------ registry operations
def entry_dir(proj, key):
    return f'{proj.genomes_dir}{key}/'


def fasta_path(proj, key):
    """Absolute path of the FASTA of a registered genome (raises if it has not been fetched)."""
    mf = entry_dir(proj, key) + 'manifest.json'
    if not os.path.exists(mf):
        raise FileNotFoundError(f'genome {key} is not registered yet: python -m yhsd.cli genomes fetch --only {key}')
    return json.load(open(mf))['fasta']


def _manifest(proj, key, entry, fasta, **extra):
    n, tot, longest = fasta_stats(fasta)
    m = dict(key=key, role=entry.get('role'), source=entry.get('source'), accession=entry.get('accession'), assembly_name=entry.get('assembly_name'),
             species=entry.get('species'), note=entry.get('note'), fasta=fasta, md5_fasta=md5_of(fasta), n_contigs=n, total_bp=tot, longest_contig=longest,
             registered=datetime.datetime.now().isoformat(timespec='seconds'))
    m.update(extra)
    json.dump(m, open(entry_dir(proj, key) + 'manifest.json', 'w'), indent=1)
    return m


def copy_verified(src, dst, expect_md5=None):
    """Real copy of src to dst (replacing a symlink), md5-verified against the source; an identical existing copy is kept. Returns the md5."""
    md5_src = expect_md5 or md5_of(src)
    if os.path.islink(dst):
        os.unlink(dst)
    if not (os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(src) and md5_of(dst) == md5_src):
        tmp = dst + '.part'
        shutil.copyfile(src, tmp)
        if md5_of(tmp) != md5_src:
            os.remove(tmp); raise RuntimeError(f'copy verification failed for {src}')
        os.replace(tmp, dst)
    fai = src + '.fai'
    if os.path.exists(fai):
        shutil.copyfile(fai, dst + '.fai')                      # identical sequence -> identical index
    else:
        sh(['samtools', 'faidx', dst])
    return md5_src


def register_filemanager(proj, key, entry):
    """Genome managed by the lab FileManager (on Dropbox): either a genome folder (genome_version + filemanager_file) or any file under the master directory
    (master_path). It is COPIED into genomes/<key>/ (the project folder is the archive), md5-verified against the source, the md5 in config (`md5`) and the md5sum.txt next to it."""
    d = entry_dir(proj, key); os.makedirs(d, exist_ok=True)
    if entry.get('master_path'):
        src = proj.fm.localMasterDir + entry['master_path']; fetch_target = src
    else:
        gdir = f"{proj.fm.localGenomesDir}{entry['genome_version']}/"; src = gdir + entry['filemanager_file']; fetch_target = gdir
    if not os.path.exists(src):
        print(f'[genomes] {key}: not on this machine, downloading with the FileManager ...')
        proj.fm.downloadData(fetch_target)
    if not os.path.exists(src):
        raise FileNotFoundError(f'{src} not found (FileManager download failed?)')
    verified = verify_against_md5sum_file(src)
    md5_src = md5_of(src)
    if entry.get('md5'):
        if md5_src != entry['md5']:
            raise RuntimeError(f"md5 mismatch for {src}: {md5_src} != expected {entry['md5']}")
        verified = True
    dst = d + os.path.basename(src)
    copy_verified(src, dst, md5_src)
    return _manifest(proj, key, entry, dst, copied_from=src, md5_verified=verified)


def materialize(proj, key):
    """Turn a genome registered by an earlier version as a symlink into a real copy inside the project folder."""
    mf = entry_dir(proj, key) + 'manifest.json'
    m = json.load(open(mf))
    src = m.get('linked_from')
    if not src:
        return m
    dst = entry_dir(proj, key) + os.path.basename(src)
    print(f'[genomes] {key}: converting the link to a real copy ({os.path.getsize(src) / 1e9:.2f} GB) ...')
    md5 = copy_verified(src, dst, m.get('md5_fasta'))
    m.update(fasta=dst, copied_from=src, md5_fasta=md5); m.pop('linked_from', None)
    json.dump(m, open(mf, 'w'), indent=1)
    return m


def verify_against_md5sum_file(path):
    """If an md5sum.txt (NCBI dataset package or `md5sum` output) sits next to `path`, check the file against it. Returns True / False / None (no listing)."""
    d = os.path.dirname(os.path.abspath(path)); mf = os.path.join(d, 'md5sum.txt')
    if not os.path.exists(mf):
        return None
    want = parse_md5checksums(open(mf).read()).get(os.path.basename(path))
    if want is None:
        return None
    got = md5_of(path)
    if got != want:
        raise RuntimeError(f'md5 mismatch for {path}: {got} != {want} (md5sum.txt)')
    return True


def register_local(proj, key, entry):
    d = entry_dir(proj, key); os.makedirs(d, exist_ok=True)
    src = os.path.expanduser(entry['path'])
    if not os.path.exists(src):
        raise FileNotFoundError(f'{src} not found')
    verified = verify_against_md5sum_file(src)
    dst = d + os.path.basename(src)
    if not os.path.exists(dst):
        shutil.copy2(src, dst)
    if dst.endswith('.gz'):
        out = dst[:-3]; gunzip(dst, out); dst = out
    sh(['samtools', 'faidx', dst])
    return _manifest(proj, key, entry, dst, copied_from=src, md5_verified_against_md5sum_txt=verified)


def register_ncbi(proj, key, entry):
    d = entry_dir(proj, key)
    rec, gz, url, md5_gz, md5_checked = fetch_ncbi(entry, d)
    fasta = gz[:-3]
    if not os.path.exists(fasta):
        gunzip(gz, fasta)
    sh(['samtools', 'faidx', fasta])
    entry = dict(entry, accession=rec['accession'], assembly_name=rec['assembly_name'])
    return _manifest(proj, key, entry, fasta, ncbi_url=url, md5_download=md5_gz, md5_verified_against_ncbi=md5_checked, organism_ncbi=rec['organism'])


# ------------------------------------------------------------------ repeat annotation (RepeatMasker) registered next to the genome
def read_fai(path):
    return {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(path)}


def annotation_names(path):
    """Counter-like dict of sequence names used in a RepeatMasker .out (column 5) or GFF (column 1)."""
    names = {}
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as f:
        for line in f:
            if line.startswith('#'):
                continue
            if path.replace('.gz', '').endswith(('.gff', '.gff3')):
                sp = line.split('\t', 1)
                n = sp[0]
            else:
                sp = line.split()
                if len(sp) < 11 or not sp[0].isdigit():
                    continue
                n = sp[4]
            names[n] = names.get(n, 0) + 1
    return names


def fasta_name_lengths(path):
    """dict name -> length of a (possibly gzipped) FASTA, streaming."""
    out, name, n = {}, None, 0
    with (gzip.open(path, 'rt') if path.endswith('.gz') else open(path)) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None: out[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name is not None: out[name] = n
    return out


_CHROM_LIKE = re.compile(r'(?:chr|Chr|CHR|LG|lg)?_?0*(\d+)$')


def validate_annotation(names, fai_lengths, lg_to_chrom=None, ann_lengths=None, big=5_000_000, min_fraction=0.95):
    """Map annotation sequence names to genome sequences and decide whether the annotation is usable for chromosome-scale work.
    1) identity;  2) LG label (LG5 / chr5 / 5 -> the genome's LG5);  3) a UNIQUE length match (needs ann_lengths, e.g. from the .masked FASTA), which maps scaffold01 -> NW_...
    Usable = every chromosome-scale name (length >= big, or LG-like) is mapped, no mapped sequence has a different length, and at least `min_fraction` of all annotation
    records land on a mapped sequence. Unmapped small names (scaffolds, a mitochondrion absent from the genome) are reported but do not block use."""
    mapping, used = {}, set()
    for n in names:
        if n in fai_lengths:
            mapping[n] = n; used.add(n)
    for n in names:
        if n in mapping: continue
        m = _CHROM_LIKE.fullmatch(n)
        if m and lg_to_chrom and f'LG{int(m.group(1))}' in lg_to_chrom and lg_to_chrom[f'LG{int(m.group(1))}'] not in used:
            mapping[n] = lg_to_chrom[f'LG{int(m.group(1))}']; used.add(mapping[n])
    if ann_lengths:
        by_len = {}
        for s, L in fai_lengths.items():
            by_len.setdefault(L, []).append(s)
        for n in names:
            if n in mapping: continue
            cands = [s for s in by_len.get(ann_lengths.get(n), []) if s not in used]
            if len(cands) == 1:
                mapping[n] = cands[0]; used.add(cands[0])
    unmapped = [n for n in names if n not in mapping]
    mismatches = [n for n, s in mapping.items() if ann_lengths and n in ann_lengths and ann_lengths[n] != fai_lengths[s]]
    chromlike = [n for n in names if (ann_lengths or {}).get(n, 0) >= big or fai_lengths.get(n, 0) >= big or (_CHROM_LIKE.fullmatch(n) and not n.lower().startswith('scaffold'))]
    total = sum(names.values()); mapped = sum(names[n] for n in mapping)
    usable = bool(chromlike) and all(n in mapping for n in chromlike) and not mismatches and (mapped / total if total else 0.0) >= min_fraction
    return dict(n_sequence_names=len(names), n_mapped=len(mapping), mapped_record_fraction=mapped / total if total else 0.0, n_unmapped=len(unmapped),
                unmapped_names=unmapped[:20], length_mismatches=mismatches[:20], n_chromosome_scale=len(chromlike), usable=usable,
                identity=all(k == v for k, v in mapping.items()) and not unmapped, mapping=mapping)


def parse_tbl(path):
    """Provenance lines of a RepeatMasker .tbl: versions, library name, sequence count and total length."""
    txt = open(path).read()
    g = lambda pat: (re.search(pat, txt) or [None, None])[1]
    lib = g(r'compared to classified sequences in "([^"]*)"') or g(r'compared to [^"\n]*"([^"]*)"')
    return dict(repeatmasker_version=g(r'RepeatMasker version\s+([^\s,]+)'), rmblastn_version=g(r'rmblastn version\s+([^\s,]+)'), library=lib,
                library_name=os.path.basename(lib) if lib else None, n_sequences=int(g(r'sequences:\s+(\d+)') or 0), total_length=int((g(r'total length:\s+(\d+)') or '0')))


def register_annotation(proj, key, spec, upload=True):
    """Copy the annotation files of genome `key` (RepeatMasker .out/.gff/.tbl/.cat.gz/.masked) into genomes/<key>/annotation/, md5 them, map the sequence names to the
    registered FASTA (identity / LG label / unique length from `lengths_from`), record the RepeatMasker version and library from the .tbl and compare them with the
    reference genome's annotation, write a manifest + name map, and upload to Dropbox. Sources are never modified."""
    from .runner import sync
    mf = entry_dir(proj, key) + 'manifest.json'
    if not os.path.exists(mf):
        raise FileNotFoundError(f'register the genome {key} first (genomes fetch --only {key})')
    fasta = json.load(open(mf))['fasta']
    dst = entry_dir(proj, key) + 'annotation/'; os.makedirs(dst, exist_ok=True)
    src_dir = os.path.expanduser(spec.get('source_dir') or dst)
    files = {}
    for fn in spec['files']:
        s, d = os.path.join(src_dir, fn), dst + fn
        if not os.path.exists(s):
            raise FileNotFoundError(f'{s} not found')
        if os.path.abspath(s) != os.path.abspath(d) and not (os.path.exists(d) and os.path.getsize(d) == os.path.getsize(s)):
            shutil.copy2(s, d)
        files[fn] = dict(size=os.path.getsize(d), md5=md5_of(d), source=s)
    main = dst + spec['main']
    fai = fasta + '.fai'
    if not os.path.exists(fai):
        sh(['samtools', 'faidx', fasta])
    lg_to_chrom = {lg: c for c, lg in proj.lg.items()} if key == proj.genome else None
    names = annotation_names(main)
    ann_lengths = fasta_name_lengths(dst + spec['lengths_from']) if spec.get('lengths_from') else None
    v = validate_annotation(names, read_fai(fai), lg_to_chrom, ann_lengths)
    with open(dst + 'name_map.tsv', 'w') as f:
        f.write('annotation_name\tgenome_sequence\n' + ''.join(f'{a}\t{b}\n' for a, b in sorted(v['mapping'].items())))
    tbl = next((fn for fn in spec['files'] if fn.endswith('.tbl')), None)
    prov = parse_tbl(dst + tbl) if tbl else {}
    warnings = []
    ref = spec.get('compare_to')
    if ref and prov:
        rm = entry_dir(proj, ref) + 'annotation/annotation_manifest.json'
        if os.path.exists(rm):
            rp = json.load(open(rm)).get('provenance', {})
            for k in ('repeatmasker_version', 'rmblastn_version', 'library_name'):
                if rp.get(k) and prov.get(k) != rp.get(k):
                    warnings.append(f'{k}: {prov.get(k)} here vs {rp.get(k)} in {ref}')
    man = dict(key=key, main=spec['main'], files=files, validation={k: v[k] for k in v if k != 'mapping'}, usable=v['usable'], provenance=prov, comparability_warnings=warnings,
               registered=datetime.datetime.now().isoformat(timespec='seconds'))
    json.dump(man, open(dst + 'annotation_manifest.json', 'w'), indent=1)
    print(f"[annotation] {key}: {v['n_sequence_names']} sequence names, {v['n_mapped']} mapped ({v['mapped_record_fraction']:.2%} of records), {v['n_chromosome_scale']} chromosome-scale, "
          f"unmapped: {v['unmapped_names'][:5]}, length mismatches: {v['length_mismatches'][:5]}" + ('' if v['usable'] else '  <-- NOT USABLE until the chromosome names/lengths are reconciled'))
    if prov: print(f"[annotation] {key}: RepeatMasker {prov.get('repeatmasker_version')}, rmblastn {prov.get('rmblastn_version')}, library {prov.get('library_name')}")
    for w in warnings:
        print(f'[annotation] WARNING {key} differs from {ref}: {w}')
    if upload and proj.cfg.get('runner', {}).get('upload', True):
        sync(proj, dst)
    return man


def excluded_regions(proj):
    """{reference sequence name: [(start, end), ...]} from `known_misassemblies` of the project's reference genome in config/genomes.yaml (LG labels resolved through proj.chrom_for_lg)."""
    entry = load_registry().get(proj.genome, {})
    out = {}
    for m in entry.get('known_misassemblies') or []:
        out.setdefault(proj.chrom_for_lg(m['lg']), []).append((int(m['start']), int(m['end'])))
    return out


def annotation_for(proj, key):
    """(path of the main annotation file, {annotation name: genome sequence}) of a registered genome; raises FileNotFoundError if none / unusable."""
    d = entry_dir(proj, key) + 'annotation/'
    mf = d + 'annotation_manifest.json'
    if not os.path.exists(mf):
        raise FileNotFoundError(f'no registered annotation for {key}: run `genomes annotate --only {key}`')
    man = json.load(open(mf))
    if not man.get('usable'):
        raise FileNotFoundError(f'annotation of {key} is registered but its sequence names do not match the genome (see {mf})')
    nm = dict(l.rstrip('\n').split('\t') for l in open(d + 'name_map.tsv').readlines()[1:])
    return d + man['main'], nm


def fetch(proj, only=None, upload=True):
    """Register / download every genome (or those in `only`). One failing entry does not stop the others; entries still needing input are listed at the end."""
    from .runner import sync
    reg = load_registry()
    rows, problems = [], []
    for key, entry in reg.items():
        if only and key not in only:
            continue
        mf = entry_dir(proj, key) + 'manifest.json'
        if os.path.exists(mf):
            m = json.load(open(mf))
            if m.get('linked_from'):
                try: m = materialize(proj, key)
                except Exception as e: problems.append((key, f'could not convert the link to a copy: {str(e)[:200]}'))
            print(f'[genomes] {key}: already registered ({m["fasta"]})'); rows.append(m); continue
        if entry['source'] == 'local' and not entry.get('path'):
            problems.append((key, 'no path set: fill `path:` in config/genomes.yaml')); print(f'[genomes] {key}: skipped (no path set)'); continue
        try:
            os.makedirs(entry_dir(proj, key), exist_ok=True)
            fn = dict(filemanager=register_filemanager, ncbi=register_ncbi, local=register_local)[entry['source']]
            m = fn(proj, key, entry)
            print(f"[genomes] {key}: {m['n_contigs']} contigs, {m['total_bp'] / 1e6:.1f} Mb, md5 {m['md5_fasta']}")
            rows.append(m)
        except Exception as e:                      # noqa
            problems.append((key, str(e)[:300])); print(f'[genomes] {key}: FAILED: {str(e)[:300]}')
    write_table(proj)
    if upload and proj.cfg.get('runner', {}).get('upload', True):
        for m in rows:
            sync(proj, entry_dir(proj, m['key']))
        sync(proj, proj.genomes_dir + 'genomes.tsv')
    if problems:
        print('\nNeeds attention:'); [print(f'  {k}: {why}') for k, why in problems]
    return rows


def write_table(proj):
    import pandas as pd
    rows = []
    for key in load_registry():
        mf = entry_dir(proj, key) + 'manifest.json'
        if os.path.exists(mf):
            m = json.load(open(mf)); rows.append({k: m.get(k) for k in ('key', 'role', 'accession', 'assembly_name', 'species', 'n_contigs', 'total_bp', 'md5_fasta', 'fasta', 'registered')})
    pd.DataFrame(rows).to_csv(proj.genomes_dir + 'genomes.tsv', sep='\t', index=False)
    return rows
