"""
make_intervals.py  (Oct 2026; replaces the per-genome make_<genome>_intervals.py scripts)

Builds the GATK interval lists that callVariants.py uses for GenomicsDBImport / GenotypeGVCFs:

    <folder containing this script>/Intervals/<genome_version>/1.interval_list ... 96.interval_list

<genome_version> is the same name callVariants.py and nikesh_file_manager.py use (e.g. Mzebra_GT3_NCBI,
Mconophoros_GT1_NCBI). The reference fasta is located and downloaded through nikesh_file_manager.FileManager.
This script only downloads from Dropbox; it never uploads anything.

What it does
  1. Checks the genome name: it must be configured in nikesh_file_manager.py (_createMasterDirs) and its fasta must
     be on Dropbox. Misspelled names get "did you mean ..." suggestions.
  2. Downloads the fasta (and its .fai / .dict if they're on Dropbox; rclone skips files that are already up to date).
     A missing .fai / .dict is created locally next to the fasta, since GATK needs both anyway.
  3. Sorts the sequences by length (from pyfaidx), using the RefSeq naming convention (NC_ = chromosome or mitogenome,
     NW_ = unplaced scaffold) only as a cross-check:
       - chromosome: >= --min_chrom_length (default 30 Mb; cichlid LGs are ~30-85 Mb)       -> KEPT
       - mitogenome: 15-20 kb (cichlid mtDNA is ~16.6 kb) and not an NW_/NT_ scaffold        -> EXCLUDED
       - everything else (unplaced scaffolds/contigs, e.g. NW_*, *_unaligned)               -> EXCLUDED
     A name that disagrees with its length (e.g. an NC_ sequence below the cut-off) is printed as a WARNING.
  4. Splits the kept chromosomes into exactly --num_intervals (default 96) pieces of near-equal size. Each extra
     interval goes to whichever chromosome currently has the largest pieces, which keeps the largest interval (= the
     slowest GenomicsDBImport / GenotypeGVCFs job) as small as possible. Each file holds one interval and the full
     sequence dictionary header from the genome's .dict.
  5. Writes interval_manifest.tsv (which file covers which region) and excluded_contigs.tsv (what was left out and
     why) next to the interval files, then re-reads everything to verify full, non-overlapping coverage.

Usage (run from anywhere; output always goes next to this script unless --output_root is given):
    python3 make_intervals.py Mconophoros_GT1_NCBI
"""
import argparse, os, sys, shutil, subprocess, difflib, datetime, heapq
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser(usage = 'python3 make_intervals.py [genome_version] [options]',
                                 description = 'Generate Intervals/<genome_version>/1..N.interval_list for callVariants.py')
parser.add_argument('genome_version', nargs = '?', help = 'genome name as used by callVariants.py / nikesh_file_manager.py (e.g. Mconophoros_GT1_NCBI). Asked for if left off.')
parser.add_argument('--num_intervals', type = int, default = 96, help = 'number of interval files to make (default 96)')
parser.add_argument('--output_root', default = SCRIPT_DIR + '/Intervals/', help = 'folder that holds one sub-folder per genome (default: Intervals/ next to this script)')
parser.add_argument('--overwrite', action = 'store_true', help = 'replace existing interval files for this genome')
parser.add_argument('--min_chrom_length', type = int, default = 30000000, help = 'sequences at least this long are treated as chromosomes (default 30,000,000 bp; cichlid LGs are ~30-85 Mb, unplaced scaffolds <= ~1 Mb)')
args = parser.parse_args()

# How sequences are classified (cichlid genomes). Lengths come from pyfaidx; names only confirm or flag.
MITO_LENGTH_RANGE = (15000, 20000)   # vertebrate/cichlid mitogenomes are ~16-17 kb (M. zebra NC_027944.1 = 16,582 bp)
CHROMOSOME_PREFIX = 'NC_'            # RefSeq accession for chromosomes - and for the mitogenome
UNPLACED_PREFIXES = ('NW_', 'NT_')   # RefSeq accessions for unplaced / unlocalized scaffolds


def resolve_genome(genome_version):
    # Make sure the name is one nikesh_file_manager.py knows and that its fasta is on Dropbox
    if not genome_version:
        genome_version = input('Reference genome name (as used by callVariants.py, e.g. Mconophoros_GT1_NCBI): ').strip()
    if not genome_version:
        sys.exit('No genome name given.')

    fm = FM(genome_version)
    cloud_genomes = fm.returnGenomeVersions()

    if not hasattr(fm, 'localGenomeFile'):
        case_match = [g for g in cloud_genomes if g.lower() == genome_version.lower()]
        suggestions = case_match or difflib.get_close_matches(genome_version, cloud_genomes, n = 3, cutoff = 0.5)
        message = f"'{genome_version}' has no reference fasta configured in nikesh_file_manager.py (_createMasterDirs)."
        if suggestions:
            message += f" Did you mean: {', '.join(suggestions)}?"
        message += f"\nGenome folders on Dropbox: {', '.join(cloud_genomes) if cloud_genomes else '(none found - check rclone)'}"
        message += '\nIf this is a new genome, add an elif branch for it in nikesh_file_manager._createMasterDirs first.'
        sys.exit(message)

    if genome_version not in cloud_genomes:
        sys.exit(f"No '{genome_version}' folder in the Genomes directory on Dropbox. Genome folders found: {', '.join(cloud_genomes)}")
    if not fm.checkCloudFile(fm.localGenomeFile):
        sys.exit(f"The fasta configured in nikesh_file_manager.py ({os.path.basename(fm.localGenomeFile)}) isn't in the Dropbox folder "
                 f"Genomes/{genome_version}/. Files there: {', '.join(fm.returnCloudFiles(fm.localGenomeDir)) or '(folder is empty)'}")
    return fm


def fetch_reference_files(fm):
    # Download the fasta (+ .fai / .dict when they're on Dropbox) through the FileManager; rclone skips files that are
    # already up to date locally. A missing .fai / .dict is created locally only (nothing is uploaded).
    fasta = fm.localGenomeFile
    fai = fasta + '.fai'
    dict_file = os.path.splitext(fasta)[0] + '.dict' # GATK's naming: reference.fna -> reference.dict

    print(f"Getting {os.path.basename(fasta)} from Dropbox (skipped by rclone if the local copy is current) ...")
    fm.downloadData(fasta)
    for index_file in [fai, dict_file]:
        if fm.checkCloudFile(index_file):
            fm.downloadData(index_file)

    if not os.path.exists(fai):
        print(f"Creating {os.path.basename(fai)} locally ...")
        Fasta(fasta) # pyfaidx writes the .fai when it opens a fasta without one
    if not os.path.exists(dict_file):
        print(f"Creating {os.path.basename(dict_file)} locally ...")
        if shutil.which('gatk'):
            output = subprocess.run(['gatk', 'CreateSequenceDictionary', '-R', fasta, '-O', dict_file], capture_output = True, encoding = 'utf-8')
        elif shutil.which('samtools'):
            output = subprocess.run(['samtools', 'dict', '-o', dict_file, fasta], capture_output = True, encoding = 'utf-8')
        else:
            sys.exit(f"{dict_file} doesn't exist and neither gatk nor samtools is on the PATH to create it.")
        if output.returncode != 0 or not os.path.exists(dict_file):
            sys.exit(f"Failed to create {dict_file}:\n{output.stderr}")
    return fasta, fai, dict_file


def read_dict(dict_file):
    # Header lines for the interval lists (@HD + every @SQ) and the (name, length) of each sequence, in dictionary order
    header, sequences = [], []
    with open(dict_file) as fh:
        for line in fh:
            if not line.startswith('@'):
                continue
            header.append(line.rstrip('\n'))
            if line.startswith('@SQ'):
                fields = dict(field.split(':', 1) for field in line.rstrip('\n').split('\t')[1:] if ':' in field)
                sequences.append((fields['SN'], int(fields['LN'])))
    if not sequences:
        sys.exit(f"No @SQ lines found in {dict_file}")
    return header, sequences


def read_fasta_lengths(fasta):
    # (name, length) of every sequence in the fasta, in file order, via pyfaidx (uses the .fai index; creates it if missing)
    index = Fasta(fasta)
    return [(name, len(index[name])) for name in index.keys()]


def classify_sequence(name, length, min_chrom_length):
    # Returns (category, reason, warning). category is 'chromosome', 'mito' or 'unplaced'; warning is None unless the
    # NC_/NW_ name disagrees with what the length says. Length decides.
    is_unplaced_name = name.startswith(UNPLACED_PREFIXES)
    if length >= min_chrom_length:
        warning = f"{name} is named like an unplaced scaffold (NW_/NT_) but is {length:,} bp - kept as a chromosome because of its length" if is_unplaced_name else None
        return 'chromosome', f"{length:,} bp >= {min_chrom_length:,} bp", warning
    if MITO_LENGTH_RANGE[0] <= length <= MITO_LENGTH_RANGE[1] and not is_unplaced_name:
        return 'mito', f"{length:,} bp is mitogenome-sized ({MITO_LENGTH_RANGE[0]:,}-{MITO_LENGTH_RANGE[1]:,} bp) and not an NW_/NT_ scaffold", None
    warning = f"{name} has a chromosome-style (NC_) name but is only {length:,} bp - excluded because it's below --min_chrom_length" if name.startswith(CHROMOSOME_PREFIX) else None
    return 'unplaced', f"{length:,} bp < {min_chrom_length:,} bp", warning


def choose_chromosomes(sequences):
    # Returns (kept [(name, length)], excluded [(name, length, category, reason)]) and prints a summary + checks
    kept, excluded, warnings = [], [], []
    for name, length in sequences:
        category, reason, warning = classify_sequence(name, length, args.min_chrom_length)
        if category == 'chromosome':
            kept.append((name, length))
        else:
            excluded.append((name, length, category, reason))
        if warning:
            warnings.append(warning)

    if not kept:
        sys.exit(f"No sequences >= {args.min_chrom_length:,} bp found. Lower --min_chrom_length.")

    # Summary
    print(f"\nKEEPING {len(kept)} chromosomes ({sum(l for _, l in kept):,} bp):")
    for name, length in kept:
        print(f"  {name:<28}{length:>14,} bp")
    mito = [x for x in excluded if x[2] == 'mito']
    unplaced = sorted((x for x in excluded if x[2] == 'unplaced'), key = lambda x: -x[1])
    print(f"\nEXCLUDING {len(excluded)} sequences ({sum(x[1] for x in excluded):,} bp):")
    print(f"  mitogenome: " + (', '.join(f"{x[0]} ({x[1]:,} bp)" for x in mito) if mito else 'none found'))
    if unplaced:
        prefixes = {}
        for x in unplaced:
            prefix = x[0][:3] if x[0][:3] in ('NW_', 'NT_', 'NC_') else 'other names'
            prefixes[prefix] = prefixes.get(prefix, 0) + 1
        print(f"  unplaced/other: {len(unplaced)} sequences ({', '.join(f'{v} {k}' for k, v in prefixes.items())}), "
              f"{sum(x[1] for x in unplaced):,} bp. Largest: " + ', '.join(f"{x[0]} ({x[1]:,} bp)" for x in unplaced[:3]))

    # Checks
    smallest_kept = min(kept, key = lambda x: x[1])
    print(f"\nLength gap: smallest kept chromosome {smallest_kept[0]} = {smallest_kept[1]:,} bp; "
          + (f"largest unplaced {unplaced[0][0]} = {unplaced[0][1]:,} bp" if unplaced else 'no unplaced sequences'))
    if len(mito) == 0:
        print(f"Note: no {MITO_LENGTH_RANGE[0]:,}-{MITO_LENGTH_RANGE[1]:,} bp mitogenome-sized sequence found (fine if this assembly doesn't include the mitogenome).")
    elif len(mito) > 1:
        warnings.append(f"{len(mito)} mitogenome-sized sequences found (all excluded): {', '.join(x[0] for x in mito)}")
    for warning in warnings:
        print('WARNING: ' + warning)
    return kept, excluded


def allocate_intervals(kept, num_intervals):
    # Every chromosome gets at least 1 interval. Each remaining interval goes to the chromosome whose pieces are currently
    # the largest, so the biggest interval ends up as small as possible.
    if num_intervals < len(kept):
        sys.exit(f"--num_intervals {num_intervals} is less than the {len(kept)} chromosomes kept (each needs at least one interval).")
    allocation = {name: 1 for name, _ in kept}
    lengths = dict(kept)
    heap = [(-length, name) for name, length in kept]
    heapq.heapify(heap)
    for _ in range(num_intervals - len(kept)):
        _, name = heapq.heappop(heap)
        allocation[name] += 1
        heapq.heappush(heap, (-lengths[name] / allocation[name], name))
    return allocation


def split_contig(length, pieces):
    # 1-based, inclusive (start, end) pairs covering 1..length in `pieces` near-equal parts (sizes differ by at most 1 bp)
    base, extra = divmod(length, pieces)
    start = 1
    for i in range(pieces):
        end = start + base + (1 if i < extra else 0) - 1
        yield start, end
        start = end + 1


def write_intervals(genome_version, header, kept, excluded, allocation, fasta, dict_file):
    output_dir = os.path.join(args.output_root, genome_version) + '/'
    existing = [f for f in os.listdir(output_dir) if f.endswith('.interval_list')] if os.path.isdir(output_dir) else []
    if existing and not args.overwrite:
        sys.exit(f"{output_dir} already has {len(existing)} interval files. Use --overwrite to replace them.")

    # Write into a temporary folder first, then swap it in, so a failed run never leaves a half-written interval set
    temp_dir = output_dir.rstrip('/') + f".tmp_{os.getpid()}/"
    os.makedirs(temp_dir)
    manifest = []
    file_num = 1
    for name, length in kept:
        for start, end in split_contig(length, allocation[name]):
            with open(f"{temp_dir}{file_num}.interval_list", 'w') as fh:
                fh.write('\n'.join(header) + '\n')
                fh.write(f"{name}\t{start}\t{end}\t+\t.\n")
            manifest.append((file_num, name, start, end, end - start + 1))
            file_num += 1

    stamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M')
    with open(temp_dir + 'interval_manifest.tsv', 'w') as fh:
        fh.write(f"# make_intervals.py {stamp} | genome {genome_version} | fasta {fasta} | dict {dict_file} | {len(manifest)} intervals\n")
        fh.write('file\tcontig\tstart\tend\tlength_bp\n')
        for file_num, name, start, end, size in manifest:
            fh.write(f"{file_num}.interval_list\t{name}\t{start}\t{end}\t{size}\n")
    with open(temp_dir + 'excluded_contigs.tsv', 'w') as fh:
        fh.write('contig\tlength_bp\tcategory\treason\n')
        for name, length, category, reason in excluded:
            fh.write(f"{name}\t{length}\t{category}\t{reason}\n")

    if os.path.isdir(output_dir):
        shutil.rmtree(output_dir)
    os.rename(temp_dir, output_dir)
    return output_dir, manifest


def verify_intervals(output_dir, kept, header, num_intervals):
    # Re-read what was written: right number of files, correct header, one interval each, and every kept chromosome
    # covered from 1 to its end with no gaps or overlaps
    files = sorted((f for f in os.listdir(output_dir) if f.endswith('.interval_list')), key = lambda f: int(f.split('.')[0]))
    expected = [f"{i}.interval_list" for i in range(1, num_intervals + 1)]
    assert files == expected, f"expected files 1..{num_intervals}.interval_list, found {len(files)} files"
    covered = {}
    for f in files:
        with open(output_dir + f) as fh:
            lines = fh.read().splitlines()
        assert lines[:len(header)] == header, f"{f}: header doesn't match the .dict"
        body = [line for line in lines[len(header):] if line.strip()]
        assert len(body) == 1, f"{f}: expected 1 interval, found {len(body)}"
        name, start, end = body[0].split('\t')[:3]
        covered.setdefault(name, []).append((int(start), int(end)))
    for name, length in kept:
        position = 1
        for start, end in sorted(covered.pop(name)):
            assert start == position, f"{name}: gap or overlap at {position}"
            position = end + 1
        assert position == length + 1, f"{name}: coverage ends at {position - 1}, contig length is {length}"
    assert not covered, f"intervals on unexpected contigs: {list(covered)}"


if __name__ == '__main__':
    fm = resolve_genome(args.genome_version)
    genome_version = fm.genome_version
    fasta, fai, dict_file = fetch_reference_files(fm)

    # Names/lengths from the fasta (pyfaidx); the interval-list header comes from the .dict, so the two must agree
    sequences = read_fasta_lengths(fasta)
    header, dict_sequences = read_dict(dict_file)
    if sequences != dict_sequences:
        dict_lengths = dict(dict_sequences)
        mismatched = [name for name, length in sequences if dict_lengths.get(name) != length]
        sys.exit(f"{os.path.basename(dict_file)} doesn't match the fasta ({len(dict_sequences)} vs {len(sequences)} sequences; "
                 f"differing names/lengths e.g. {mismatched[:3]}). One of the index files is stale: delete "
                 f"{os.path.basename(fai)} and {os.path.basename(dict_file)} and re-run so they're recreated.")

    kept, excluded = choose_chromosomes(sequences)
    allocation = allocate_intervals(kept, args.num_intervals)

    print(f"\n{args.num_intervals} intervals: {sum(l for _, l in kept) // args.num_intervals:,} bp on average, "
          f"largest {max(-(-length // allocation[name]) for name, length in kept):,} bp, "
          f"smallest {min(length // allocation[name] for name, length in kept):,} bp")
    print('Intervals per chromosome: ' + ', '.join(f"{name}={allocation[name]}" for name, _ in kept))

    output_dir, manifest = write_intervals(genome_version, header, kept, excluded, allocation, fasta, dict_file)
    verify_intervals(output_dir, kept, header, args.num_intervals)
    print(f"\nWrote and verified {len(manifest)} interval files in {output_dir}")
    print(f"  interval_manifest.tsv: which region each file covers")
    print(f"  excluded_contigs.tsv: {len(excluded)} sequences left out and why")
