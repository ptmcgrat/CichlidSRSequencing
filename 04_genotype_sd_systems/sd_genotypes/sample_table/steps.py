"""sample_table analysis. Step `inversions`: add Blumer et al. 2025 inversion genotypes (chr2, 9, 11, 13; chr10 is validated) to the
SampleDatabase. Inputs (Table_S1.tsv, Table_S7.tsv, Data_S1.nwk, Data_S2.nwk) are copied into the run folder, so the exact external
data used are stored on Dropbox with the result. The Google sheet is changed ONLY with --write-sheet, and only by appending new columns."""
import os, shutil
import numpy as np
import pandas as pd
from .inversions import annotate, NEW_COLS
from . import bc_inversions as bci
import yaml

EXTERNAL = ['Table_S1.tsv', 'Table_S7.tsv', 'Data_S1.nwk', 'Data_S2.nwk']
SHEET_COLS = ['BlumerSequenceID', 'BlumerSpecimenID'] + NEW_COLS + ['InversionSource']


def col_letter(n):                                    # 1-based -> A, B, ..., AA
    s = ''
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def write_new_columns(proj, table, run_dir, execute=False):
    """Append SHEET_COLS to the SampleDatabase sheet, matched by SampleID. Never touches existing columns. Dry run unless execute."""
    ws = proj.fm.gc.open_by_key(proj.fm.s_ID).worksheet('SampleDatabase')
    vals = ws.get_all_values()
    pd.DataFrame(vals).to_csv(run_dir + 'SampleDatabase_sheet_backup.tsv', sep='\t', index=False, header=False)     # backup first
    header = vals[0]
    present = [c for c in SHEET_COLS if c in header]
    if present:
        raise RuntimeError(f'columns already exist in the sheet: {present}; refusing to overwrite')
    first_blank = next((i for i, h in enumerate(header) if h == ''), len(header))
    # every cell in the target block must be empty
    for row in vals[1:]:
        if any(c != '' for c in row[first_blank:first_blank + len(SHEET_COLS)]):
            raise RuntimeError('target columns are not empty in at least one row; refusing to write')
    id_col = header.index('SampleID')
    lookup = table.set_index('SampleID')[SHEET_COLS]
    block = [SHEET_COLS]
    for row in vals[1:]:
        sid = row[id_col] if id_col < len(row) else ''
        if sid in lookup.index:
            r = lookup.loc[sid]
            block.append(['' if (isinstance(v, float) and np.isnan(v)) else (int(v) if isinstance(v, (float, np.floating)) else v) for v in r.tolist()])
        else:
            block.append([''] * len(SHEET_COLS))
    rng = f'{col_letter(first_blank + 1)}1:{col_letter(first_blank + len(SHEET_COLS))}{len(vals)}'
    print(f'[sheet] would write {len(SHEET_COLS)} columns x {len(vals) - 1} rows to {rng}')
    if execute:
        ws.update(range_name=rng, values=block, value_input_option='RAW')
        print('[sheet] written (backup saved in the run folder)')


LEGACY_MC_ALL = {'LG2': 'NC_036781.1', 'LG9': 'NC_036788.1', 'LG10': 'NC_036789.1', 'LG11': 'NC_036790.1', 'LG13': 'NC_036792.1'}
LEGACY_MC_CHROMS = {k: LEGACY_MC_ALL[k] for k in ('LG2', 'LG13')}       # Mconophoros_GT1 naming used by the first BC genotype file


def fill_blank_cells(proj, table, column, run_dir, execute=False):
    """Fill ONLY blank cells of an existing sheet column (e.g. Inversion10) from `table`; never overwrites a non-empty cell."""
    ws = proj.fm.gc.open_by_key(proj.fm.s_ID).worksheet('SampleDatabase')
    vals = ws.get_all_values()
    header = vals[0]
    if column not in header:
        raise RuntimeError(f'column {column} not found in the sheet')
    c, id_col = header.index(column), header.index('SampleID')
    lookup = table.set_index('SampleID')[column]
    data = []
    for k, row in enumerate(vals[1:], start=2):                      # sheet row number (row 1 is the header)
        sid = row[id_col] if id_col < len(row) else ''
        cur = row[c] if c < len(row) else ''
        if cur == '' and sid in lookup.index and pd.notna(lookup[sid]):
            data.append({'range': f'{col_letter(c + 1)}{k}', 'values': [[int(lookup[sid])]]})
    print(f'[sheet] {column}: {len(data)} blank cells would be filled')
    pd.DataFrame([(d['range'], d['values'][0][0]) for d in data], columns=['cell', 'value']).to_csv(run_dir + f'{column}_cells_to_fill.tsv', sep='\t', index=False)
    if execute and data:
        ws.batch_update(data, value_input_option='RAW')
        print(f'[sheet] {column}: filled {len(data)} cells')
    return len(data)


def resolve_bc_genotypes(proj, spec):
    """'legacy' -> <localNikeshDir>/OutputGenotypes.csv. If the file is missing locally and lies under the FileManager's local master directory,
    download it from Dropbox with the FileManager (the same rclone path the original GenotypeBC code used)."""
    path = proj.fm.localNikeshDir + 'OutputGenotypes.csv' if spec == 'legacy' else os.path.expanduser(spec)
    if not os.path.exists(path):
        if not path.startswith(proj.fm.localMasterDir):
            raise FileNotFoundError(f'{path} not found, and it is outside the FileManager master directory ({proj.fm.localMasterDir}), so it cannot be downloaded')
        print(f'[sample_table] {path} not on utaka; downloading with the FileManager ...')
        proj.ensure_local(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f'download via FileManager did not produce {path}')
    return path


def inversions_step(proj, blumer_dir, write_sheet=False, bc_genotypes=None, bc_legacy_mc=False):
    with proj.run('sample_table', 'inversions', params=dict(blumer_dir=blumer_dir, write_sheet=write_sheet, bc_genotypes=bc_genotypes, bc_legacy_mc=bc_legacy_mc)) as run:
        for f in EXTERNAL:                                         # permanent copy of the external inputs, with md5 in the manifest
            src = os.path.join(blumer_dir, f)
            run.add_input(src)                                       # copied into <run>/inputs/ by the runner (small files), md5 in the manifest
        proj.fm.readSampleDatabase()
        db = proj.fm.sample_dt.copy()
        s1 = pd.read_csv(os.path.join(blumer_dir, 'Table_S1.tsv'), sep='\t')
        s7 = pd.read_csv(os.path.join(blumer_dir, 'Table_S7.tsv'), sep='\t')
        table, rep = annotate(db, s1, s7)
        if bc_genotypes:                                         # infer chr2 / chr13 for YH x MC BC fish from the haplotype calls
            bcfg = proj.cfg['bc_inversions']
            bc_genotypes = resolve_bc_genotypes(proj, bc_genotypes)
            run.add_input(bc_genotypes)
            grid = bci.load_grid(bc_genotypes)
            chroms = LEGACY_MC_CHROMS if bc_legacy_mc else {lg: proj.chrom_for_lg(lg) for lg in bcfg['lgs']}
            calls, blocks = bci.infer(grid, chroms, bcfg['mc_inverted'], bcfg['max_rate'], bcfg['min_block_mb'] * 1e6, bcfg['min_consistency'])
            calls.to_csv(run.dir + 'bc_inversion_calls.tsv', sep='\t')
            blocks.to_csv(run.dir + 'bc_inversion_blocks.tsv', sep='\t', index=False)
            print(blocks.to_string(index=False))
            table = bci.apply_to_table(table, calls, bcfg['lgs'])
            exp_chroms = {lg: (LEGACY_MC_ALL[lg] if bc_legacy_mc else proj.chrom_for_lg(lg)) for lg in bcfg['expected_lgs'].values()}
            no_block = bci.check_no_block(grid, exp_chroms, bcfg['max_rate'], bcfg['min_block_mb'] * 1e6)
            print('no recombination-suppressed block (supports the pedigree expectation):', no_block)
            table, filled = bci.apply_expected(table, bcfg['expected'], no_block, bcfg['expected_lgs'])
            print('pedigree-expected columns filled for BC fish:', filled)
            print('BC fish with inferred genotypes:', int(calls.filter(like='Inversion').notna().any(axis=1).sum()), 'of', len(calls))
        table.to_csv(run.dir + 'SampleDatabase_with_inversions.tsv', sep='\t', index=False)
        rep['unmatched'].to_csv(run.dir + 'unmatched_samples.tsv', sep='\t', index=False)
        rep['s1_not_in_db'].to_csv(run.dir + 'blumer_samples_not_in_sampledatabase.tsv', sep='\t', index=False)
        summary = {k: v for k, v in rep.items() if k not in ('unmatched', 's1_not_in_db')}
        open(run.dir + 'join_report.txt', 'w').write('\n'.join(f'{k}: {v}' for k, v in summary.items()) + '\n')
        print(open(run.dir + 'join_report.txt').read())
        if rep['inv10_discordant'] or rep['sex_discordant']:
            print('WARNING: discordances with Table S7/S1; see join_report.txt')
        write_new_columns(proj, table, run.dir, execute=write_sheet)
        if bc_genotypes:
            fill_blank_cells(proj, table, 'Inversion10', run.dir, execute=write_sheet)
