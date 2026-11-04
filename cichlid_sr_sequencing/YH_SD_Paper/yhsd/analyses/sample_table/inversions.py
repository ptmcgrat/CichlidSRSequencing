"""Add Blumer et al. 2025 inversion genotypes (Table S7) to the SampleDatabase.
Join key: SampleDatabase.SampleID (BioSample, e.g. SAMEA117805487)  <->  Table S1 `Biosample`  ->  `sequence_id`  <->  Table S7 `sample_id`.
Table S7 coding: 0 = homozygous non-inverted, 1 = heterozygous, 2 = homozygous inverted (relative to fAstCal1.2).
Pure functions (no FileManager / network) so they can be tested on local files."""
import numpy as np
import pandas as pd

INV_CHROMS = ['2', '9', '10', '11', '13']
NEW_COLS = ['Inversion2', 'Inversion9', 'Inversion11', 'Inversion13']


def drop_blank_columns(db):
    """The Google-sheet export carries empty 'Unnamed: N' columns; drop them (never write these back to the sheet)."""
    keep = [c for c in db.columns if not (str(c).startswith('Unnamed') and db[c].isna().all())]
    return db[keep]


def annotate(db, s1, s7):
    """Return (annotated table, report dict). Inversion10 is NOT overwritten; it is validated against Table S7."""
    db = drop_blank_columns(db).copy()
    assert db.SampleID.is_unique, 'SampleDatabase.SampleID is not unique'
    assert s1.Biosample.is_unique, 'Table S1 Biosample is not unique'
    m = s1[['sequence_id', 'specimen_id', 'Biosample', 'sex']].merge(s7, left_on='sequence_id', right_on='sample_id', how='left', validate='1:1')
    m = m.drop(columns='sample_id')
    out = db.merge(m, left_on='SampleID', right_on='Biosample', how='left', validate='m:1', suffixes=('', '_S1'))
    for c in INV_CHROMS:
        out[f'S7_chr{c}'] = out[f'chr{c}']
    out = out.rename(columns={'sequence_id': 'BlumerSequenceID', 'specimen_id': 'BlumerSpecimenID'})
    for c in ['2', '9', '11', '13']:
        out[f'Inversion{c}'] = out[f'chr{c}']
    out['InversionSource'] = np.where(out.chr10.notna(), 'Blumer2025_TableS7', '')
    matched = out.BlumerSequenceID.notna()
    rep = {}
    both = out[out.Inversion10.notna() & out.chr10.notna()]
    rep['n_db'] = len(db); rep['n_matched'] = int(matched.sum())
    rep['inv10_compared'] = len(both); rep['inv10_discordant'] = int((both.Inversion10 != both.chr10).sum())
    rep['inv10_discordant_ids'] = both[both.Inversion10 != both.chr10].SampleID.tolist()
    sx = out[matched & out.Sex.notna() & out.sex.notna()]
    rep['sex_compared'] = len(sx); rep['sex_discordant'] = int((sx.Sex != sx.sex).sum())
    rep['inv10_filled_from_S7'] = int((out.Inversion10.isna() & out.chr10.notna()).sum())       # DB blank, S7 has a value
    rep['inv10_without_S7'] = int((out.Inversion10.notna() & out.chr10.isna()).sum())          # lab-genotyped, not in S7
    rep['unmatched'] = out[~matched]
    rep['s1_not_in_db'] = s1[~s1.Biosample.isin(db.SampleID)]
    keep_cols = list(db.columns) + ['BlumerSequenceID', 'BlumerSpecimenID'] + NEW_COLS + ['InversionSource']
    return out[keep_cols], rep


def fill_inversion10(table):
    """Optional: where the DB has no Inversion10 but Table S7 does, fill it (kept separate so it is a deliberate choice)."""
    return table
