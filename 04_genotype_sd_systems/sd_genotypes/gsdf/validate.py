"""Validation of gsdf calls: (1) agreement with the values already in the sheet (`gsdf number`), (2) sex x copy number per species. In a species that uses the duplication as a sex determiner the
extra copies sit on the Y: females have two copies (1 on each X), males three (2 on the Y, 1 on the X), as in the MC pedigree (46 of 46 unflagged calls). Four or five copies appear only in males or
unsexed fish, which points to extra tandem copies on the Y, not to dup/dup. Depth gives total copy number only."""
import numpy as np
import pandas as pd


def sheet_frame(df):
    """SampleID, Sex (M / F; uncertain 'M?' / 'F?' and unknown are NaN), Species, Ecogroup, sheet_gsdf (the existing `gsdf number`)."""
    out = pd.DataFrame(dict(SampleID=df['SampleID'].astype(str), Species=df.get('Species'), Ecogroup=df.get('Ecogroup')))
    sex = df['Sex'].astype(str).str.strip() if 'Sex' in df else pd.Series('', index=df.index)
    out['Sex'] = sex.where(sex.isin(['M', 'F']))
    out['sheet_gsdf'] = pd.to_numeric(df['gsdf number'], errors='coerce') if 'gsdf number' in df else np.nan
    return out


def compare_with_sheet(calls, sheet, only_unflagged=True):
    """Agreement of the new calls with the existing column. Returns (summary dict, mismatches DataFrame)."""
    m = calls.merge(sheet[['SampleID', 'sheet_gsdf', 'Species', 'Sex']], on='SampleID', how='inner')
    both = m[m.sheet_gsdf.notna() & m.gsdf_number.notna()]
    use = both[both['call_flags'] == ''] if only_unflagged else both
    agree = use[use.gsdf_number == use.sheet_gsdf]
    summ = dict(n_calls=len(calls), n_with_sheet_value=len(both), n_compared=len(use), n_agree=len(agree), frac_agree=len(agree) / len(use) if len(use) else np.nan,
                n_flagged_among_with_sheet=int((both['call_flags'] != '').sum()))
    return summ, use[use.gsdf_number != use.sheet_gsdf].sort_values('SampleID')


def sex_table(df, min_per_sex=3, male_min_copies=3, frac=0.9):
    """df: SampleID, Sex, Species, Ecogroup, copies. One row per species with >= min_per_sex sexed males AND females."""
    rows, disc = [], []
    for sp, g in df.dropna(subset=['copies', 'Species']).groupby('Species'):
        f, m = g[g.Sex == 'F'], g[g.Sex == 'M']
        if len(f) < min_per_sex or len(m) < min_per_sex:
            continue
        ff, mf = float((f.copies == 2).mean()), float((m.copies >= male_min_copies).mean())
        any_dup = bool((g.copies >= 3).any())
        pattern = 'no duplication seen' if not any_dup else ('male-specific (XY-like)' if ff >= frac and mf >= frac else 'duplication, not male-specific')
        cnt = lambda s: ';'.join(f'{int(k)}:{v}' for k, v in s.copies.astype(int).value_counts().sort_index().items())
        rows.append(dict(Species=sp, Ecogroup=g.Ecogroup.mode().iloc[0] if g.Ecogroup.notna().any() else '', n_F=len(f), n_M=len(m), F_copies=cnt(f), M_copies=cnt(m),
                         frac_F_two_copies=ff, frac_M_ge_min=mf, pattern=pattern))
        if pattern == 'male-specific (XY-like)':                                       # exceptions in an XY-like species
            disc.append(g[((g.Sex == 'F') & (g.copies >= 3)) | ((g.Sex == 'M') & (g.copies == 2))].assign(why=lambda x: np.where(x.Sex == 'F', 'female with the male-specific copy', 'male without it')))
    t = pd.DataFrame(rows, columns=['Species', 'Ecogroup', 'n_F', 'n_M', 'F_copies', 'M_copies', 'frac_F_two_copies', 'frac_M_ge_min', 'pattern'])
    return t.sort_values(['pattern', 'Species']).reset_index(drop=True), (pd.concat(disc) if disc else pd.DataFrame(columns=list(df.columns) + ['why']))
