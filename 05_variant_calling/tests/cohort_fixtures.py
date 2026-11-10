"""Fake SampleDatabase and AlignmentDatabase tables for the cohort tests."""
import pandas as pd

G = 'Mzebra_GT3_NCBI'
COLS = ['SampleID', 'Sex', 'Species', 'BroodID', 'Parents', 'Ecogroup', 'Subgroup', 'Category', 'LabReared', 'Project']
YH = 'Aulonocara Nhkata Yellow Head'


def fish(sid, sex='M', species=YH, brood='', parents='', eco='Shallow_Benthic', sub='', cat='', lab='1', project='P1'):
    return dict(SampleID=sid, Sex=sex, Species=species, BroodID=brood, Parents=parents, Ecogroup=eco, Subgroup=sub, Category=cat, LabReared=lab, Project=project)


def brood(b='B1', n_m=12, n_f=12, sire='S_m', dam='D_f', species=YH, parents_in_sheet=True):
    rows = []
    if parents_in_sheet:
        rows += [fish(sire, 'M', species), fish(dam, 'F', species)]
    for i in range(n_m):
        rows.append(fish(f'{b}_m{i}', 'M', species, b, f'{sire},{dam}'))
    for i in range(n_f):
        rows.append(fish(f'{b}_f{i}', 'F', species, b, f'{sire},{dam}'))
    return rows


def sheet(rows):
    return pd.DataFrame(rows, columns=COLS)


def ledger(ids, cov=20.0, genome=G, **overrides):
    """one row per id; overrides: {id: coverage or None}"""
    rows = []
    for s in ids:
        c = overrides.get(s, cov)
        rows.append(dict(SampleID=s, GenomeVersion=genome, Coverage='' if c is None else c))
    return pd.DataFrame(rows, columns=['SampleID', 'GenomeVersion', 'Coverage'])


def family_spec(**kw):
    spec = dict(name='fam', kind='family', genome=G, caller='bcftools', groups='per_sample', broods=['B1'], include_parents=True, exclude_ecogroups=['Error', 'Exclude'])
    spec.update(kw)
    return spec


def sexset_spec(**kw):
    spec = dict(name='ss', kind='sex_set', genome=G, caller='bcftools', groups='per_sample', species=YH, exclude_ecogroups=['Error', 'Exclude'])
    spec.update(kw)
    return spec


def bc_spec(samples, **kw):
    spec = dict(name='bc', kind='bc_parents', genome=G, caller='bcftools', groups='per_sample', samples=list(samples), exclude_ecogroups=['Error', 'Exclude'])
    spec.update(kw)
    return spec
