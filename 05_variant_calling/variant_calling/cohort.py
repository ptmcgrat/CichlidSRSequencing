"""Cohort definitions and the members they select (module 05).

A cohort is an entry of config/callsets.yaml, the file the GATK `callset build` also reads (it ignores the keys it does not know):

    yh_broods:
      kind: family               # family | sex_set | bc_parents: decides the rules of `cohort check` (SPEC decision 61)
      genome: Mzebra_GT3_NCBI    # the reference the BAMs were aligned to (default)
      caller: bcftools           # bcftools | gatk (decision 59)
      groups: per_sample         # Hardy-Weinberg prior groups of `bcftools call -G`: per_sample | ecogroup | species (decision 63)

Selectors (any mix; the result is the union, in this order, without repeats):
    samples: [ID, ...]                      explicit SampleIDs
    broods: [BroodID, ...]                  the brood's offspring, and (include_parents, default true) their listed parents
    species / ecogroups / subgroups / categories
                                            every sample matching ALL the given criteria. `species` alone selects by species only when
                                            there are no `broods` and no `samples`; next to them it is the species the cohort DECLARES
                                            (the check blocks members of another species).
    exclude_samples: [ID, ...]              removed at the end
    exclude_ecogroups: [Error, Exclude]     skipped from criteria-selected samples (as the alignment pipeline does); explicit IDs are kept and flagged
    require_sexed: true                     criteria-selected samples must have Sex M or F (default true for kind sex_set)
"""
import os
import re

import pandas as pd
import yaml

from cichlid_common.common_light import config_path

KINDS = ('family', 'sex_set', 'bc_parents')
CALLERS = ('bcftools', 'gatk')
GROUPS = ('per_sample', 'ecogroup', 'species')
SELECTORS = {'samples', 'broods', 'include_parents', 'species', 'ecogroups', 'subgroups', 'categories', 'exclude_samples', 'exclude_ecogroups', 'require_sexed'}
OPTIONS = {'kind', 'genome', 'caller', 'groups', 'description', 'name'}
DEFAULTS = dict(genome='Mzebra_GT3_NCBI', caller='bcftools', groups='per_sample', include_parents=True, exclude_ecogroups=['Error', 'Exclude'])
SHEET_COLUMNS = ('SampleID', 'Sex', 'Species', 'BroodID', 'Parents', 'Ecogroup', 'Subgroup', 'Category')


# ------------------------------------------------------------------ small helpers shared with the rules
def blank(x):
    return x is None or (isinstance(x, float) and x != x) or str(x).strip().lower() in ('', 'nan', 'none')


def parse_sex(x):
    """1 for M, 0 for F (any case), None for everything else (U, F?, M?, blank): the strictness of sexloci's preflight."""
    if blank(x):
        return None
    s = str(x).strip().upper()
    return 1 if s == 'M' else 0 if s == 'F' else None


_UNKNOWN = {'', '?', 'na', 'nan', 'none', 'unknown', 'unk', '-', '.'}


def parse_parents(x, known=None):
    """'sire,dam' -> ('ok', sire, dam); blank -> ('none', None, None); 'sire,' or 'sire,?' -> ('partial', sire, None); ',dam' -> ('partial', None, dam);
    a single ID that is a known sample -> ('partial_single', id, None); anything else -> ('malformed', None, None). The rules of sexloci's `parse_parents`."""
    if blank(x):
        return 'none', None, None
    toks = [t.strip() for t in str(x).split(',')]
    if len(toks) == 2:
        a, b = [None if t.lower() in _UNKNOWN else t for t in toks]
        if a and b:
            return 'ok', a, b
        if a or b:
            return 'partial', a, b
        return 'none', None, None
    if len(toks) == 1 and known is not None and toks[0] in known:
        return 'partial_single', toks[0], None
    return 'malformed', None, None


def norm_species(x):
    return '' if blank(x) else ' '.join(str(x).split()).lower()


def norm_key(s):
    """Letters and digits only, lower case: finds `YH-008-m` for `YH_008_m`."""
    return re.sub(r'[^a-z0-9]', '', str(s).lower())


# ------------------------------------------------------------------ the definitions
def load_cohorts(path=None):
    """{name: spec}: every entry of `callsets:` with the `defaults: cohort:` section and the module defaults filled in."""
    raw = yaml.safe_load(open(path or config_path('callsets.yaml'))) or {}
    base = dict(DEFAULTS)
    base.update(((raw.get('defaults') or {}).get('cohort')) or {})
    out = {}
    for name, spec in (raw.get('callsets') or {}).items():
        merged = dict(base)
        merged.update(spec or {})
        merged['name'] = name
        out[name] = merged
    return out


def validate(spec):
    """The problems of a definition, one sentence each (empty when it is usable)."""
    bad = []
    if spec.get('kind') not in KINDS:
        bad.append(f"kind is {spec.get('kind')!r}: it must be one of {', '.join(KINDS)}")
    if spec.get('caller') not in CALLERS:
        bad.append(f"caller is {spec.get('caller')!r}: it must be one of {', '.join(CALLERS)}")
    if spec.get('groups') not in GROUPS:
        bad.append(f"groups is {spec.get('groups')!r}: it must be one of {', '.join(GROUPS)}")
    unknown = sorted(set(spec) - SELECTORS - OPTIONS)
    if unknown:
        bad.append('unknown keys (a typo?): ' + ', '.join(unknown))
    for key in ('samples', 'broods', 'ecogroups', 'subgroups', 'categories', 'exclude_samples', 'exclude_ecogroups'):
        if key in spec and spec[key] is not None and not isinstance(spec[key], (list, tuple)):
            bad.append(f'{key} must be a list')
    if not any(spec.get(k) for k in ('samples', 'broods', 'species', 'ecogroups', 'subgroups', 'categories')):
        bad.append('it selects nothing: give samples, broods, species, ecogroups, subgroups or categories')
    return bad


# ------------------------------------------------------------------ the members
def resolve_members(dt, spec):
    """(ordered unique SampleIDs, notes) for a cohort definition, from the SampleDatabase table `dt`. IDs that are not in the table are kept (the check reports them)."""
    notes = []
    d = dt.copy()
    for c in SHEET_COLUMNS:
        if c not in d.columns:
            d[c] = ''
    d = d.drop_duplicates('SampleID', keep='first')
    ids = list(spec.get('samples') or [])
    if spec.get('broods'):
        b = d[d.BroodID.isin(spec['broods'])]
        ids += b.SampleID.tolist()
        if spec.get('include_parents', True):
            for ps in b.Parents:
                _, sire, dam = parse_parents(ps)
                ids += [p for p in (sire, dam) if p]
    by_criteria = any(spec.get(k) for k in ('ecogroups', 'subgroups', 'categories')) or (spec.get('species') and not spec.get('broods') and not spec.get('samples'))
    if by_criteria:
        m = pd.Series(True, index=d.index)
        if spec.get('species'):
            m &= d.Species.map(norm_species) == norm_species(spec['species'])
        for col, key in (('Ecogroup', 'ecogroups'), ('Subgroup', 'subgroups'), ('Category', 'categories')):
            if spec.get(key):
                m &= d[col].isin(spec[key])
        rows = d[m]
        flagged = rows[rows.Ecogroup.isin(spec.get('exclude_ecogroups') or [])]
        if len(flagged):
            notes.append(f"{len(flagged)} samples skipped because their Ecogroup is {'/'.join(spec['exclude_ecogroups'])}: " + ', '.join(flagged.SampleID.head(6)) + (' ...' if len(flagged) > 6 else ''))
            rows = rows.drop(flagged.index)
        want_sexed = spec.get('require_sexed')
        if want_sexed is None:
            want_sexed = spec.get('kind') == 'sex_set'
        if want_sexed:
            unsexed = rows[rows.Sex.map(parse_sex).isna()]
            if len(unsexed):
                notes.append(f'{len(unsexed)} samples skipped because Sex is not M or F: ' + ', '.join(unsexed.SampleID.head(6)) + (' ...' if len(unsexed) > 6 else ''))
                rows = rows.drop(unsexed.index)
        ids += rows.SampleID.tolist()
    drop = set(spec.get('exclude_samples') or [])
    seen, out = set(), []
    for s in ids:
        if s not in seen and s not in drop:
            seen.add(s)
            out.append(s)
    return out, notes
