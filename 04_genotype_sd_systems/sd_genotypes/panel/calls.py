"""PILOT calls from STORED measurements: no BAM, no Project, no network. One uniform layer over the four modules so that a writer needs no module-specific code.

    from sd_genotypes.panel import calls
    panel = json.load(open('S.panel.json'))                   # what `panel measure` wrote (schema 1 or 2)
    r = calls.call_panel(panel)                               # {'sample', 'bam', 'bam_size_bytes', 'measured', 'panel_version', 'calls': [...]}
    for c in r['calls']: c['column'], c['value_to_write'], c['flags'], c['confident'], c['call_tag']
    values, flags = calls.sheet_row(r)                        # {column: value} and {flags column: 'pilot-v1;flag;flag'}

A call is a dict: module, column (the sheet column of the value), flags_column, value (what the data say; None if there is no estimate), flags (list), confident (True when there is a value and none of the
module's `blocking_flags` is set and the value is not in `sheet.unvalidated_values`), value_to_write (value if confident else None), call_tag ('pilot-v1': status and version of the call settings), definition (the key of the measurement) and details.
Values: gsdf an integer number of copies; the three transposon sites '0/0', '0/1' or '1/1'; the LG5 inversion 0, 1 or 2 (inverted chromosomes); the B chromosome a float score (no B+/B- call).
Everything that can change lives in config/*.yaml (thresholds, `status`, `version`, `sheet`): changing it never needs a BAM. A measurement made under another DEFINITION than the repository's current one (the key
differs) is refused: call it with the commit that holds that definition, or measure it again."""
import argparse, importlib, json, math
import pandas as pd
from . import identity as ID, measure as pm


class DefinitionMismatch(ValueError):
    pass


class ConfigError(ValueError):
    pass


def _nan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def _num(x, nd=4):
    return None if _nan(x) else (round(float(x), nd) if isinstance(x, float) or hasattr(x, 'dtype') else x)


def call_tag(cfg):
    """'pilot-v1': the status and version of the call settings; it is written first in every flags cell, so a cell can be traced to the settings that made it."""
    for k in ('status', 'version', 'sheet'):
        if k not in cfg:
            raise ConfigError(f"the config has no '{k}': every panel config needs status, version and sheet (see config/*.yaml)")
    return f"{cfg['status']}-v{int(cfg['version'])}"


def _call(module, column, flags_column, value, flags, blocking, cfg, definition, details):
    """One call. A value listed in the config's `sheet.unvalidated_values` (a value never seen on a real fish, e.g. an Inversion5 of 2) keeps its value but gets the flag `unvalidated_value` and is never confident, so it is
    never written: the guard lives here, so that it protects every consumer of the calls and not only one writer (which must take `value_to_write`, never `value`)."""
    flags = list(flags)
    if value is not None and value in cfg['sheet'].get('unvalidated_values', []) and 'unvalidated_value' not in flags:
        flags.append('unvalidated_value')
    confident = value is not None and not (set(flags) & set(blocking)) and 'unvalidated_value' not in flags
    return dict(module=module, column=column, flags_column=flags_column, value=value, flags=flags, confident=bool(confident), value_to_write=value if confident else None,
                call_tag=call_tag(cfg), definition=definition, details=details)


def call_gsdf(block, cfg, context=None):
    from sd_genotypes.gsdf import call as GC, measure as GM
    m = GM.metrics(block, cfg); c = GC.call_one(m, cfg); s = cfg['sheet']
    value = None if _nan(c['gsdf_number']) else int(c['gsdf_number'])
    details = dict(copies_estimate=_num(c['copies_estimate']), genotype=c['genotype'], ratio=_num(m['ratio']), flank_depth=_num(m['flank_depth'], 2), unit_depth=_num(m['unit_depth'], 2), junction_ratio=_num(m['junction_ratio']))
    return [_call('gsdf', s['column'], s['flags_column'], value, list(c['flags']) + list(c['check_flags']), s['blocking_flags'], cfg, block.get('definition'), details)]


def call_te(block, cfg, context=None):
    """The three sites, from the stored read-pair counts. context['gsdf_copies'] (the unflagged gsdf copy number of the same sample) switches on the chr7ins dilution flag."""
    from sd_genotypes.te import call as TC
    copies = (context or {}).get('gsdf_copies')
    row = TC.call_table([dict(SampleID='x', **block)], cfg, copies=None if copies is None else pd.Series({'x': copies})).iloc[0]
    out = []
    for site, spec in cfg['sheet']['sites'].items():
        gt, gq, status = row[f'{site}_GT'], row[f'{site}_GQ'], row[f'{site}_status']
        flags, value = [], None
        if status != 'OK':
            flags.append('site_error')
        elif gt not in ('0/0', '0/1', '1/1'):
            flags.append('no_call')
        else:
            value = gt
            if not gq >= cfg['call']['min_gq']:
                flags.append('low_gq')
        flags += [f for f in str(row.get(f'{site}_flags', '')).split(';') if f]
        details = dict(ref_pairs=_num(row[f'{site}_ref'], 0), alt_pairs=_num(row[f'{site}_alt'], 0), GQ=_num(gq, 1), GT=gt, status=status)
        out.append(_call('te', spec['column'], spec['flags_column'], value, flags, cfg['sheet']['blocking_flags'], cfg, block.get('definition'), details))
    return out


def call_lg5inv(block, cfg, context=None):
    from sd_genotypes.lg5inv import call as LC
    m = LC.metrics(block, cfg); c = LC.call_one(m, cfg); s = cfg['sheet']
    value = None if _nan(c['n_inverted']) else int(c['n_inverted'])
    details = dict(genotype=c['genotype'], median_site_fraction=_num(m['median_site_fraction']), aggregate_fraction=_num(m['aggregate_fraction']), n_usable=m['n_usable'], mean_depth=_num(m['mean_depth'], 1),
                   junction_A=m['junction_A'], junction_B=m['junction_B'], junction_C=m['junction_C'], junction_fraction=_num(m['junction_fraction']))
    return [_call('lg5inv', s['column'], s['flags_column'], value, c['flags'], s['blocking_flags'], cfg, block.get('definition'), details)]


def call_bchrom(block, cfg, context=None):
    from sd_genotypes.bchrom import measure as BM, score as BS
    blocks, _, _ = BM.load_definition(cfg)
    m = BS.metrics(block, cfg, blocks); c = BS.call_one(m, cfg); s = cfg['sheet']
    details = dict(scale_depth=_num(m['scale_depth'], 2), n_blocks_used=m['n_blocks_used'], frac_ge2=_num(m['frac_ge2']), frac_ge3=_num(m['frac_ge3']), p10_scaled=_num(m['p10_scaled']))
    return [_call('bchrom', s['column'], s['flags_column'], _num(m['score']), c['flags'], s['blocking_flags'], cfg, block.get('definition'), details)]


CALLERS = {'gsdf': call_gsdf, 'te': call_te, 'lg5inv': call_lg5inv, 'bchrom': call_bchrom}


def call_module(module, block, cfg=None, context=None, check_definition=True):
    """Calls of ONE module from its stored measurement `block` (panel['modules'][module]); a list (the te module makes three)."""
    if module not in CALLERS:
        raise ValueError(f'no calls for module {module!r}; available: {sorted(CALLERS)}')
    cfg = cfg or pm.module_cfg(module)
    call_tag(cfg)
    if check_definition:
        key = importlib.import_module(pm.MODULES[module][0]).definition_key(cfg)
        if block.get('definition') != key:
            raise DefinitionMismatch(f"{module}: measured under definition {block.get('definition')}, the current definition is {key}: use the commit that holds the definition of the measurement, or measure again")
    return CALLERS[module](block, cfg, context)


def call_panel(panel, modules=None, cfgs=None, check_definition=True):
    """All calls of one panel JSON (a parsed dict, or a path). gsdf is called first: its copy number (when the call is confident) feeds the chr7ins dilution flag."""
    panel = json.load(open(panel)) if isinstance(panel, str) else panel
    name, size = ID.normalize(panel); context, calls = {}, []
    for m in [x for x in ('gsdf', 'te', 'lg5inv', 'bchrom') if x in panel['modules'] and (modules is None or x in modules)]:
        got = call_module(m, panel['modules'][m], (cfgs or {}).get(m), context, check_definition)
        if m == 'gsdf':
            context['gsdf_copies'] = got[0]['value_to_write']
        calls += got
    return dict(sample=panel['sample'], bam=name, bam_size_bytes=size, measured=panel.get('measured'), panel_version=panel.get('panel_version'), calls=calls)


def sheet_row(result, write_unconfident=False):
    """({column: value}, {flags column: 'tag;flag;flag'}) of a call_panel result. A value is included only if the call is confident (or write_unconfident); the flags cell always, starting with the call tag."""
    values, flags = {}, {}
    for c in result['calls']:
        v = c['value'] if write_unconfident else c['value_to_write']
        if v is not None:
            values[c['column']] = v
        flags[c['flags_column']] = ';'.join([c['call_tag']] + c['flags'])
    return values, flags


def main(argv=None):
    ap = argparse.ArgumentParser(prog='python -m sd_genotypes.panel.calls', description='Print the calls of a stored panel JSON (no BAM needed).')
    ap.add_argument('--json', required=True); ap.add_argument('--modules', default='all')
    a = ap.parse_args(argv)
    r = call_panel(a.json, None if a.modules == 'all' else a.modules.split(','))
    print(f"{r['sample']}  (measured on {r['bam']}, {r['bam_size_bytes']} bytes; panel {r['panel_version']})")
    for c in r['calls']:
        print(f"  {c['column']:12s} value {str(c['value']):>8s}  {'confident' if c['confident'] else 'NOT confident':13s} flags: {';'.join(c['flags']) or '-':40s} [{c['call_tag']}]")


if __name__ == '__main__':
    main()
