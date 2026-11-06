"""python -m cichlid_common chroms | genomes fetch|annotate|list | selftest"""
import os
from . import cli_base


def _chroms(proj, a, cfg):
    print(proj.chrom_table.to_string(index=False))


def _genomes(proj, a, cfg):
    from . import genomes_registry as gr
    if a.step == 'fetch':
        gr.fetch(proj, a.only.split(',') if a.only else None, upload=cfg['runner']['upload'])
    elif a.step == 'annotate':
        reg = gr.load_registry()
        for key in (a.only.split(',') if a.only else [k for k, v in reg.items() if (v.get('annotation') or {}).get('source_dir')]):
            spec = (reg[key].get('annotation') or {})
            if not spec.get('source_dir'):
                print(f'[annotation] {key}: no source_dir in genomes.yaml; skipped'); continue
            gr.register_annotation(proj, key, spec, upload=cfg['runner']['upload'])
    elif a.step == 'list':
        import pandas as pd
        gr.write_table(proj); print(pd.read_csv(proj.genomes_dir + 'genomes.tsv', sep='\t').to_string(index=False))
    else:
        raise SystemExit('usage: python -m cichlid_common genomes fetch [--only KEY,...] | genomes annotate [--only KEY,...] | genomes list')


def main(argv=None):
    return cli_base.main('python -m cichlid_common', ['chroms', 'genomes'], {'chroms': _chroms, 'genomes': _genomes}, tests_dir=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'tests'), argv=argv)
