"""Alignment of one registered assembly to another (minimap2, with cigar), shared by the landmarks (01) and the B-chromosome comparison (04)."""
import os
from .common_light import sh
from .genomes_registry import fasta_path


def ensure_paf(proj, qkey, ref_key, mm_cfg):
    """minimap2 alignment (with cigar) of a registered assembly to the reference, kept in the project's tmp folder (not uploaded); reused if it already exists."""
    full = f"{proj.dirs['tmp']}{qkey}_vs_{ref_key}.paf"
    if not os.path.exists(full):
        extra = ' '.join(mm_cfg.get('extra', []))
        print(f'[align] {qkey} -> {ref_key} (minimap2 -x {mm_cfg["preset"]}, {proj.cores} threads) ...')
        sh(['bash', '-o', 'pipefail', '-c', f"minimap2 -x {mm_cfg['preset']} -c -t {proj.cores} {extra} {fasta_path(proj, ref_key)} {fasta_path(proj, qkey)} > {full}.part && mv {full}.part {full}"])
    return full
