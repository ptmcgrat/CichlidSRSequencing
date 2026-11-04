"""OPTIONAL: keep a cloud copy of the old MC-based (Mconophoros_GT1) results under YH_SD_Paper/legacy_Mconophoros_GT1/ for
reference. Sources are never modified or deleted. Default is a dry run; pass --execute to copy."""
from ...common_light import sh

FILES = ['OutputGenotypes.csv', 'RecombinantsByChromosome.pdf', 'AllBroodParents_Filtered.vcf.gz', 'BC_Pileup_Samples.txt']


def main(proj, execute=False):
    src = proj.fm.cloudMasterDir + '2026_Nikesh/'
    dst = proj.fm.cloudMasterDir + proj.cfg['project_dir_name'] + '/legacy_Mconophoros_GT1/'
    flag = [] if execute else ['--dry-run']
    for f in FILES:
        r = sh(['rclone', 'copyto', src + f, dst + f, '-v'] + flag, check=False)
        print(('COPIED ' if execute else 'WOULD COPY ') + f, '' if r.returncode == 0 else 'ERROR: ' + r.stderr[-200:])
