import argparse
import shutil
import pathlib
import pdb

bad_samples = [
    'MCYHBC1-548-1', 'MCYHBC1-591-1', 'MCYHBC1-605-1', 'MCYHBC1-622-1',
    'MCYHBC1-646-1', 'MCYHBC1-663-1', 'MCYHBC1-669-1', 'MCYHBC1-671-1',
    'MCYHBC1-711-1', 'MCYHBC1-712-1', 'MCYHBC1-717-1', 'MCYHBC1-718-1',
    'MCYHBC1-720-1', 'MCYHBC1-731-1', 'MCYHBC1-734-1', 'MCYHBC1-751-1',
    'MCYHBC1-773-1', 'MCYHBC1-775-1', 'MCYHBC1-790-1', 'MCYHBC1-795-1',
    'MCYHBC1-799-1', 'MCYHBC1-801-1', 'MCYHBC1-802-1', 'MCYHBC1-803-1',
    'MCYHBC1-806-1', 'MCYHBC1-810-1', 'MCYHBC1-813-1', 'MCYHBC1-835-1',
    'MCYHBC1-839-1', 'MCYHBC1-841-1', 'MCYHBC1-855-1', 'MCYHBC1-918-1',
    'MCYHBC1-922-1', 'MCYHBC1-935-1'
]


temp_base    = pathlib.Path('/Data/mcgrath-lab/Temp/CichlidSequencingData/Temp')
bamfile_base = pathlib.Path('/Data/mcgrath-lab/Temp/CichlidSequencingData/Bamfiles/Mconophoros_GT1')

parser = argparse.ArgumentParser()
parser.add_argument('--dry-run', action='store_true', help='Print what would be deleted without deleting')
args = parser.parse_args()

for sample in bad_samples:
    for base in [temp_base, bamfile_base]:
        sample_dir = base / sample
        if sample_dir.exists():
            if args.dry_run:
                print(f'DRY RUN - would delete: {sample_dir}')
            else:
                shutil.rmtree(sample_dir)
                print(f'Deleted: {sample_dir}')
        else:
            print(f'NOT FOUND: {sample_dir}')