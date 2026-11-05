"""python unit_scripts/split_bamfile_single_pass.py <sample.all.bam> <output prefix> [--threads 2]
Writes <prefix>.{unmapped,discordant,inversion,duplication,clipped,chimeric}.bam from a coordinate-sorted, indexed BAM in ONE pass, plus <prefix>.split_summary.json.
Index the outputs afterwards (samtools index)."""
import argparse, json, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from helper_modules.bam_split import split_bam

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('Bamfile'); parser.add_argument('Prefix'); parser.add_argument('--threads', type=int, default=2)
args = parser.parse_args()
counts = split_bam(args.Bamfile, args.Prefix, args.threads)
with open(args.Prefix + '.split_summary.json', 'w') as f:
    json.dump(counts, f, indent=1)
print(json.dumps(counts))
