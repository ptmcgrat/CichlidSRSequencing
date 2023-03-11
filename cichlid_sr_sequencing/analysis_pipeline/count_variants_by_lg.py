import argparse, os, pdb
import subprocess as sp


parser = argparse.ArgumentParser(usage='count variants per linkage group in input PCA directory')
parser.add_argument('input_dir', help = 'location of dir containing the split lgs')
parser.add_argument('-o', '--out', help = 'name of file to which the variant counts are written', default='lg_variants_counts.txt')
args = parser.parse_args()


class Counter:
    def __init__(self, input_directory, output_file):
        self.in_dir = input_directory
        self.out = output_file
        self.linkage_groups = os.listdir(self.in_dir)
        self.linkage_groups.sort()
        self._unzip_input_vcf()

    def _unzip_input_vcf(self):
        with open(self.out, 'w') as f:
            for lg in self.linkage_groups:
                self.input_file = self.in_dir + '/' + lg + '/' + lg + '.vcf.gz'
                count = sp.run(f"bcftools view {self.input_file} | wc -l", shell=True, capture_output=True, text=True)
                f.write(f"{lg}\t{count.stdout.strip()}\n")
                print(f"{lg}\t{count.stdout.strip()}")

Counter(args.input_dir, args.out)