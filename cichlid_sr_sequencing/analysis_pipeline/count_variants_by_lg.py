import argparse, os, pdb
import subprocess as sp

parser = argparse.ArgumentParser(usage='count variants per linkage group in input PCA directory')
parser.add_argument('-i', '--input_dir', help = 'location of dir containing the split lgs')
parser.add_argument('-o', '--out', help = 'name of file to which the variant counts are written', default='lg_variants_counts.tsv')
parser.add_argument('--single', help = 'path to a single gzipped file on which to operate and count linkage groups', nargs = 1, default=False)
args = parser.parse_args()

class Counter:
    def __init__(self, input_directory, output_file, single_file):
        self.in_dir = input_directory
        self.out = output_file
        self.linkage_groups = os.listdir(self.in_dir)
        self.linkage_groups.sort()
        self.single = single_file
        if self.single:
             self.out = 'single_file_counts.tsv'
             self._count_lgs_in_file(self.single[0])
        else:
            self._count_variants(self.linkage_groups)

    def _count_variants(self, linkage_group_list):
        with open(self.out, 'w') as f:
            for lg in linkage_group_list:
                print(f'bcftools is reading in and counting variants for {lg}')
                self.input_file = self.in_dir + '/' + lg + '/' + lg + '.vcf.gz'
                count = sp.run(f"bcftools view {self.input_file} | wc -l", shell=True, capture_output=True, text=True)
                f.write(f"{lg}\t{count.stdout.strip()}\n")
                print(f"{lg}\t{count.stdout.strip()}")

    def _count_lgs_in_file(self, single_input_file):
        lgs = ['NC_036780.1', 'NC_036781.1', 'NC_036782.1', 'NC_036783.1', 'NC_036784.1', 'NC_036785.1', 'NC_036786.1', 'NC_036787.1', 'NC_036788.1', 'NC_036789.1', 'NC_036790.1', 'NC_036791.1', 'NC_036792.1', 'NC_036793.1', 'NC_036794.1', 'NC_036795.1', 'NC_036796.1', 'NC_036797.1', 'NC_036798.1', 'NC_036799.1', 'NC_036800.1', 'NC_036801.1']
        with open(self.out, 'w') as f:
            for lg in lgs:
                print(f'bcftools is reading in and counting variants for {lg}')
                count = sp.run(f'bcftools view {single_input_file} -r {lg} | wc -l', shell=True, capture_output=True, text=True)
                f.write(f"{lg}\t{count.stdout.strip()}\n")
                print(f"{lg}\t{count.stdout.strip()}")

    def _count_variants_in_parallel(self):
        #### IN DEVELOPEMNT
        processes = []
        with open(self.out, 'w') as f:
                    for lg in self.linkage_groups:
                        print(f'bcftools is reading in and counting variants for {lg}')
                        self.input_file = self.in_dir + '/' + lg + '/' + lg + '.vcf.gz'
                        p1 = sp.Popen(f"bcftools view {self.input_file} | wc -l", shell=True, stdout=sp.PIPE, text=True)
                        processes.append(p1)
                        if len(processes) == 4:
                            # pdb.set_trace()
                            for proc in processes:
                                count = proc.communicate()
                                f.write(f"{lg}\t{count[0].strip()}\n")
                                print(f"{lg}\t{count[0].strip()}")
                            processes = []

Counter(args.input_dir, args.out, args.single)