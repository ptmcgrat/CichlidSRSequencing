from pyfaidx import Fasta
import os, pdb

fasta_dir = '/Users/kmnike/Desktop/S3_genomes/'
yh_good_mito = fasta_dir + 'A_spYH_GT1.fasta'
fasta = Fasta(yh_good_mito)
with open('good_yh_mito.fasta', 'w') as fh:
    fh.write(fasta['NC_027944.1'][:].__str__())


umd2a = '/Users/kmnike/Data/CichlidSequencingData/Genomes/Mzebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna'
fasta_umd2a = Fasta(umd2a)

with open('umd2a_mito.fasta', 'w') as fh:
    fh.write(fasta_umd2a['NC_027944.1'][:].__str__())
# pdb.set_trace()
# with open(fasta_dir + 'mito_summary.tsv', 'a+') as f:
#     for file in os.listdir(fasta_dir):
#         if not file.endswith('.fasta') or file == 'P_nyererei_v2.fasta':
#             continue
#         fasta_obj = Fasta(os.path.join(fasta_dir, file))
#         if os.path.getsize(fasta_dir + 'mito_summary.tsv') == 0:
#             f.write(f"file\tmito_length\n")
#         f.write(f"{file}\t{len(fasta_obj['NC_027944.1'][:].__str__())}\n")