from helper_modules.file_manager_Replacement import FileManager as FM
import pdb, subprocess, sys

# Create fm_obj and grab sample file
fm_obj = FM(genome_version = 'Mzebra_GT3')
fm_obj.downloadData(fm_obj.localGenomeFile)

soft_dir = '/home/ad.gatech.edu/bio-mcgrath-dropbox/Patrick/RepeatAnnotation/RepeatModeler-2.0.5'

subprocess.run([soft_dir + 'BuildDatabase','-name',fm_obj.localGenomeDir + 'Mzebra_GT3', fm_obj.localGenomeFile])
subprocess.run([soft_dir + 'RepeatModeler','-database', 'Mzebra_GT3','-threads','12','-LTRStruct'])
