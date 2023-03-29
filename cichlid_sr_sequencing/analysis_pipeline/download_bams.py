from helper_modules.nikesh_file_manager import FileManager as FM
import subprocess as sp, shlex

fm_obj = FM('Mzebra_UMD2a')
samples = ['MC_1_m', 'SAMEA2661294', 'SAMEA2661322', 'SAMEA4032100', 'SAMEA4033261']
for sample in samples:
    fm_obj.createSampleFiles(sample)
    sp.run(shlex.split(f'rclone copy ptm_dropbox:{fm_obj.locaBamFile} .'))
    sp.run(shlex.split(f'rclone copy ptm_dropbox:{fm_obj.locaBamIndexFile} .'))