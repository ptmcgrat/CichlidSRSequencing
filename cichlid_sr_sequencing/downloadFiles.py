from helper_modules.file_manager_Replacement import FileManager as FM

fm_obj = FM(genome_version = 'Mzebra_GT3')
fm_obj.downloadData(fm_obj.localGenomesDir)
