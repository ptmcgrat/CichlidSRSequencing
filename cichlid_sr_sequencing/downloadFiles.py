from helper_modules.file_manager_Replacement import FileManager as FM

fm_obj = FM(genome_version = 'Mzebra_GT3')
fm_obj.createSampleFiles('YH_1_m')
fm_obj.downloadData(fm_obj.localSampleBamDir)
fm_obj.createSampleFiles('MZ_1_m')
fm_obj.downloadData(fm_obj.localSampleBamDir)
# Read and download essential data
main_vcf = fm_obj.localMasterDir + 'Outputs/FilteredFiles/Mzebra_GT3/FilteredFilesGT3Cohort/gt3_cohort_pass_variants.vcf.gz'
fm_obj.downloadData(main_vcf)
fm_obj.downloadData(fm_obj.localSampleFile_v2)