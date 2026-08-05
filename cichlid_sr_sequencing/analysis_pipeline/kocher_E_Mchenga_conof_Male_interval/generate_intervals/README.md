This directory contains the make_kocher_E_Mchenga_conof_Male_intervals script.
This was used to split the 22 linkage groups of Mzebra_GT3a into 96 roughly equal intervals so that each one is ~10Mbp long
Note that the header.csv is created using the dict file for the genome that was generated usign gatk CreateSequenceDictionary. Esnure that the extra tabs are present in the @HD line along with the 'SO:coordinate' part

dependencies: pyfaidx
python make_kocher_E_Mchenga_conof_Male_intervals.py