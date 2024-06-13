This directory contains the make_GT3_intervals script.
This was used to split the 22 linkage groups of Mzebra_GT3a into 96 roughly equal intervals so that each one is ~10Mbp long
Note that the header.csv is created using the dict file for the genome that was generated usign gatk CreateSequenceDictionary. Esnure that the extra tabs are present in the @HD line along with the 'SO:coordinate' part

Note that the BAM files have been generated on the Mzebra_GT3.fasta genome with the old LG names from UMD2a. Therefore, the scripts will need to run based on these LG names

dependencies: pyfaidx
python make_GT3_intervals.py