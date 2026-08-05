import pdb
import argparse, os, subprocess, datetime, time, multiprocessing, pathlib, random
import pandas as pd
from pyfaidx import Fasta
from helper_modules.nikesh_file_manager import FileManager as FM

def get_sampleids():
    """
    Get sample IDs from the input directory, filtering for files ending in .fa or .fasta.
    """
    sampleIDs = []
    for file in os.listdir(args.fasta_input_dir):
        if file.endswith('.fa') or file.endswith('.fasta'):  # Filter for .fa and .fasta files
            sampleIDs.append(file.split('.')[0])  # Extract the sample ID (filename without extension)
    return sampleIDs

def minimapAlign(sample):
    """
    Align a single sample using minimap2 and convert SAM to BAM if necessary.
    """
    # Create the pipeline directory within alignments
    alignment_dir = fm_obj.localOutputDir + 'alignments/pipeline/' + sample
    pathlib.Path(alignment_dir).mkdir(parents=True, exist_ok=True)

    # Output SAM, BAM, and RG BAM file paths
    sam_file_path = alignment_dir + '/' + sample + '.sam'
    bam_file_path = alignment_dir + '/' + sample + '.bam'
    sorted_bam_file_path = alignment_dir + '/' + sample + '.sorted.bam'
    rg_bam_file_path = alignment_dir + '/' + sample + '.rg.bam'

    # Check if the SAM file exists and has a non-zero size
    if os.path.exists(sam_file_path) and os.path.getsize(sam_file_path) > 0:
        print(f"SAM file for sample {sample} already exists and is non-zero in size. Skipping minimap2.")
    else:
        # Run minimap2 to generate the SAM file
        with open(sam_file_path, 'w') as f:
            try:
                subprocess.run(
                    ['minimap2', fm_obj.localGenomeFile, args.fasta_input_dir + '/' + sample + '.fa', '-a', '-t', '48'],  # Use 48 threads per sample
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True
                )
            except subprocess.CalledProcessError as e:
                print(f"Error running minimap2 for sample {sample}: {e.stderr.decode()}")
                raise

    # Check if the BAM file exists; if not, create it using samtools
    if not os.path.exists(bam_file_path):
        print(f"BAM file for sample {sample} does not exist. Creating it from the SAM file.")
        try:
            subprocess.run(
                ['samtools', 'view', '-bS', sam_file_path, '-o', bam_file_path],
                stderr=subprocess.PIPE,
                check=True
            )
            subprocess.run(
                ['samtools', 'sort', '-o', sorted_bam_file_path, bam_file_path],
                stderr=subprocess.PIPE,
                check=True
            )
            subprocess.run(
                ['samtools', 'index', sorted_bam_file_path],
                stderr=subprocess.PIPE,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error creating BAM file for sample {sample}: {e.stderr.decode()}")
            raise
    else:
        print(f"BAM file for sample {sample} already exists. Skipping BAM creation.")

    # Add read groups to the sorted BAM file
    if not os.path.exists(rg_bam_file_path):
        print(f"Adding read groups to BAM file for sample {sample}.")
        try:
            subprocess.run(
                [
                    'gatk', 'AddOrReplaceReadGroups',
                    '-I', sorted_bam_file_path,
                    '-O', rg_bam_file_path,
                    '-RGID', sample,
                    '-RGLB', 'lib1',
                    '-RGPL', 'pacbio',
                    '-RGPU', 'unit1',
                    '-RGSM', sample
                ],
                stderr=subprocess.PIPE,
                check=True
            )
            # Index the .rg.bam file
            subprocess.run(
                ['samtools', 'index', rg_bam_file_path],
                stderr=subprocess.PIPE,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error adding read groups or indexing for sample {sample}: {e.stderr.decode()}")
            raise
    else:
        print(f"Read group BAM file for sample {sample} already exists. Skipping read group addition.")

def runHaplotypeCaller(sampleIDs):
    for sample in sampleIDs:
        # Create the pipeline directory within alignments
        alignment_dir = fm_obj.localOutputDir + 'alignments/pipeline/' + sample
        pathlib.Path(alignment_dir).mkdir(parents=True, exist_ok=True)

        # Output GVCF file path
        gvcf_file_path = alignment_dir + '/' + sample + '.g.vcf'

        # Ensure the .rg.bam file exists before running HaplotypeCaller
        rg_bam_file_path = alignment_dir + '/' + sample + '.rg.bam'
        if not os.path.exists(rg_bam_file_path):
            print(f"Error: Read group BAM file for sample {sample} does not exist. Skipping HaplotypeCaller.")
            continue

        # Run GATK HaplotypeCaller to generate GVCF
        with open(gvcf_file_path, 'w') as f:
            try:
                subprocess.run(
                    [
                        'gatk', 'HaplotypeCaller',
                        '-R', fm_obj.localGenomeFile,
                        '-I', rg_bam_file_path,
                        '-O', gvcf_file_path,
                        '-ERC', 'GVCF'  # Emit GVCF format
                    ],
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True
                )
            except subprocess.CalledProcessError as e:
                print(f"Error running HaplotypeCaller for sample {sample}: {e.stderr.decode()}")
                raise


def jointGenotyping(sampleIDs):
    # Create the pipeline directory within alignments
    alignment_dir = fm_obj.localOutputDir + 'alignments/pipeline/'
    pathlib.Path(alignment_dir).mkdir(parents=True, exist_ok=True)

    # Output joint VCF file path
    joint_vcf_path = alignment_dir + 'jointGenotyping.vcf'

    with open(joint_vcf_path, 'w') as f:
        try:
            subprocess.run(
                ['gatk', 'GenotypeGVCFs', '-R', fm_obj.localGenomeFile, '-V', alignment_dir + sampleIDs[0] + '/' + sampleIDs[0] + '.vcf', '-O', joint_vcf_path],
                stdout=f,
                stderr=subprocess.PIPE,
                check=True
            )
        except subprocess.CalledProcessError as e:
            print(f"Error running joint genotyping: {e.stderr.decode()}")
            raise

def multiprocess(function):
    """
    Multiprocess function to utilize available cores for processing samples.
    """
    inputs = get_sampleids()
    num_samples = len(inputs)
    concurrent_processes = min(96, num_samples)  # Use up to 96 cores or the number of samples, whichever is smaller

    print(f"Running {function.__name__} with {concurrent_processes} processes for {num_samples} samples.")

    try:
        with multiprocessing.Pool(processes=concurrent_processes) as pool:
            pool.map(function, inputs)
    except Exception as e:
        print(f"Error occurred during multiprocessing: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser('YH lg 10 analysis script for aligning with minimap2 and calling variants using haplotypecaller')
    parser.add_argument("Genome", help='reference genome')
    parser.add_argument('-f', '--fasta_input_dir', help='directory with fasta input files to align')
    parser.add_argument('-a', '--align', help='Call flag to run minimap alignment on input samples', action='store_true')
    parser.add_argument('-g', '--gatk', help='Call flag to run GATK HaplotypeCaller on input samples', action='store_true')
    parser.add_argument('--local_test', help='when this flag is called, variables will be preset to test the code locally', action='store_true')
    parser.add_argument('-t', '--threads', type=int, default=96, help='Number of threads to use')
    args = parser.parse_args()

    fm_obj = FM(args.Genome)
    if args.local_test:
        fm_obj.localGenomeFile = '/Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/alignments/target/small_target.fa'

    # Calling the functions
    if args.align:
        sampleIDs = get_sampleids()
        print(f"Found {len(sampleIDs)} samples: {sampleIDs}")

        # Use multiprocessing to run minimapAlign concurrently for all samples
        with multiprocessing.Pool(processes=min(len(sampleIDs), args.threads // 48)) as pool:
            pool.map(minimapAlign, sampleIDs)
    
    if args.gatk:
        sampleIDs = get_sampleids()
        print(f"Found {len(sampleIDs)} samples: {sampleIDs}")

        # Temporarily disable multiprocessing for debugging
        for sample in sampleIDs:
            runHaplotypeCaller([sample])  # Pass a single sample as a list to match the function signature

        # Run joint genotyping (if needed)
        jointGenotyping(sampleIDs)


# Run alignment and GATK Haplotypecaller
    # time python analyzeYH10.py Mzebra_GT3 -f /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/alignments/queries -a -g
# Run alignment only
    # time python analyzeYH10.py Mzebra_GT3 -f /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/alignments/queries -a
# Run GATK Haplotypecaller and GenotypeGVCFs only
    # time python analyzeYH10.py Mzebra_GT3 -f /Data/mcgrath-lab/Data/CichlidSequencingData/Outputs/alignments/queries -g