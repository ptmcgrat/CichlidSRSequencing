import argparse, pdb, pysam, subprocess, os
from helper_modules.file_manager import FileManager as FM

from helper_modules.alignment_worker import AlignmentWorker as AW
from helper_modules.Timer import Timer

fm_obj = FM()

# Need to make SampleIDs and ProjectIDs mutually exclusive
parser = argparse.ArgumentParser(usage = 'This script will download fastq data the McGrath lab dropbox and align it to the Genome version of choice. It will also create gvcf files')
parser.add_argument('Genome', type = str, choices = fm_obj.returnOptions('Genomes'), help = 'Version of the genome to align to')
parser.add_argument('-n', '--NumberParallel', type = int, default = 48, help = 'Specify the number of samples run in parallel.')
parser.add_argument('-s', '--SampleIDs', nargs = '+', metavar = '', choices = fm_obj.returnOptions('Samples'), help = 'Restrict analysis to the listed sampleIDs')
parser.add_argument('-c', '--Species', nargs = '+', metavar = '', choices = fm_obj.returnOptions('Species'), help = 'Restrict analysis to the following species: ' + ','.join(fm_obj.returnOptions('Species')))
parser.add_argument('-p', '--ProjectIDs', nargs = '+', metavar = '', choices = fm_obj.returnOptions('ProjectIDs'), help = 'Restrict analysis to a specific ProjectIDs: ' + ','.join(fm_obj.returnOptions('ProjectIDs')))
parser.add_argument('-r', '--Rerun', action = 'store_true', help = 'Default behavior is to not rerun alignment if already completed. Use this to force realignment')
args = parser.parse_args()

# Create FileManager object to keep track of filenames
fm_obj.setGenome(args.Genome)

# Create timer object to keep track of time
timer = Timer()

# This command identifies all the samples that will need to be run based upon user input and stores in self.samples and self.s_dt
fm_obj.setSamples(projectIDs = args.ProjectIDs, sampleIDs = args.SampleIDs, species = args.Species, rerun = args.Rerun)

# Download genome data necessary for analysis
timer.start('Downloading genome')		
#fm_obj.uploadData(fm_obj.localGenomeDir)
fm_obj.downloadData(fm_obj.localGenomeDir)
timer.stop()

# Create alignment worker object:
aw_obj = AW(args.Genome, fm_obj, check_size = True)

timer.start('  Parallel Downloading uBams files')
aw_obj.downloadReadData()
timer.stop()

timer.start('  Aligning Reads to created sorted Bamfiles')
aw_obj.alignData()
timer.stop()

timer.start('  Marking duplicates for bamfiles')
aw_obj.markDuplicates()
timer.stop()

timer.start('  Splitting reads based upon their alignment')
aw_obj.splitBamfiles()
timer.stop()

timer.start('  Calling haplotypes to create gvcf files')
aw_obj.createGVCF()
timer.stop()

timer.start('  Uploading and updating database')
aw_obj.uploadAndUpdateDatabase()
timer.stop()

