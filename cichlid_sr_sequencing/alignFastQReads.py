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
#fm_obj.downloadData(fm_obj.localGenomeDir)
timer.stop()

# Create alignment worker object:
aw_obj = AW(args.Genome, fm_obj)

timer.start('  Parallel Downloading uBams files')
#aw_obj.downloadReadData()
timer.stop()

timer.start('  Aligning Reads to created sorted Bamfiles')
#aw_obj.alignData()
timer.stop()

timer.start('  Marking duplicates for bamfiles')
#aw_obj.markDuplicates()
timer.stop()

timer.start('  Splitting reads based upon their alignment')
#aw_obj.splitBamfiles()
timer.stop()

timer.start('  Calling haplotypes to create gvcf files')
#aw_obj.createGVCF(parallel = True)
timer.stop()

processes = []
for sample in fm_obj.samples:
	fm_obj = FM(args.Genome, sample)
	#fm_obj.uploadData(fm_obj.localSampleBamDir)

	stats = aw_obj.calculateStats(sample)
	s_dt = fm_obj.s_dt
	read_length = s_dt[s_dt['SampleID'] == sample]['ReadLength'].values[0]/2
	reference_size = sum(pysam.FastaFile(fm_obj.localGenomeFile).lengths)
	coverage = stats['all'] * read_length / reference_size
	stats = {k:v/stats['all'] if k!= 'all' else v for k,v in stats.items()}
	sample_data = {'SampleID':sample, 'GenomeVersion': args.Genome, 'RunIDs':',,'.join(list(s_dt[s_dt['SampleID'] == sample].RunID)), 
			   'Coverage':coverage, 'TotalReads':stats['all'], 'UnmappedReads':stats['unmapped'], 'DiscordantReads':stats['discordant'], 'InversionReads':stats['inversion'],
			   'DuplicationReads':stats['duplication'], 'ClippedReads':stats['clipped'], 'ChimericReads':stats['chimeric']}

	output = subprocess.run(['conda', 'list'], capture_output = True)
	sample_data['bwa_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('bwa')][0]
	sample_data['gatk_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('gatk4')][0]
	sample_data['pysam_version'] = [x.split()[1] for x in output.stdout.decode('utf-8').split('\n') if x.startswith('pysam')][0]
	sample_data['BamSize'] = os.path.getsize(fm_obj.localBamFile)
	pdb.set_trace()
	if sample not in fm_obj.a_dt[fm_obj.a_dt.GenomeVersion == args.Genome].SampleID.to_list():
		fm_obj.a_dt.loc[len(fm_obj.a_dt)] = sample_data
	else:
		for key, value in sample_data.items():
			fm_obj.a_dt.loc[(fm_obj.a_dt.GenomeVersion == args.Genome) & (fm_obj.a_dt.SampleID == sample), key] = value
		pdb.set_trace()
	# Upload data and delete
	#subprocess.run(['rm','-rf', fm_obj.localSampleBamDir])
	#subprocess.run(['rm','-rf', fm_obj.localTempDir])
	fm_obj._setDatabase('AlignmentDatabase', fm_obj.a_dt)
	#print(' Finished with sample ' + sample + ': ' + str(datetime.datetime.now()))
	#print()

