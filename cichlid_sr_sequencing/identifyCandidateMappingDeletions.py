import subprocess,pysam,pdb
from helper_modules.file_manager import FileManager as FM
from helper_modules.ChimericData import ChimericRead
from collections import defaultdict

fm_obj = FM(genome_version = 'Mzebra_GT3_NCBI')

yh_samples = ['YH_003_m','YH_004_m','YH_005_m','YH_008_m','YH_009_m','YH_011_m','YH_012_m','YH_014_m','YH_015_m']
yh_f_samples = ['Kocher_YH2f','YH_006_f','YH_007_f','YH_010_f','YH_013_f']
mc_samples = ['MC-5B6B-f','MC-2G8G-f','MC-3P6P-f','MC-4O5O-f','MC-2B4B-f']
f1_f_samples = ['MCYHF1_016_f','MCYHF1_015_f','MCYHF1_014_f','MCYHF1_017_f']
f1_m_samples = ['MCYHF1_012_m','MCYHF1_010_m','MCYHF1_013_m','MCYHF1_011_m']

discoveryChimeras = defaultdict(int)
yh_data = {}
exclusionChimeras = defaultdict(int)
femaleChimeras = defaultdict(int)
f1femaleChimeras = defaultdict(int)
f1maleChimeras = defaultdict(int)

for sampleID in yh_samples:
    yh_data[sampleID] = defaultdict(int)
    fm_obj.createSampleFiles(sampleID)
    fm_obj.downloadData(fm_obj.localChimericBamFile)
    fm_obj.downloadData(fm_obj.localChimericBamFile.replace('.bam','.bai'))

    bamObj = pysam.AlignmentFile(fm_obj.localChimericBamFile)
    for read in bamObj.fetch('NC_135176.1'):
        if not read.is_secondary and read.mapq > 10:
            try:
                SA = read.get_tag('SA')
            except KeyError:
                continue
            if read.get_tag('SA').split(',')[0] != 'NC_135176.1':
                continue
            read.set_tag('SA',SA.replace(SA.split(',')[0], str(9)))
            try:
                newRead = ChimericRead(read)
            except TypeError:
                continue
            discoveryChimeras[newRead.data] += 1
            yh_data[sampleID][newRead.data] += 1

for sampleID in mc_samples:
    fm_obj.createSampleFiles(sampleID)
    fm_obj.downloadData(fm_obj.localChimericBamFile)
    fm_obj.downloadData(fm_obj.localChimericBamFile.replace('.bam','.bai'))

    bamObj = pysam.AlignmentFile(fm_obj.localChimericBamFile)
    for read in bamObj.fetch('NC_135176.1'):
        if not read.is_secondary and read.mapq > 10:
            try:
                SA = read.get_tag('SA')
            except KeyError:
                continue
            if read.get_tag('SA').split(',')[0] != 'NC_135176.1':
                continue
            read.set_tag('SA',SA.replace(SA.split(',')[0], str(9)))
            try:
                newRead = ChimericRead(read)
            except TypeError:
                continue
            exclusionChimeras[newRead.data] += 1

for sampleID in yh_f_samples:
    fm_obj.createSampleFiles(sampleID)
    fm_obj.downloadData(fm_obj.localChimericBamFile)
    fm_obj.downloadData(fm_obj.localChimericBamFile.replace('.bam','.bai'))

    bamObj = pysam.AlignmentFile(fm_obj.localChimericBamFile)
    for read in bamObj.fetch('NC_135176.1'):
        if not read.is_secondary and read.mapq > 10:
            try:
                SA = read.get_tag('SA')
            except KeyError:
                continue
            if read.get_tag('SA').split(',')[0] != 'NC_135176.1':
                continue
            read.set_tag('SA',SA.replace(SA.split(',')[0], str(9)))
            try:
                newRead = ChimericRead(read)
            except TypeError:
                continue
            femaleChimeras[newRead.data] += 1

for sampleID in f1_f_samples:
    fm_obj.createSampleFiles(sampleID)
    fm_obj.downloadData(fm_obj.localChimericBamFile)
    fm_obj.downloadData(fm_obj.localChimericBamFile.replace('.bam','.bai'))

    bamObj = pysam.AlignmentFile(fm_obj.localChimericBamFile)
    for read in bamObj.fetch('NC_135176.1'):
        if not read.is_secondary and read.mapq > 10:
            try:
                SA = read.get_tag('SA')
            except KeyError:
                continue
            if read.get_tag('SA').split(',')[0] != 'NC_135176.1':
                continue
            read.set_tag('SA',SA.replace(SA.split(',')[0], str(9)))
            try:
                newRead = ChimericRead(read)
            except TypeError:
                continue
            f1femaleChimeras[newRead.data] += 1

for sampleID in f1_m_samples:
    fm_obj.createSampleFiles(sampleID)
    fm_obj.downloadData(fm_obj.localChimericBamFile)
    fm_obj.downloadData(fm_obj.localChimericBamFile.replace('.bam','.bai'))

    bamObj = pysam.AlignmentFile(fm_obj.localChimericBamFile)
    for read in bamObj.fetch('NC_135176.1'):
        if not read.is_secondary and read.mapq > 10:
            try:
                SA = read.get_tag('SA')
            except KeyError:
                continue
            if read.get_tag('SA').split(',')[0] != 'NC_135176.1':
                continue
            read.set_tag('SA',SA.replace(SA.split(',')[0], str(9)))
            try:
                newRead = ChimericRead(read)
            except TypeError:
                continue
            f1maleChimeras[newRead.data] += 1

print('Chrom,Start,Stop,Length,AllHits,YHFemHits,F1MaleHits,F1FemHits' + ','.join(yh_samples))
for info,hits in sorted(discoveryChimeras.items(), key=lambda item: item[0][1]):
    if hits < 10 or info[7] != 'del':
        continue
    length = info[4] - info[1]
    if length < 40 or length > 300:
        continue
    if exclusionChimeras[info] > 0:
        continue
    if info[1] < 23700000 or info[1] > 25700000:
        continue
    outdata = ['NC_135176.1',info[1],info[4],length,hits,femaleChimeras[info],f1maleChimeras[info],f1femaleChimeras[info]] + [yh_data[x][info] for x in yh_samples]
    print(','.join([str(x) for x in outdata]))