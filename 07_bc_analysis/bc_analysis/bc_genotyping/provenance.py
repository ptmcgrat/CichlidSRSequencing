"""Pre-flight report for the chosen genome: chromosome table, alignment status, missing BAM/GVCF files, pedigree."""
from cichlid_common.common import sh


def main(proj):
    cfg, out = proj.cfg, []
    out.append(f'## genome {proj.genome}: {proj.genome_fa}')
    out.append('## chromosome table (check LG labels!)\n' + proj.chrom_table.to_string(index=False))
    dt = proj.sample_table()
    chk = dt[dt.SampleID.isin(['MCYHF1-002', 'MCYHF1-003', 'MCYHF1_010_m', 'MCYHF1_014_f'])][['SampleID', 'Parents']]
    out.append('## F1 pedigree from SampleDatabase\n' + chk.to_string(index=False))
    bc = proj.bc_table()
    out.append('## BC brood / sire / dam\n' + bc.groupby(['BroodID', 'Sire', 'Dam']).size().to_string())
    need_deep = [cfg['yh_sire'], cfg['f1_dam']] + cfg['bc_dams'] + cfg['mc_panel'] + list(cfg['f1_sire_names'].values())
    need_deep += [cfg['extra_f1_family']['yh_sire']] + cfg['extra_f1_family']['f1']
    need_deep = list(dict.fromkeys(need_deep))
    proj.fm.readAlignmentDatabase()
    adb = proj.fm.alignment_dt
    have = set(adb[adb.GenomeVersion == proj.genome].SampleID)
    allneed = need_deep + bc.SampleID.tolist()
    out.append(f'## alignment database ({proj.genome}): {len([s for s in allneed if s in have])}/{len(allneed)} needed samples aligned\n'
               'NOT in alignment DB: ' + (', '.join(s for s in allneed if s not in have) or 'none'))
    miss = []
    for s in allneed:
        files = proj.fm.returnCloudFiles(proj.sample_dir(s))
        if s + '.all.bam' not in files:
            miss.append(s + '.all.bam')
        if s in need_deep and s + '.g.vcf.gz' not in files:
            miss.append(s + '.g.vcf.gz (will be called from the BAM)')
    out.append('## missing in cloud\n' + (', '.join(miss) or 'none'))
    out.append('## tools\n' + '\n'.join(f'{t}: ' + sh(f'{t} --version 2>&1 | head -1', check=False).stdout.strip()
                                        for t in ('samtools', 'bcftools', 'gatk', 'bedtools', 'rclone')))
    txt = '\n\n'.join(out)
    open(proj.dirs['qc'] + 'provenance.txt', 'w').write(txt)
    print(txt)
