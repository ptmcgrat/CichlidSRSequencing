"""Centromeres of a published near-complete assembly (Table 3 of the Chromosome Research 2026 paper, 10577_2026_9819_MOESM2_ESM.xlsx) lifted into the reference.
The published assembly is aligned to the reference; each centromere interval (and centroid core) is mapped through the alignment, bounded by the nearest alignable flanks
when the satellite array itself does not align. Pure functions, tested on synthetic data."""
import numpy as np
import pandas as pd
from ..inversion_breakpoints import paf as P


def read_published(path, sheet, genome):
    """Rows of the sheet for one genome and chromosomes named chrN, with Size/Distance computed when Excel did not cache the formulas."""
    import openpyxl
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    rows = list(wb[sheet].iter_rows(values_only=True))
    hdr = [str(h).strip() if h is not None else '' for h in rows[0]]
    df = pd.DataFrame(rows[1:], columns=hdr)
    df = df[df['genome'] == genome].copy()
    m = df['chr'].astype(str).str.extract(r'^chr(\d+)$')[0]
    df = df[m.notna()].copy(); df['lg'] = 'LG' + m[m.notna()].astype(int).astype(str)
    for c in ('start', 'end', 'Centroid1', 'Centroid2'):
        df[c] = pd.to_numeric(df[c], errors='coerce')
    out = pd.DataFrame(dict(chr=df['chr'].values, lg=df['lg'].values, start=df.start.values, end=df.end.values, size=(df.end - df.start).values,
                            centroid1=df.Centroid1.values, centroid2=df.Centroid2.values, gap_in_centromere=df['Gap cen'].values, t2t=df['T2T'].values))
    return out.reset_index(drop=True)


def map_published_chr(contig_map, lg_to_chrom, lgs):
    """lg -> the published-assembly sequence that belongs to it: the sequence whose dominant reference chromosome is that LG's chromosome, with the most aligned bases.
    contig_map columns: contig, reference_chrom, aligned_bp, share."""
    out = {}
    for lg in lgs:
        target = lg_to_chrom.get(lg)
        c = contig_map[contig_map.reference_chrom == target].sort_values('aligned_bp', ascending=False)
        out[lg] = dict(seq=c.contig.iloc[0], aligned_bp=int(c.aligned_bp.iloc[0]), share=float(c.share.iloc[0]), n_candidates=int(len(c))) if len(c) else None
    return out


def lift_all(paf_df, tbl, mapping, lg_to_chrom, window=3_000_000):
    """One row per published centromere with its interval in reference coordinates (bounding interval between alignable flanks) and the centroid cores, where alignable."""
    rows = []
    for r in tbl.itertuples():
        m = mapping.get(r.lg)
        row = dict(lg=r.lg, published_chr=r.chr, published_start=int(r.start), published_end=int(r.end), centroid1=r.centroid1, centroid2=r.centroid2, t2t=r.t2t,
                   gap_in_centromere=r.gap_in_centromere, reference_chrom=lg_to_chrom.get(r.lg), assembly_sequence=None, sequence_share=np.nan,
                   gt3_start=np.nan, gt3_end=np.nan, gt3_width=np.nan, strand='', lift_status='', left_flank_distance=np.nan, right_flank_distance=np.nan,
                   centroid1_ref=np.nan, centroid2_ref=np.nan, left_flank_ref=np.nan, left_flank_strand='', right_flank_ref=np.nan, right_flank_strand='', bounded_start=np.nan, bounded_end=np.nan)
        if m is None:
            row['lift_status'] = 'no assembly sequence aligns to this chromosome'; rows.append(row); continue
        row.update(assembly_sequence=m['seq'], sequence_share=m['share'])
        # only alignments to the chromosome this LG should lift to: every centromere carries the same satellite family, so the nearest aligned block next to a satellite array
        # is often a repeat hit on ANOTHER chromosome ('flanks map to different chromosomes'); those are ignored
        own = paf_df[paf_df.tname == row['reference_chrom']]
        lf = P.lift_interval(own, m['seq'], int(r.start), int(r.end), max_scan=window)
        if lf is None:
            row['lift_status'] = 'no alignable flanks within the search window'
        else:
            row['lift_status'] = lf['status']
            row.update(left_flank_ref=lf.get('left_ref', np.nan), left_flank_strand=lf.get('left_strand', ''), right_flank_ref=lf.get('right_ref', np.nan), right_flank_strand=lf.get('right_strand', ''))
            if lf.get('tname', lf.get('tname_left')) == row['reference_chrom'] == lf.get('tname', lf.get('tname_right')) and 'left_ref' in lf:     # both flanks on the expected chromosome
                row.update(bounded_start=min(lf['left_ref'], lf['right_ref']), bounded_end=max(lf['left_ref'], lf['right_ref']))        # the region BETWEEN the flank landings, even if they differ in strand
            if lf['status'] == 'ok':                                    # an interval whose flanks lie on opposite strands straddles a breakpoint: no coordinates
                if lf['tname'] != row['reference_chrom']:
                    row['lift_status'] = f"lifts to {lf['tname']}, expected {row['reference_chrom']}"
                else:
                    row.update(gt3_start=lf['tstart'], gt3_end=lf['tend'], gt3_width=lf['tend'] - lf['tstart'], strand=lf['strand'],
                               left_flank_distance=lf['left_flank_distance'], right_flank_distance=lf['right_flank_distance'])
        for k in ('centroid1', 'centroid2'):
            v = getattr(r, k)
            if pd.notna(v):
                hit = [h for h in P.lift_point(paf_df, m['seq'], int(v)) if h[3] == 'aligned' and h[0] == row['reference_chrom']]
                if hit:
                    row[k + '_ref'] = hit[0][1]
        rows.append(row)
    return pd.DataFrame(rows)


def compare_with_scan(lifted, scan, agree_distance=1_000_000):
    """Add the satellite-scan call (reference genome) next to each lifted interval and say whether they agree (intervals overlap or lie within agree_distance)."""
    out = lifted.copy()
    for col in ('scan_start', 'scan_end', 'scan_confidence', 'scan_vs_published'):
        out[col] = np.nan if col != 'scan_confidence' and col != 'scan_vs_published' else ''
    s = scan.set_index('lg') if scan is not None and len(scan) else None
    for i, r in out.iterrows():
        if s is None or r.lg not in s.index or pd.isna(r.gt3_start) or pd.isna(s.loc[r.lg, 'start']):
            continue
        a0, a1 = float(s.loc[r.lg, 'start']), float(s.loc[r.lg, 'end'])
        gap = max(0.0, max(a0, r.gt3_start) - min(a1, r.gt3_end))
        out.loc[i, ['scan_start', 'scan_end']] = [a0, a1]
        out.loc[i, 'scan_confidence'] = s.loc[r.lg, 'confidence']
        out.loc[i, 'scan_vs_published'] = 'agrees' if gap <= agree_distance else f'disagrees: {gap / 1e6:.1f} Mb apart'
    return out
