import os, pandas as pd
from sexloci.config import Config
from sexloci.labdata import LabData
from sexloci.tables import load_tables_labdata

al = load_tables_labdata(LabData.from_config(Config()), "Mzebra_GT3_NCBI").alignment
al = al.drop_duplicates("SampleID").set_index("SampleID")
for sp in ("gertrudae", "saulosi"):
    qc = pd.read_csv(os.path.expanduser(f"~/Temp/SexDeterminationAnalysis/{sp}_pedigree_extra_coverage/sample_qc.tsv"), sep="\t")
    qc["sheet_coverage"] = qc["sample"].map(al["Coverage"])
    qc["vcf_over_sheet"] = (qc["mean_depth"] / qc["sheet_coverage"]).round(2)
    print(f"\n{sp}: VCF depth is {qc.vcf_over_sheet.median():.2f} x the alignment-sheet coverage (median)")
    print(qc[["sample", "sex", "sheet_coverage", "mean_depth", "vcf_over_sheet", "missing_rate"]].sort_values("vcf_over_sheet").round(2).to_string(index=False))
