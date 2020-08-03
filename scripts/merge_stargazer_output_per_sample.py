import glob
import sys

import pandas as pd


OUTPUT_COLUMNS = [
    "chr",
    "target_gene",
    "status",
    "hap1_main",
    "hap2_main",
    "hap1_cand",
    "hap2_cand",
    "hap1_score",
    "hap2_score",
    "dip_score",
    "phenotype",
    "dip_sv",
    "hap1_sv",
    "hap2_sv",
    "ssr",
    "dip_cand",
    "hap1_main_core",
    "hap2_main_core",
    "hap1_main_tag",
    "hap2_main_tag",
    "hap1_af_mean_gene",
    "hap2_af_mean_gene",
    "hap1_af_mean_main",
    "hap2_af_mean_main",
]


TARGET_GENES = [
    "CACNA1S",
    "CFTR",
    "CYP1A1",
    "CYP1A2",
    "CYP1B1",
    "CYP2A6",
    "CYP2A7",
    "CYP2A13",
    "CYP2B6",
    "CYP2B7",
    "CYP2C8",
    "CYP2C9",
    "CYP2C19",
    "CYP2D6",
    "CYP2D7",
    "CYP2E1",
    "CYP2F1",
    "CYP2J2",
    "CYP2R1",
    "CYP2S1",
    "CYP2W1",
    "CYP3A4",
    "CYP3A5",
    "CYP3A7",
    "CYP3A43",
    "CYP4B1",
    "CYP26A1",
    "CYP4F2",
    "CYP19A1",
    "DPYD",
    "G6PD",
    "GSTM1",
    "GSTP1",
    "GSTT1",
    "IFNL3",
    "NAT1",
    "NAT2",
    "NUDT15",
    "POR",
    "RYR1",
    "SLC15A2",
    "SLC22A2",
    "SLCO1B1",
    "SLCO1B3",
    "SLCO2B1",
    "SULT1A1",
    "TBXAS1",
    "TPMT",
    "UGT1A1",
    "UGT1A4",
    "UGT2B7",
    "UGT2B15",
    "UGT2B17",
    "VKORC1",
]


sample = sys.argv[1]
result_files = glob.glob("*.stargazer-genotype.txt")

# Create a pandas dataframe for each found Stargazer result file
result_dfs = []
for result_file in result_files:
    gene = result_file.split(".")[1]
    df = pd.read_csv(result_file, sep="\t")
    df["gene"] = gene
    result_dfs.append(df)

# Concat all previously created dataframes into one, and find the processed samples
full_results_df = pd.concat(result_dfs)

# Separate Stargazer result per sample
sample_df = full_results_df[full_results_df["name"] == sample].reset_index(drop=True)

# Create empty row for genes with no Stargazer output
target_genes_df = pd.DataFrame(columns=TARGET_GENES).transpose().reset_index()
target_genes_df.rename(columns={"index": "target_gene"}, inplace=True)
sample_df = sample_df.merge(
    target_genes_df, how="right", left_on="gene", right_on="target_gene"
)
sample_df.sort_values("target_gene", inplace=True)

# Read intervals.json file and create dataframe from it
intervals_df = pd.read_json(sys.argv[2]).transpose().reset_index()

# Get genomic coordinates from intervals_df
sample_df = sample_df.merge(intervals_df, left_on="target_gene", right_on="index")

sample_df.to_csv(
    "{}.haplotypes.tsv".format(sample),
    sep="\t",
    na_rep=".",
    columns=OUTPUT_COLUMNS,
    index=False,
)
