import glob
import sys

import pandas as pd


OUTPUT_COLUMNS = [
    "chr",
    "gene",
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

# Read intervals.json file and create dataframe from it
intervals_df = pd.read_json(sys.argv[2]).transpose().reset_index()

# Get genomic coordinates from intervals_df
full_results_df = full_results_df.merge(intervals_df, left_on="gene", right_on="index")

# Separate Stargazer result per sample and create its respective output file
sample_df = full_results_df[full_results_df["name"] == sample].reset_index(
    drop=True
).sort_values("gene")
sample_df.to_csv(
    "{}.haplotypes.tsv".format(sample),
    sep="\t",
    columns=OUTPUT_COLUMNS,
    index=False,
)
