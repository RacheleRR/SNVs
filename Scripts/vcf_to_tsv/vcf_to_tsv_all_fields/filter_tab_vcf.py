import pandas as pd
import argparse
import re

# Argument parser
parser = argparse.ArgumentParser(description='Filter TSV file.')
parser.add_argument('-i', '--input', required=True, help='Input TSV file')
parser.add_argument('-o', '--output', required=True, help='Output TSV file')
args = parser.parse_args()

input_tsv = args.input
output_tsv = args.output

info_columns = [
    "#CHROM", "POS", "ID", "REF", "ALT", "Allele", "pLI_gene_value", "IMPACT", "SYMBOL", 
    "Consequence", "MPC_score", "MVP_score", "gnomADg_AF", "gnomADe_AF", 
    "am_pathogenicity", "CADD_phred", "GERP___RS",
    "CLIN_SIG", "SOMATIC", "PHENO", "TSSDistance", "am_class", "am_pathogenicity", 
    "MPC", "Enformer_SAD", "Enformer_SAR", "pLI_gene_value", "Aloft_Confidence", 
    "Aloft_Fraction_transcripts_affected", "Aloft_pred", "BayesDel_addAF_pred", 
    "BayesDel_addAF_score", "BayesDel_noAF_pred", "BayesDel_noAF_score", "CADD_phred", 
    "CADD_raw", "ClinPred_score", "DANN_score", "DEOGEN2_pred", "DEOGEN2_score", 
    "ESM1b_pred", "ESM1b_score", "EVE_Class10_pred", "EVE_score", "Eigen-PC-raw_coding", 
    "Eigen-raw_coding", "FATHMM_pred", "FATHMM_score", "GERP++_RS", "GERP_91_mammals", 
    "GTEx_V8_eQTL_gene", "GTEx_V8_eQTL_tissue", "GTEx_V8_sQTL_gene", "GTEx_V8_sQTL_tissue", 
    "GenoCanyon_score", "LINSIGHT", "LIST-S2_pred", "LIST-S2_score", "LRT_pred", 
    "LRT_score", "M-CAP_pred", "M-CAP_score", "MPC_score", "MVP_score", "MetaLR_pred", 
    "MetaLR_score", "MetaRNN_pred", "MetaRNN_score", "MetaSVM_pred", "MetaSVM_score", 
    "MutFormer_rankscore", "MutFormer_score", "MutPred_score", "MutScore_score", 
    "MutationAssessor_pred", "MutationTaster_pred", "MutationTaster_score", 
    "PHACTboost_rankscore", "PHACTboost_score", "PROVEAN_pred", "PROVEAN_score", 
    "Polyphen2_HDIV_pred", "Polyphen2_HDIV_score", "Polyphen2_HVAR_pred", 
    "Polyphen2_HVAR_score", "PrimateAI_pred", "PrimateAI_score", "REVEL_score", 
    "Reliability_index", "SIFT_pred", "SIFT_score", "VARITY_ER_score", "VARITY_R_score", 
    "VEST4_score", "clinvar_clnsig", "clinvar_id", "fathmm-MKL_coding_group", 
    "fathmm-MKL_coding_pred", "fathmm-MKL_coding_rankscore", "fathmm-MKL_coding_score", 
    "fathmm-XF_coding_pred", "fathmm-XF_coding_rankscore", "fathmm-XF_coding_score", 
    "gMVP_score", "integrated_confidence_value", "integrated_fitCons_score", 
    "phastCons100way_vertebrate", "phastCons17way_primate", "phastCons470way_mammalian", 
    "phyloP100way_vertebrate", "phyloP17way_primate", "phyloP470way_mammalian", "SIFT"
]


# Load the TSV file into a pandas DataFrame
df = pd.read_csv(input_tsv, sep='\t', low_memory=False)

# Check if column names contain indices and remove them if they do
if any(re.search(r'\[\d+\]', col) for col in df.columns):
    df.columns = [re.sub(r'\[\d+\]', '', col).strip() for col in df.columns]

# Filter the DataFrame to keep only the desired columns that exist in the DataFrame
existing_columns = [col for col in info_columns if col in df.columns]
filtered_df = df[existing_columns].copy()

# Rename the '#CHROM' column to 'CHROM'
filtered_df.rename(columns={'#CHROM': 'CHROM'}, inplace=True)

# Keep only unique rows
filtered_df.drop_duplicates(inplace=True)

# Save the filtered DataFrame to a new TSV file
filtered_df.to_csv(output_tsv, sep='\t', index=False)

print("TSV data successfully filtered and saved to new TSV file.")