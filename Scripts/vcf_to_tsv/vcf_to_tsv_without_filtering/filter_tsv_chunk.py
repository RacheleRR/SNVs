import pandas as pd
import argparse
import re
import os
from tqdm import tqdm  # For progress bars

def main():
    parser = argparse.ArgumentParser(description='Filter large TSV file in chunks.')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # Define columns exactly once
    INFO_COLUMNS = [
    "#CHROM", "POS", "ID", "REF", "ALT", "Allele", "pLI_gene_value", "IMPACT", "SYMBOL", 
    "Consequence", "MPC_score", "MVP_score", "gnomADg_AF", "gnomADe_AF", 
    "am_pathogenicity", "CADD_phred", "GERP___RS",
    "CLIN_SIG", "SOMATIC", "PHENO", "TSSDistance", "am_class", "am_pathogenicity", 
    "MPC", "pLI_gene_value",  "CADD_phred", 
    "CADD_raw", "ClinPred_score", "DANN_score","GERP_91_mammals", "MPC_score", "MVP_score", "MetaLR_pred", 
    "MetaLR_score", "MetaRNN_pred", "MetaRNN_score", "MetaSVM_pred", "MetaSVM_score", 
    "MutFormer_rankscore", "MutFormer_score", "MutPred_score", "MutScore_score", 
    "MutationAssessor_pred", "MutationTaster_pred", "MutationTaster_score",  "SIFT_pred", "SIFT_score", "VARITY_ER_score", "VARITY_R_score", 
    "VEST4_score", "clinvar_clnsig", "clinvar_id",
    "gMVP_score", "SIFT"
    ]

    # Clean existing output
    if os.path.exists(args.output):
        os.remove(args.output)

    # Get total rows for progress bar (optional)
    total_rows = sum(1 for _ in open(args.input)) - 1  # Subtract header
    
    # Process in chunks
    seen_hashes = set()
    first_chunk = True
    
    for chunk in tqdm(
        pd.read_csv(args.input, sep='\t', chunksize=100000, low_memory=False),
        total=total_rows//100000 + 1,
        desc="Processing"
    ):
        # Clean columns
        chunk.columns = [re.sub(r'\[\d+\]', '', col).strip() for col in chunk.columns]
        
        # Standardize CHROM column
        # if '#CHROM' in chunk.columns:
        #     chunk = chunk.rename(columns={'#CHROM': 'CHROM'})
        
        # Select existing columns
        cols_to_keep = [col for col in INFO_COLUMNS if col in chunk.columns]
        filtered = chunk[cols_to_keep].copy()
        
        # Deduplicate using row hashes
        row_hashes = filtered.apply(lambda x: hash(tuple(x)), axis=1)
        filtered = filtered[~row_hashes.isin(seen_hashes)]
        seen_hashes.update(row_hashes.unique())
        
        # Save to disk
        filtered.to_csv(
            args.output,
            sep='\t',
            index=False,
            mode='a',
            header=first_chunk
        )
        first_chunk = False

    print(f"\n✅ Success! Output saved to {args.output}")

if __name__ == "__main__":
    main()