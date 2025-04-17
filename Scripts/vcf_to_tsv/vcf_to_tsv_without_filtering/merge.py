import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(description='Merge filtered and samples TSV files.')
    parser.add_argument('-f', '--filtered', required=True, help='Filtered TSV file')
    parser.add_argument('-s', '--samples', required=True, help='Samples TSV file')
    parser.add_argument('-o', '--output', required=True, help='Output TSV file')
    args = parser.parse_args()

    # Load data with explicit dtypes and low_memory=False
    filtered_df = pd.read_csv(
        args.filtered,
        sep='\t',
        dtype={
            '#CHROM': str,
            'CHROM': str,
            'POS': 'Int64',  # Handles missing values
            'ID': str,
            'REF': str,
            'ALT': str
        },
        low_memory=False
    )
    
    samples_df = pd.read_csv(
        args.samples,
        sep='\t',
        dtype={
            'CHROM': str,
            'POS': 'Int64',
            'ID': str,
            'REF': str,
            'ALT': str,
            'SAMPLES': str,
            'GENOTYPES': str
        }
    )

    # Handle chromosome column naming
    chrom_col = '#CHROM' if '#CHROM' in filtered_df.columns else 'CHROM'
    
    # Ensure consistent types
    filtered_df[chrom_col] = filtered_df[chrom_col].astype(str)
    samples_df['CHROM'] = samples_df['CHROM'].astype(str)
    samples_df['POS'] = pd.to_numeric(samples_df['POS'], errors='coerce').fillna(0).astype('int64')

    # Merge data
    merged_df = pd.merge(
        filtered_df,
        samples_df[['CHROM', 'POS', 'SAMPLES', 'GENOTYPES']],
        left_on=[chrom_col, 'POS'],
        right_on=['CHROM', 'POS'],
        how='inner'  # or 'left' depending on your needs
    )

    # Clean up
    merged_df.drop_duplicates(inplace=True)
    
    # Save output
    merged_df.to_csv(args.output, sep='\t', index=False)
    print(f"Successfully merged data. Output saved to {args.output}")

if __name__ == "__main__":
    main()