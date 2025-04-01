#!/bin/bash

# Ensure the input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 -i input_file"
    exit 1
fi

input_file=""
while getopts "i:" opt; do
    case $opt in
        i) input_file=$OPTARG ;;
        *) echo "Usage: $0 -i input_file"; exit 1 ;;
    esac
done

# Define the output file
output_file="${input_file%.vcf}_all_fields_output.tsv"

# Run the bcftools command and check for errors
if ! bcftools +split-vep "$input_file" -c "CSQ" -HH -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%CSQ\n" -d -A tab -o "$output_file"; then
    echo "Error: bcftools command failed."
    exit 1
fi

echo "VCF file has been transformed to a tab-separated file with all fields."
