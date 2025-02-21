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
output_file="${input_file%.vcf}_acmg_fields_output.tsv"

# Define the header
header="CHROM\tPOS\tID\tREF\tALT\tacmg_classification_base\tacmg_criteria_base\tacmg_score_base"

# Run bcftools to extract only ACMG-related fields and add the header
{
    echo -e "$header"
    bcftools query -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%acmg_classification_base\t%acmg_criteria_base\t%acmg_score_base\n" "$input_file"
} > "$output_file"

if [ $? -ne 0 ]; then
    echo "Error: bcftools query failed."
    exit 1
fi

echo "ACMG-related fields extracted successfully into $output_file."
