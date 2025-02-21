#!/bin/bash
#SBATCH --job-name=SNV_Rachele    # Job name
#SBATCH --output=snv_%j.out       # Output file name
#SBATCH --error=snv_%j.err        # Error file name
#SBATCH --ntasks=1                # Number of tasks (1 task in this case)
#SBATCH --nodelist=macgyver  
#SBATCH --cpus-per-task=8         # Number of CPUs per task (8 CPUs)
#SBATCH --mem=200G                # Memory allocation
#SBATCH --time=INFINITE           # Unlimited execution time
#SBATCH --partition=core          # Partition to use (default "core" partition)

# Create necessary directories for each step of the pipeline
mkdir -p step1_index_sort step2_merge step3_normalize step4_VR step5_VQSR step6_filtered step7_split step8_VAR_Filter step9_filtered step10_merge step11_clean step12_normalize step13_rs step14_VEP stepX_log

# STEP 1: Index and sort
# Define the base NAS folder where the input data is stored
NAS_FOLDER="/media/NAS_KREBS_1/WGS/DRAGEN_pipeline_results"

# Input file containing the list of sample names
SAMPLES_FILE="wgs.list"

# Log files for tracking progress and failures
FAILED_SAMPLES_FILE="failed_samples.txt"
PROCESSED_LOG_FILE="processed_files.log"

# Clear previous logs to start fresh
> "$FAILED_SAMPLES_FILE"
> "$PROCESSED_LOG_FILE"

# Loop through each sample in the sample list
while IFS= read -r sample; do
    # Construct the folder path for the current sample
    folder="$NAS_FOLDER/${sample}/dragen/"
    
    count=0
    vcf_found=false

    # Use process substitution to avoid subshell issues
    while IFS= read -r -d '' vcf; do
        vcf_found=true
        ((count++))
        
        # Extract base name and define sorted file name
        base_name=$(basename "$vcf" .vcf.gz)
        sorted_file="${base_name}_sort.vcf.gz"

        # Sort and index the VCF file
        if bcftools sort \
            "$vcf" \
            -O z \
            -o "step1_index_sort/$sorted_file" >> stepX_log/output.log 2>> stepX_log/error.log; then
            bcftools index \
                -t "step1_index_sort/$sorted_file" \
                -o "step1_index_sort/${sorted_file}.tbi" \
                --threads 8 >> stepX_log/output.log 2>> stepX_log/error.log || {
                echo "Failed to index: $sorted_file" >> "$PROCESSED_LOG_FILE"
                continue
            }
        else
            echo "Failed to sort: $vcf" >> "$PROCESSED_LOG_FILE"
            continue
        fi
    done < <(find "$folder" -name "*hard-filtered.vcf.gz" -print0)

    # Log results
    if ! $vcf_found; then
        echo "$sample" >> "$FAILED_SAMPLES_FILE"
    else
        echo "Processed $count files for sample: $sample" >> "$PROCESSED_LOG_FILE"
    fi
done < "$SAMPLES_FILE"

echo "Processing complete. Samples with no VCF files logged in $FAILED_SAMPLES_FILE."

# STEP 2: Merge all VCF files together
echo "Merging all VCF files"
# Merge all sorted VCF files into a single file
bcftools merge \
    step1_index_sort/*sorted_*.vcf.gz \
    -Oz \
    -o step2_merge/merged.vcf.gz \
    --threads 8 >> stepX_log/output.log 2>> stepX_log/error.log
echo "All files merged successfully."

# STEP 3: Normalize
# Normalize the merged VCF file and save the output to the step3_normalize directory
echo "Normalizing merged VCF file"
bcftools norm \
    -m -any \
    -Oz \
    -o step3_normalize/normalized_merged.vcf.gz \
    step2_merge/merged.vcf.gz \
    --threads 8 >> stepX_log/output.log 2>> stepX_log/error.log

# STEP 4: Variant Recalibrator
echo "Processing normalized merged VCF file..."
# Index the normalized merged VCF file
if ! bcftools index \
    -t step3_normalize/normalized_merged.vcf.gz \
    -o step3_normalize/normalized_merged.vcf.gz.tbi \
    --threads 8 >> stepX_log/output.log 2>> stepX_log/error.log; then
    echo "Error indexing merged file. Exiting."
    exit 1
fi

# Run GATK VariantRecalibrator
# This step uses GATK (Genome Analysis Toolkit) to recalibrate the variant quality scores.
# It takes the normalized VCF file and several known variant resources as input.
# The output includes a recalibration file, a tranches file, and an R script for plotting.
if ! gatk VariantRecalibrator \
    -R Homo_sapiens_assembly38.fasta \
    -V step3_normalize/normalized_merged.vcf.gz \
    --resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
    --resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf.gz \
    --resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf.gz \
    --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 Homo_sapiens_assembly38.dbsnp138.vcf \
    -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
    -mode BOTH \
    -O step4_VR/normalized_merged.recal \
    --tranches-file step4_VR/normalized_merged.tranches \
    --rscript-file step4_VR/normalized_merged.plots.R \
    --dont-run-rscript >> stepX_log/output.log 2>> stepX_log/error.log; then
    echo "Error processing merged file. Exiting."
    exit 1
fi
echo "Finished processing merged file."

# STEP 5: Apply VQSR
# Apply VQSR (Variant Quality Score Recalibration) using GATK
# This step applies the recalibration to the VCF file using the recalibration file and tranches file generated in the previous step.
if ! gatk --java-options "-Xmx15g -Xms15g" ApplyVQSR \
    -R Homo_sapiens_assembly38.fasta \
    -V step3_normalize/normalized_merged.vcf.gz \
    -O step5_VQSR/recalibrated.vcf.gz \
    --truth-sensitivity-filter-level 99.0 \
    --tranches-file step4_VR/normalized_merged.tranches \
    --recal-file step4_VR/normalized_merged.recal \
    -mode BOTH >> stepX_log/output.log 2>> stepX_log/error.log; then
    echo "Error applying VQSR. Exiting."
    exit 1
fi
echo "Finished applying VQSR."

# STEP 6: Filter out variants that didn't pass
# Filter out variants that didn't pass using GATK SelectVariants
# This step removes variants that did not pass the VQSR filtering.
if ! gatk SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V step5_VQSR/recalibrated.vcf.gz \
    -O step6_filtered/filtered.vcf.gz \
    --exclude-filtered >> stepX_log/output.log 2>> stepX_log/error.log; then
    echo "Error filtering variants. Exiting."
    exit 1
fi
echo "Finished filtering variants."

# STEP 7: Split file into INDELs/SNVs
# Split the filtered VCF file into separate files for SNPs and INDELs
# This step uses GATK SelectVariants to create separate VCF files for SNPs and INDELs.
gatk --java-options '-Xmx64g' SelectVariants \
    --variant step6_filtered/filtered.vcf.gz \
    --select-type-to-include SNP \
    --output step7_split/filtered_SNV.vcf.gz >> stepX_log/output.log 2>> stepX_log/error.log
gatk --java-options '-Xmx64g' SelectVariants \
    --variant step6_filtered/filtered.vcf.gz \
    --select-type-to-include INDEL \
    --output step7_split/filtered_INDEL.vcf.gz >> stepX_log/output.log 2>> stepX_log/error.log

# STEP 8: Variant filtration
# Apply variant filtration to SNPs
# This step uses GATK VariantFiltration to apply several filters to the SNPs.
# Filters include QD, QUAL, SOR, FS, MQ, MQRankSum, and ReadPosRankSum.
gatk --java-options '-Xmx64g' VariantFiltration \
    --variant step7_split/filtered_SNV.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "SOR > 3.0" --filter-name "SOR3" \
    --filter-expression "FS > 60.0" --filter-name "FS60" \
    --filter-expression "MQ < 40.0" --filter-name "MQ40" \
    --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    --output step8_VAR_Filter/filtered_SNV.filtered.vcf.gz >> stepX_log/output.log 2>> stepX_log/error.log

# Apply variant filtration to INDELs
# This step uses GATK VariantFiltration to apply several filters to the INDELs.
# Filters include QD, QUAL, FS, and ReadPosRankSum.
gatk --java-options '-Xmx64g' VariantFiltration \
    --variant step7_split/filtered_INDEL.vcf.gz \
    --filter-expression "QD < 2.0" --filter-name "QD2" \
    --filter-expression "QUAL < 30.0" --filter-name "QUAL30" \
    --filter-expression "FS > 200.0" --filter-name "FS200" \
    --filter-expression "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" \
    --output step8_VAR_Filter/filtered_INDEL.filtered.vcf.gz >> stepX_log/output.log 2>> stepX_log/error.log

# STEP 9: Filter out flagged variants
# Filter out flagged variants using GATK SelectVariants
# This step removes variants that were flagged by the previous filtration step.
gatk --java-options '-Xmx64g' SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V step8_VAR_Filter/filtered_SNV.filtered.vcf.gz \
    -O step9_filtered/filtered_SNV.filtered.PASS.vcf.gz \
    --exclude-filtered >> stepX_log/output.log 2>> stepX_log/error.log

gatk --java-options '-Xmx64g' SelectVariants \
    -R Homo_sapiens_assembly38.fasta \
    -V step8_VAR_Filter/filtered_INDEL.filtered.vcf.gz \
    -O step9_filtered/filtered_INDEL.filtered.PASS.vcf.gz \
    --exclude-filtered >> stepX_log/output.log 2>> stepX_log/error.log

# STEP 10: Merge them
# Merge the filtered SNP and INDEL VCF files
# This step uses GATK MergeVcfs to combine the filtered SNP and INDEL VCF files into a single VCF file.
gatk --java-options '-Xmx64g' MergeVcfs \
    -I step9_filtered/filtered_INDEL.filtered.PASS.vcf.gz \
    -I step9_filtered/filtered_SNV.filtered.PASS.vcf.gz \
    -O step10_merge/merged.filtered.PASS.vcf.gz >> stepX_log/output.log 2>> stepX_log/error.log

# STEP 11: Apply sed
# Remove 'chr' prefix from chromosome names
# This step uses sed to remove the 'chr' prefix from chromosome names in the VCF file.
zcat step10_merge/merged.filtered.PASS.vcf.gz | sed '/^#/! s/^chr//' > step11_clean/clean.vcf

# STEP 12: Normalize
# Normalize the merged VCF file
# This step uses bcftools norm to normalize the VCF file, ensuring consistent representation of variants.
bcftools norm \
    --threads 8 \
    -m -any \
    -Oz \
    -o step12_normalize/normalized_merged.filtered.PASS.vcf.gz \
    step10_merge/merged.filtered.PASS.vcf.gz \
    --threads 8 >> stepX_log/output.log 2>> stepX_log/error.log

# STEP 13: Annotate with rs IDs
# Annotate the normalized VCF files with rs IDs
# This step uses bcftools annotate to add rs IDs to the VCF file based on a reference VCF file.
bcftools index \
    --threads 8 \
    -t step12_normalize/normalized_merged.filtered.PASS.vcf.gz \
    -o step12_normalize/normalized_merged.filtered.PASS.vcf.gz.tbi \
    --threads 8 >> stepX_log/output.log 2>> stepX_log.error.log

bcftools annotate \
    -a 00-ALL.vcf.gz \
    -c ID \
    -o step13_rs/rs_normalized_merged.filtered.PASS.vcf.gz \
    step12_normalize/normalized_merged.filtered.PASS.vcf.gz \
    --threads 8 >> stepX_log/output.log 2>> stepX_log.error.log

# STEP 14: VEP (Variant Effect Predictor)
# ...existing code for VEP...


