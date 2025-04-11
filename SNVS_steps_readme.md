🧬 SNVS STEPS — A Guide to Small but Mighty Variants
🌟 Overview

Welcome to the SNV (Single Nucleotide Variant) workflow! This guide walks you through the full journey — from raw VCF files to insightful statistical, enrichment, and network analyses. You'll prep, clean, annotate, and analyze your SNVs with precision and care. 🧹🔬📊
🎯 Objectives

Here's what this pipeline helps you achieve:

    📥 Process & filter SNVs from raw VCF files

    ✅ Apply quality control + variant recalibration

    🏷️ Annotate variants using public databases

    🔄 Convert VCF → TSV for downstream fun

    📈 Perform statistical + enrichment analyses

🧪 Pipeline Steps
🗂 Step 1: Index & Sort

    Reads sample names from wgs.list

    Locates *hard-filtered.vcf.gz files

    Sorts/indexes using bcftools

    Logs successes + errors

🧬 Step 2: Merge VCF Files

    Combines all sorted VCFs into merged.vcf.gz

🧹 Step 3: Normalize

    Uses bcftools norm for structural consistency

📊 Step 4: Variant Recalibration

    Runs GATK VariantRecalibrator with known resources

🔧 Step 5: Apply VQSR

    Applies variant quality score recalibration using GATK

🚫 Step 6: Filter Variants

    Removes low-quality variants (SelectVariants)

🔍 Step 7: Split Variants

    Separates SNPs vs INDELs

🎛️ Step 8: Variant Filtration

    Applies quality filters to each file type

🧼 Step 9: Remove Filtered Variants

    Excludes flagged variants from files

🧩 Step 10: Merge SNP + INDEL

    Combines cleaned SNP + INDEL into one VCF

🧹 Step 11: Clean Chromosome Names

    Removes chr prefix via sed

🧬 Step 12: Normalize Again

    One last normalization sweep using bcftools

🧾 Step 13: Add rs IDs

    Annotates with reference SNP IDs

🔮 Step 14: Variant Effect Prediction (VEP)

    Adds functional prediction scores

🧠 Step 15: GeneBe

    Annotates variants with ACMG + ClinGene scores

✂️ Step 16: Filtering

    Filters for relevant variant types (e.g., pathogenic, PTVs)

📄 Step 17: VCF → TSV

    Converts to TSV for easier use in R Studio

🧠 Step 18: Statistical Analysis

    Runs Wilcoxon, Fisher, and t-tests to assess significance

🧬 Step 19: Enrichment Analysis

    Uses tools like Enrichr + Enrichment Map

🌐 Step 20: Network Analysis

    Explores interactions with GeneMANIA, STRING, HumanBase

💻 Commands
🛠 VCF Preparation (Steps 1–13)

sbatch SNV_all.sh

🔮 VEP Annotation

sbatch VEP_long.sh
sbatch VEP_short.sh

🧠 GeneBe Annotation (Run Locally)

java -jar GeneBeClient-0.1.0-a.4.jar vcf annotate \
--input-vcf ___.vcf.gz \
--output-vcf try.vcf \
--genome hg38 \
--api-key [KEY] \
--username [USERNAME]

📄 VCF to TSV Conversion

transform_VCF.sh -i path/to/your/filtered_file.vcf.gz

Repeat the above for each relevant file:

    pred_PTV_VEP_filtered.vcf.gz

    pred_PTV_pLI_VEP_filtered.vcf.gz

    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered.vcf.gz

    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered.vcf.gz

    pred_NON_patho_MPC_VEP_filtered.vcf.gz

    pred_NON_patho_MPC_pLI_VEP_filtered.vcf.gz

📊 Statistical Analysis Scripts

    fisher_stat_scz_detail.r

    statistics_fisher_new_correct.r

    Statistics_fisher_scz_3_groups.r

    statistics_wilcoxon_SCZ_3_groups.r

    statistics_wilxocoxon_right.r

    statstics_wilcoxon_SCZ.r

✨ Enrichment & Network Scripts

    GENE_lists.r

🧡 Final Notes

This pipeline is built with love, logic, and lots of bcftools. Whether you're identifying critical SNVs or just getting started in variant analysis — this guide has your back. Happy analyzing! 🧬✨
