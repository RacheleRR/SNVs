## SNVS STEPS

### OVERVIEW
This document outlines the step-by-step workflow for processing, annotating, and analyzing SNVs. It includes variant preparation, quality control, annotation, statistical analysis, and enrichment/network analysis.

### OBJECTIVES
- Process and filter SNVs from raw VCF files
- Apply quality control and variant recalibration
- Annotate variants with known databases
- Convert VCF to TSV for downstream analysis
- Perform statistical and enrichment analyses

## STEPS 

### Pipeline Steps

#### Step 1: Index and Sort
- Reads sample names from `wgs.list`.
- Finds `*hard-filtered.vcf.gz` files in the NAS directory.
- Sorts and indexes VCF files using `bcftools`.
- Logs failed and processed samples.

#### Step 2: Merge VCF Files
- Merges all sorted VCF files into `step2_merge/merged.vcf.gz`.

#### Step 3: Normalize
- Normalizes the merged VCF file using `bcftools norm`.

#### Step 4: Variant Recalibration
- Runs GATK `VariantRecalibrator` using known variant resources.

#### Step 5: Apply Variant Quality Score Recalibration (VQSR)
- Applies recalibration to variants using GATK `ApplyVQSR`.

#### Step 6: Filter Variants
- Removes low-quality variants using GATK `SelectVariants`.

#### Step 7: Split Variants
- Separates SNPs and INDELs into different files using GATK `SelectVariants`.

#### Step 8: Variant Filtration
- Applies quality filters to SNPs and INDELs using GATK `VariantFiltration`.

#### Step 9: Remove Filtered Variants
- Excludes flagged variants from SNP and INDEL files.

#### Step 10: Merge SNP and INDEL Files
- Merges filtered SNP and INDEL files into one VCF.

#### Step 11: Clean Chromosome Names
- Removes `chr` prefix from chromosome names using `sed`.

#### Step 12: Normalize
- Ensures consistent representation of variants using `bcftools norm`.

#### Step 13: Annotate with rs IDs
- Adds reference SNP (rs) IDs using `bcftools annotate`.

#### Step 14: Variant Effect Prediction (VEP)
- Adds prediction scores etc. to each variant.

#### Step 15: GeneBe
- Adds ACMG prediction scores and ClinGene scores to variants.

#### Step 16: Filtering
- Filters out variants that don't meet the criteria (pathogenic, non-pathogenic, PTV).

#### Step 17: VCF to TSV
- Converts the VCF file to TSV file for easier visualization and R studio utilization.

#### Step 18: Statistic Applied
- Applies Wilcoxon, Fisher, and t-test to see significance.

#### Step 19: Enrichment Analysis
- Performs enrichment analysis (Enrichr, Enrichment Map, etc.).

#### Step 20: Network Analysis
- Conducts network analysis (GeneMANIA, STRING, HumanBase).






## COMMANDS RUN 

### VCF Preparation (Step 1 - 13) (on server)
```bash
sbatch SNV_all_modifed.sh
```

### VEP Annotation (on server)
```bash
sbatch VEP_long.sh
sbatch VEP_short.sh
```

### GeneBe (on my local computer)
```bash
java -jar GeneBeClient-0.1.0-a.4.jar vcf annotate --input-vcf pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered.vcf.gz --output-vcf try.vcf --genome hg38 --api-key ak-lIYIgqyez5moX9NUkADnekytk --username rubiu.1997351@studenti.uniroma1.it
```

### VCF to TSV
```bash
transform_VCF.sh -i /media/rachele/NAS_InfoGene/Rachele_Rubiu/final_vcf_to_csv/pred_PTV_VEP_filtered.vcf.gz
transform_VCF.sh -i /media/rachele/NAS_InfoGene/Rachele_Rubiu/final_vcf_to_csv/pred_PTV_pLI_VEP_filtered.vcf.gz
transform_VCF.sh -i /media/rachele/NAS_InfoGene/Rachele_Rubiu/final_vcf_to_csv/pred_patho_MPC_ALPHAMISSENSE_VEP_filtered.vcf.gz
transform_VCF.sh -i /media/rachele/NAS_InfoGene/Rachele_Rubiu/final_vcf_to_csv/pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered.vcf.gz
transform_VCF.sh -i /media/rachele/NAS_InfoGene/Rachele_Rubiu/final_vcf_to_csv/pred_NON_patho_MPC_VEP_filtered.vcf.gz
transform_VCF.sh -i /media/rachele/NAS_InfoGene/Rachele_Rubiu/final_vcf_to_csv/pred_NON_patho_MPC_pLI_VEP_filtered.vcf.gz
```

### Statistics
- `#ALLA STATISTICS.r`
- `#wilxoxonallrows.r`
- `plots_snvs_statitics.r`

### Enrichment and Network
- `#GET GET NAMES.r`
- `#GET GENE NAMES ASSOCIATES TO SCHIZO ETC.r`