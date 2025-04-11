---

# 🧬 SNVs: From VCF to Insights  

Welcome to **SNVs**! This project guides you through processing, annotating, and analyzing **Single Nucleotide Variants (SNVs)** using a structured pipeline that includes quality control, annotation, statistical tests, and biological interpretation through enrichment and network analysis.

---

## 🌟 Overview  

This repository walks through how to:

- Prepare and process SNV data from raw VCFs  
- Apply quality control and variant recalibration  
- Annotate variants using public databases  
- Convert VCFs into TSVs for analysis  
- Run statistical, enrichment, and network analyses  

---

## 🎯 Objectives  

✔️ Filter and prepare SNVs from raw variant data  
✔️ Apply GATK-based quality control and recalibration  
✔️ Annotate SNVs using rsIDs, VEP, and GeneBe  
✔️ Convert VCF files to TSV format for visualization and R analysis  
✔️ Run statistical and enrichment analyses on variant sets  

---

## 🛠️ Prerequisites  

Make sure you have the following tools installed and available:

- `bcftools` 🔧  
- `GATK` (v4+) 🧬  
- `VEP` (Variant Effect Predictor) 📋  
- `GeneBeClient` 🧠  
- `R` + necessary analysis scripts 📊  
- `transform_VCF.sh` script for VCF → TSV conversion 🔄  

---
![ChatGPT Image Apr 11, 2025, 05_46_40 PM](https://github.com/user-attachments/assets/c4abfe56-6aca-4d4c-9ef5-e9de605a1a79)

## 🔍 SNV Processing Steps  

### 1️⃣ Index & Sort  

- Read sample names from `wgs.list`  
- Locate and sort `*hard-filtered.vcf.gz` files using `bcftools`  
- Index and log processed/failed samples  

### 2️⃣ Merge VCF Files  

- Combine all VCFs into `step2_merge/merged.vcf.gz`  

### 3️⃣ Normalize VCF  

- Use `bcftools norm` to standardize format  

### 4️⃣ Variant Recalibration  

- Run `VariantRecalibrator` with known resources (e.g. HapMap, dbSNP)  

### 5️⃣ Apply VQSR  

- Apply recalibration using `ApplyVQSR`  

### 6️⃣ Filter Variants  

- Remove low-quality variants (`SelectVariants`)  

### 7️⃣ Split SNPs & INDELs  

- Use GATK to create separate SNP and INDEL files  

### 8️⃣ Apply Variant Filtration  

- Add filters using `VariantFiltration`  

### 9️⃣ Remove Filtered Variants  

- Exclude flagged entries from each VCF  

### 🔟 Merge Final SNP & INDELs  

- Recombine into a final cleaned VCF  

### 1️⃣1️⃣ Clean Chromosome Names  

- Strip `chr` prefix using `sed`  

### 1️⃣2️⃣ Normalize Again  

- Final normalization sweep with `bcftools`  

### 1️⃣3️⃣ Annotate with rs IDs  

- Use `bcftools annotate` to add rs identifiers  

### 1️⃣4️⃣ Predict Variant Effects (VEP)  

- Run `VEP` for functional annotation  

### 1️⃣5️⃣ Add ACMG & ClinGene Annotations (GeneBe)  

- Use `GeneBeClient` for pathogenicity predictions  

### 1️⃣6️⃣ Filter Final Set  

- Use `filter_variants_for_pathogenicity.sh` to retain only relevant variants:

    ✅ Pathogenic

    ❌ Non-pathogenic

    🧬 Protein-truncating variants (PTVs)
   

### 1️⃣7️⃣ Convert VCF to TSV  

- Use `transform_VCF.sh` to make data R-friendly  

### 1️⃣8️⃣ Statistical Tests  

- Apply Wilcoxon, Fisher, and t-tests  

### 1️⃣9️⃣ Enrichment Analysis  

- Use Enrichr, Enrichment Map, etc. for biological interpretation  

### 2️⃣0️⃣ Network Analysis  

- Visualize relationships with GeneMANIA, STRING, HumanBase  

---

## 💻 Command Summary  

### 🔹 Run Full VCF Pipeline (Steps 1–13)  
```bash
sbatch SNV_all.sh
```

### 🔹 Run VEP Annotation  
```bash
sbatch VEP_long.sh
sbatch VEP_short.sh
```

### 🔹 Run GeneBe Annotation 
```bash
java -jar GeneBeClient-0.1.0-a.4.jar vcf annotate \
--input-vcf ___.vcf.gz \
--output-vcf try.vcf \
--genome hg38 \
--api-key [YOUR_API_KEY] \
--username [YOUR_USERNAME]
```
### 🔹 Filter Variants for Pathogenicity
```bash
./filter_variants_for_pathogenicity.sh -i input.vcf.gz -o filtered_output.vcf.gz
```

### 🔹 Convert VCF to TSV  
```bash
transform_VCF.sh -i path/to/your_file.vcf.gz
```

Examples:
```bash
transform_VCF.sh -i pred_PTV_VEP_filtered.vcf.gz
transform_VCF.sh -i pred_patho_MPC_ALPHAMISSENSE_VEP_filtered.vcf.gz
```

### 🔹 Run Statistical Analyses  
```r
Rscript fisher_stat_scz_detail.r
Rscript statistics_fisher_new_correct.r
Rscript Statistics_fisher_scz_3_groups.r
Rscript statistics_wilcoxon_SCZ_3_groups.r
Rscript statistics_wilxocoxon_right.r
Rscript statstics_wilcoxon_SCZ.r
```

### 🔹 Run Enrichment & Network Analysis  
```r
Rscript GENE_lists.r
```

---

## 🎉 Final Thoughts  

This pipeline gives you an end-to-end workflow for SNV analysis — from raw variants to deep biological insights. It’s optimized for clarity, reproducibility, and ease of use. We hope it saves you time and brings those variants to life! 💡🔬

📫 Feel free to fork, contribute, or reach out with questions or ideas! 😊  

---


