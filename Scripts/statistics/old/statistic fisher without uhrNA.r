#statistic fisher without uhrNA

# ALLA STATISTICS 

# Load necessary libraries for data manipulation and analysis
library(readr)
library(data.table)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)

#! Load my data from CSV files
SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
GWAS_120 <- read.delim("~/GWAS_120.csv")       
BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")
SFARI <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SFARI.csv")

# !Clean GENESETS
# Convert gene names to a standard format for easier comparison
# BipEx_Bipolar
convert_col <- gconvert(query = BipEx_Bipolar$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
BipEx_Bipolar <- merge(BipEx_Bipolar, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
BipEx_Bipolar <- BipEx_Bipolar %>% select(Gene, name, everything())

# SCHEMA
convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
SCHEMA <- SCHEMA %>% select(Gene, name, everything())

# General clean-up of data
# Remove unnecessary columns, rename columns, and clean up gene names
BipEx_Bipolar <- BipEx_Bipolar %>% subset(select = -Gene) %>% rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.))) 
BipEx_Bipolar_p_val_PTV <- BipEx_Bipolar[!is.na(BipEx_Bipolar$PTV.Fisher.p.val) & BipEx_Bipolar$PTV.Fisher.p.val <= 0.05,]
BipEx_Bipolar_p_val_Missense <- BipEx_Bipolar[!is.na(BipEx_Bipolar$Damaging.Missense.Fisher.p.val) & BipEx_Bipolar$Damaging.Missense.Fisher.p.val <= 0.05,]

SCHEMA <- SCHEMA %>% subset(select = -Gene) %>% rename('Gene' = name) %>% mutate(Gene = gsub(' ', '', Gene)) %>% setNames(make.names(colnames(.)))  
SCHEMA_OR <- SCHEMA[SCHEMA$OR..Class.I. >= 1, ]
SCHEMA_pVal <- SCHEMA[SCHEMA$P.meta <= 0.05, ]
SFARI <- SFARI %>% subset(select = -c(status, `ensembl-id`)) %>% rename("Gene" = `gene-symbol`, "gene_name" = `gene-name`) %>% mutate(Gene = gsub(' ', '', Gene))
SFARI_syndromic <- SFARI 
SFARI_non_syndromic_lower_3 <- SFARI[SFARI$`gene-score` < 3, ]
GWAS_120 <- GWAS_120 %>% rename('Gene' = GENE_name)

# Combine all the dataframes into a single list of unique genes
genes_schema_pval <- unique(SCHEMA_pVal$Gene)
genes_bipolar <- unique(c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene))
genes_sfari_non_syndromic <- unique(SFARI_non_syndromic_lower_3$Gene)
genes_gwas <- unique(GWAS_120$Gene)
genes_schema_or <- unique(SCHEMA_OR$Gene)
genes_sfari_syndromic <- unique(SFARI_syndromic$Gene)


#!LOAD PERSONAL DATA 
# Specify the folder containing your TSV files
folder_path <- "/home/rachele/vcf_csv"

# List all the TSV files in the folder
tsv_files <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)

# Loop through each file and assign it to an object
for (file in tsv_files) {
  file_name <- tools::file_path_sans_ext(basename(file))
  assign(file_name, read_tsv(file))
}

# Get the names of all dataframes in the environment
dataframes <- ls(pattern = ".*")

#! General clean up  
# Loop through each dataframe and modify the SAMPLES column
for (df_name in dataframes) {
  df <- get(df_name)
  if ("SAMPLES" %in% names(df)) {
    df$SAMPLES <- gsub("_pool", "", df$SAMPLES)
    assign(df_name, df)
  }
}

#!ADD LABELS 
# Function to get the labels (case or control) for the outlier values
get_outlier_labels <- function(outlier_value, manifest_df) {
  outlier_values <- strsplit(outlier_value, ",")[[1]]
  labels <- sapply(outlier_values, function(val) {
    if (val %in% manifest_df$Sequencing_number) {
      return(manifest_df$new_column[manifest_df$Sequencing_number == val])
    } else {
      return(NA)
    }
  })
  return(paste(labels, collapse = ", "))
}

# Apply the get_outlier_labels function to the dataframes
pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_2 <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
pred_PTV_pLI_VEP_filtered$sample_label_2 <- sapply(pred_PTV_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
pred_PTV_VEP_filtered$sample_label_2 <- sapply(pred_PTV_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)

process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or) {
  
  # Separate rows and clean up sample IDs
  expanded_df <- df %>%
    separate_rows(SAMPLES, sample_label_2, sep = ",") %>%
    mutate(
      sample_id = trimws(SAMPLES),
      group = trimws(sample_label_2)
    ) %>%
    distinct()
  
  # Filter out UHR_NA group
  expanded_df <- expanded_df %>% filter(group != "UHR_NA")
  
  # Filter based on the genes in the brain, brain filter_ntpm, and additional gene sets
  expanded_df <- expanded_df %>%
    mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
           brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
           schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
           bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
           sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
           schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))
  
  # Calculate the number of variants, individuals, and genes for each group
  num_non_pathogenic_pLi_variants <- expanded_df %>%
    group_by(group) %>%
    summarize(num_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_non_pathogenic_pLi <- expanded_df %>%
    group_by(group) %>%
    summarize(unique_individuals = n_distinct(sample_id), .groups = "drop")
  
  num_brain_variants <- expanded_df %>%
    filter(brain == 1) %>%
    group_by(group) %>%
    summarize(num_brain_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_brain_variants <- expanded_df %>%
    filter(brain == 1) %>%
    group_by(group) %>%
    summarize(unique_individuals_brain = n_distinct(sample_id), .groups = "drop")
  
  num_brain_filter_ntpm_variants <- expanded_df %>%
    filter(brain_filter_ntpm == 1) %>%
    group_by(group) %>%
    summarize(num_brain_filter_ntpm_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_brain_filter_ntpm_variants <- expanded_df %>%
    filter(brain_filter_ntpm == 1) %>%
    group_by(group) %>%
    summarize(unique_individuals_brain_filter_ntpm = n_distinct(sample_id), .groups = "drop")
  
  num_schema_pval_variants <- expanded_df %>%
    filter(schema_pval == 1) %>%
    group_by(group) %>%
    summarize(num_schema_pval_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_schema_pval_variants <- expanded_df %>%
    filter(schema_pval == 1) %>%
    group_by(group) %>%
    summarize(unique_individuals_schema_pval = n_distinct(sample_id), .groups = "drop")
  
  num_bipolar_variants <- expanded_df %>%
    filter(bipolar == 1) %>%
    group_by(group) %>%
    summarize(num_bipolar_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_bipolar_variants <- expanded_df %>%
    filter(bipolar == 1) %>%
    group_by(group) %>%
    summarize(unique_individuals_bipolar = n_distinct(sample_id), .groups = "drop")
  
  num_sfari_non_syndromic_variants <- expanded_df %>%
    filter(sfari_non_syndromic == 1) %>%
    group_by(group) %>%
    summarize(num_sfari_non_syndromic_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_sfari_non_syndromic_variants <- expanded_df %>%
    filter(sfari_non_syndromic == 1) %>%
    group_by(group) %>%
    summarize(unique_individuals_sfari_non_syndromic = n_distinct(sample_id), .groups = "drop")
  
  num_schema_or_variants <- expanded_df %>%
    filter(schema_or == 1) %>%
    group_by(group) %>%
    summarize(num_schema_or_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_schema_or_variants <- expanded_df %>%
    filter(schema_or == 1) %>%
    group_by(group) %>%
    summarize(unique_individuals_schema_or = n_distinct(sample_id), .groups = "drop")
  
  # Ensure unique genes are counted
  num_genes <- expanded_df %>%
    group_by(group) %>%
    summarize(num_genes = n_distinct(SYMBOL), .groups = "drop")
  
  num_brain_genes <- expanded_df %>%
    filter(brain == 1) %>%
    group_by(group) %>%
    summarize(num_brain_genes = n_distinct(SYMBOL), .groups = "drop")
  
  num_brain_filter_ntpm_genes <- expanded_df %>%
    filter(brain_filter_ntpm == 1) %>%
    group_by(group) %>%
    summarize(num_brain_filter_ntpm_genes = n_distinct(SYMBOL), .groups = "drop")
  
  num_schema_pval_genes <- expanded_df %>%
    filter(schema_pval == 1) %>%
    group_by(group) %>%
    summarize(num_schema_pval_genes = n_distinct(SYMBOL), .groups = "drop")
  
  num_bipolar_genes <- expanded_df %>%
    filter(bipolar == 1) %>%
    group_by(group) %>%
    summarize(num_bipolar_genes = n_distinct(SYMBOL), .groups = "drop")
  
  num_sfari_non_syndromic_genes <- expanded_df %>%
    filter(sfari_non_syndromic == 1) %>%
    group_by(group) %>%
    summarize(num_sfari_non_syndromic_genes = n_distinct(SYMBOL), .groups = "drop")
  
  num_schema_or_genes <- expanded_df %>%
    filter(schema_or == 1) %>%
    group_by(group) %>%
    summarize(num_schema_or_genes = n_distinct(SYMBOL), .groups = "drop")
  
  # Define the row names based on the dataframe's name
  row_names <- c(
    paste("Number of", name ,"variants"),
    paste("Number of individuals with",name ,"variants"),
    paste("Number of brain", name ,"variants"),
    paste("Number of individuals with brain",name ,"variants"),
    paste("Number of brain filter_ntpm", name ,"variants"),
    paste("Number of individuals with brain filter_ntpm",name ,"variants"),
    paste("Number of schema_pval", name ,"variants"),
    paste("Number of individuals with schema_pval",name ,"variants"),
    paste("Number of bipolar", name ,"variants"),
    paste("Number of individuals with bipolar",name ,"variants"),
    paste("Number of sfari_non_syndromic", name ,"variants"),
    paste("Number of individuals with sfari_non_syndromic",name ,"variants"),
    paste("Number of schema_or", name ,"variants"),
    paste("Number of individuals with schema_or",name ,"variants"),
    paste("Number of genes in", name),
    paste("Number of brain genes in", name),
    paste("Number of brain filter_ntpm genes in", name),
    paste("Number of schema_pval genes in", name),
    paste("Number of bipolar genes in", name),
    paste("Number of sfari_non_syndromic genes in", name),
    paste("Number of schema_or genes in", name)
  )
  
  # Create the result dataframe without UHR_NA
  result_df <- data.frame(
    row_names = row_names,
    
    Case = c(
      num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Case"],
      num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "Case"],
      num_brain_variants$num_brain_variants[num_brain_variants$group == "Case"],
      num_individuals_with_brain_variants$unique_individuals_brain[num_individuals_with_brain_variants$group == "Case"],
      num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "Case"],
      num_individuals_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "Case"],
      num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "Case"],
      num_individuals_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "Case"],
      num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "Case"],
      num_individuals_with_bipolar_variants$unique_individuals_bipolar[num_individuals_with_bipolar_variants$group == "Case"],
      num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "Case"],
      num_individuals_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "Case"],
      num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "Case"],
      num_individuals_with_schema_or_variants$unique_individuals_schema_or[num_individuals_with_schema_or_variants$group == "Case"],
      num_genes$num_genes[num_genes$group == "Case"],
      num_brain_genes$num_brain_genes[num_brain_genes$group == "Case"],
      num_brain_filter_ntpm_genes$num_brain_filter_ntpm_genes[num_brain_filter_ntpm_genes$group == "Case"],
      num_schema_pval_genes$num_schema_pval_genes[num_schema_pval_genes$group == "Case"],
      num_bipolar_genes$num_bipolar_genes[num_bipolar_genes$group == "Case"],
      num_sfari_non_syndromic_genes$num_sfari_non_syndromic_genes[num_sfari_non_syndromic_genes$group == "Case"],
      num_schema_or_genes$num_schema_or_genes[num_schema_or_genes$group == "Case"]
    ),
    
    Control = c(
      num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Control"],
      num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "Control"],
      num_brain_variants$num_brain_variants[num_brain_variants$group == "Control"],
      num_individuals_with_brain_variants$unique_individuals_brain[num_individuals_with_brain_variants$group == "Control"],
      num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "Control"],
      num_individuals_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "Control"],
      num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "Control"],
      num_individuals_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "Control"],
      num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "Control"],
      num_individuals_with_bipolar_variants$unique_individuals_bipolar[num_individuals_with_bipolar_variants$group == "Control"],
      num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "Control"],
      num_individuals_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "Control"],
      num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "Control"],
      num_individuals_with_schema_or_variants$unique_individuals_schema_or[num_individuals_with_schema_or_variants$group == "Control"],
      num_genes$num_genes[num_genes$group == "Control"],
      num_brain_genes$num_brain_genes[num_brain_genes$group == "Control"],
      num_brain_filter_ntpm_genes$num_brain_filter_ntpm_genes[num_brain_filter_ntpm_genes$group == "Control"],
      num_schema_pval_genes$num_schema_pval_genes[num_schema_pval_genes$group == "Control"],
      num_bipolar_genes$num_bipolar_genes[num_bipolar_genes$group == "Control"],
      num_sfari_non_syndromic_genes$num_sfari_non_syndromic_genes[num_sfari_non_syndromic_genes$group == "Control"],
      num_schema_or_genes$num_schema_or_genes[num_schema_or_genes$group == "Control"]
    )
  )
  
  return(result_df)
}

result_df_list <- list(
  process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
  process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
  process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
  process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
  process_df(pred_PTV_VEP_filtered, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
)

# Combine all result dataframes into one
final_result_df <- bind_rows(result_df_list)

# View the final result
print(final_result_df)

# Process each dataframe and create separate result tables
result_df_NON_patho_pLI <- process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_patho_pLI <- process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_patho <- process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_PTV_pLI <- process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_PTV <- process_df(pred_PTV_VEP_filtered, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)


