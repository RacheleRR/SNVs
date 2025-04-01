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


#!CREATE DATAFRAME TO BE USED FOR STATISTICS 
# Define a function to process each dataframe and return the result
# This function processes the data to separate rows, clean up sample IDs, and filter based on specific gene sets
process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or) {
  
  expanded_df <- df %>%
    separate_rows(SAMPLES, sample_label_2, sep = ",") %>%
    mutate(
        sample_id = trimws(SAMPLES),
        group = trimws(sample_label_2)
    ) %>%
    distinct()
  
  # Filter based on the genes in the brain, brain filter_ntpm, and additional gene sets
  expanded_df <- expanded_df %>%
    mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
           brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
           schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
           bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
           sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
           schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))
  
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
  
  # Ensure all groups have the same number of rows by filling missing values with zeros
  groups <- c("Case", "Control", "UHR_NA")
  num_non_pathogenic_pLi_variants <- num_non_pathogenic_pLi_variants %>% complete(group = groups, fill = list(num_variants = 0))
  num_individuals_with_non_pathogenic_pLi <- num_individuals_with_non_pathogenic_pLi %>% complete(group = groups, fill = list(unique_individuals = 0))
  num_brain_variants <- num_brain_variants %>% complete(group = groups, fill = list(num_brain_variants = 0))
  num_individuals_with_brain_variants <- num_individuals_with_brain_variants %>% complete(group = groups, fill = list(unique_individuals_brain = 0))
  num_brain_filter_ntpm_variants <- num_brain_filter_ntpm_variants %>% complete(group = groups, fill = list(num_brain_filter_ntpm_variants = 0))
  num_individuals_with_brain_filter_ntpm_variants <- num_individuals_with_brain_filter_ntpm_variants %>% complete(group = groups, fill = list(unique_individuals_brain_filter_ntpm = 0))
  num_schema_pval_variants <- num_schema_pval_variants %>% complete(group = groups, fill = list(num_schema_pval_variants = 0))
  num_individuals_with_schema_pval_variants <- num_individuals_with_schema_pval_variants %>% complete(group = groups, fill = list(unique_individuals_schema_pval = 0))
  num_bipolar_variants <- num_bipolar_variants %>% complete(group = groups, fill = list(num_bipolar_variants = 0))
  num_individuals_with_bipolar_variants <- num_individuals_with_bipolar_variants %>% complete(group = groups, fill = list(unique_individuals_bipolar = 0))
  num_sfari_non_syndromic_variants <- num_sfari_non_syndromic_variants %>% complete(group = groups, fill = list(num_sfari_non_syndromic_variants = 0))
  num_individuals_with_sfari_non_syndromic_variants <- num_individuals_with_sfari_non_syndromic_variants %>% complete(group = groups, fill = list(unique_individuals_sfari_non_syndromic = 0))
  num_schema_or_variants <- num_schema_or_variants %>% complete(group = groups, fill = list(num_schema_or_variants = 0))
  num_individuals_with_schema_or_variants <- num_individuals_with_schema_or_variants %>% complete(group = groups, fill = list(unique_individuals_schema_or = 0))
  num_genes <- num_genes %>% complete(group = groups, fill = list(num_genes = 0))
  num_brain_genes <- num_brain_genes %>% complete(group = groups, fill = list(num_brain_genes = 0))
  num_brain_filter_ntpm_genes <- num_brain_filter_ntpm_genes %>% complete(group = groups, fill = list(num_brain_filter_ntpm_genes = 0))
  num_schema_pval_genes <- num_schema_pval_genes %>% complete(group = groups, fill = list(num_schema_pval_genes = 0))
  num_bipolar_genes <- num_bipolar_genes %>% complete(group = groups, fill = list(num_bipolar_genes = 0))
  num_sfari_non_syndromic_genes <- num_sfari_non_syndromic_genes %>% complete(group = groups, fill = list(num_sfari_non_syndromic_genes = 0))
  num_schema_or_genes <- num_schema_or_genes %>% complete(group = groups, fill = list(num_schema_or_genes = 0))

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
  
  # Create the result dataframe
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
    ),
    
    UHR_NA = c(
        num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "UHR_NA"],
        num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "UHR_NA"],
        num_brain_variants$num_brain_variants[num_brain_variants$group == "UHR_NA"],
        num_individuals_with_brain_variants$unique_individuals_brain[num_individuals_with_brain_variants$group == "UHR_NA"],
        num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "UHR_NA"],
        num_individuals_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "UHR_NA"],
        num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "UHR_NA"],
        num_individuals_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "UHR_NA"],
        num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "UHR_NA"],
        num_individuals_with_bipolar_variants$unique_individuals_bipolar[num_individuals_with_bipolar_variants$group == "UHR_NA"],
        num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "UHR_NA"],
        num_individuals_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "UHR_NA"],
        num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "UHR_NA"],
        num_individuals_with_schema_or_variants$unique_individuals_schema_or[num_individuals_with_schema_or_variants$group == "UHR_NA"],
        num_genes$num_genes[num_genes$group == "UHR_NA"],
        num_brain_genes$num_brain_genes[num_brain_genes$group == "UHR_NA"],
        num_brain_filter_ntpm_genes$num_brain_filter_ntpm_genes[num_brain_filter_ntpm_genes$group == "UHR_NA"],
        num_schema_pval_genes$num_schema_pval_genes[num_schema_pval_genes$group == "UHR_NA"],
        num_bipolar_genes$num_bipolar_genes[num_bipolar_genes$group == "UHR_NA"],
        num_sfari_non_syndromic_genes$num_sfari_non_syndromic_genes[num_sfari_non_syndromic_genes$group == "UHR_NA"],
        num_schema_or_genes$num_schema_or_genes[num_schema_or_genes$group == "UHR_NA"]
    )
  )
  
  return(result_df)
}

# Process each dataframe
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


#!statitics


#!  VARIANTS CONFRONATION
# Fisher's exact test is used to determine if there are nonrandom associations between two categorical variables.
# Chi-squared test is used to determine if there is a significant association between two categorical variables.
# It tests for associations between the presence of specific variants and the case/control status.

# Hypotheses:
# - Null Hypothesis (H0): The distribution of brain PTV variants between cases and controls is the same, 
#   meaning there is no significant association.
# - Alternative Hypothesis (H1): The distribution of brain PTV variants between cases and controls is different, 
#   indicating a potential association.

# Statistical Methods:
# - Fisher’s Exact Test: Determines if there is a nonrandom association between cases and controls with respect 
#   to variant presence.
# - Chi-Squared Test: Evaluates whether the observed distribution of variants differs significantly from the 
#   expected distribution.


# Contingency table:

#                                        case                         |    control
# ----------------------------|----------------------------------------------------------------
# variant in schizo                           x                       |      y 
# variant not in schizo          total variant in case - x            |      total variant in control -y


perform_variant_confrontation <- function(result_df, row) {
    variant_counts <- data.frame(
        Group = c("Cases", "Controls"),
        Brain_PTV_Variants = c(result_df[row, "Case"], result_df[row, "Control"]),
        Non_Brain_PTV_Variants = c(result_df[1, "Case"] - result_df[row, "Case"], result_df[1, "Control"] - result_df[row, "Control"])
    )
    
    
    fisher_result <- fisher.test(variant_counts[, 2:3])
    chi_squared_result <- chisq.test(variant_counts[, 2:3])
    
    return(list(fisher = fisher_result, chi_squared = chi_squared_result))
}

# Rename one of the perform_tests_for_all_rows functions to avoid conflict
perform_variant_confrontation_tests <- function(result_df) {
    rows <- c(3, 5, 7, 9, 11, 13)  # Rows to be tested
    row_names <- c("brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants", "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")
    
    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        res <- perform_variant_confrontation(result_df, row)
        list(row_name = row_name, result = res)
    })
    
    names(results) <- row_names
    return(results)
}

# Use the renamed function
result_var_conf <- list(
    NON_patho_pLI = perform_variant_confrontation_tests(result_df_NON_patho_pLI),
    patho_pLI = perform_variant_confrontation_tests(result_df_patho_pLI),
    patho = perform_variant_confrontation_tests(result_df_patho),
    PTV_pLI = perform_variant_confrontation_tests(result_df_PTV_pLI),
    PTV = perform_variant_confrontation_tests(result_df_PTV)
)

result_var_conf_df <- do.call(rbind, lapply(names(result_var_conf), function(test_name) {
    do.call(rbind, lapply(result_var_conf[[test_name]], function(res) {
        data.frame(
            Test = test_name,
            Row = res$row_name,
            fisher_p_value = res$result$fisher$p.value,
            fisher_odds_ratio = res$result$fisher$estimate,
            chi_squared_p_value = res$result$chi_squared$p.value,
            chi_squared_statistic = res$result$chi_squared$statistic
        )
    }))
}))


print(result_var_conf_df)

#!TO GET CONTINGENCY TABEL 
perform_variant_confrontation <- function(result_df, row) {
    variant_counts <- data.frame(
        Group = c("Cases", "Controls"),
        Brain_PTV_Variants = c(result_df[row, "Case"], result_df[row, "Control"]),
        Non_Brain_PTV_Variants = c(result_df[1, "Case"] - result_df[row, "Case"], result_df[1, "Control"] - result_df[row, "Control"])
    )
    
    fisher_result <- fisher.test(variant_counts[, 2:3])
    chi_squared_result <- chisq.test(variant_counts[, 2:3])
    
    return(list(fisher = fisher_result, chi_squared = chi_squared_result, variant_counts = variant_counts))
}

# Rename one of the perform_tests_for_all_rows functions to avoid conflict
perform_variant_confrontation_tests <- function(result_df) {
    rows <- c(3, 5, 7, 9, 11, 13)  # Rows to be tested
    row_names <- c("brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants", "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")
    
    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        res <- perform_variant_confrontation(result_df, row)
        list(row_name = row_name, result = res)
    })
    
    names(results) <- row_names
    return(results)
}

# Use the renamed function
result_var_conf <- list(
    NON_patho_pLI = perform_variant_confrontation_tests(result_df_NON_patho_pLI),
    patho_pLI = perform_variant_confrontation_tests(result_df_patho_pLI),
    patho = perform_variant_confrontation_tests(result_df_patho),
    PTV_pLI = perform_variant_confrontation_tests(result_df_PTV_pLI),
    PTV = perform_variant_confrontation_tests(result_df_PTV)
)

# Collect all contingency tables into a single data frame
contingency_tables_df <- do.call(rbind, lapply(names(result_var_conf), function(test_name) {
    do.call(rbind, lapply(result_var_conf[[test_name]], function(res) {
        variant_counts <- res$result$variant_counts
        variant_counts$Test <- test_name
        variant_counts$Row <- res$row_name
        return(variant_counts)
    }))
}))

# Print the contingency tables data frame
print(contingency_tables_df)


#! Perform Fisher's exact test for individuals carries at least 1 varaint
# This function calculates the Fisher's exact test and odds ratio for each variant type.
# This script tests whether cases are more likely than controls to carry at least one of the filtered variants.
# It does this using Fisher’s exact test, which determines whether there is a nonrandom association between 
# two categorical variables (case/control status and variant presence).

# Hypotheses:
# - Null Hypothesis (H0): Cases and controls have the same probability of carrying at least one filtered variant.
# - Alternative Hypothesis (H1): Cases are more or less likely than controls to carry at least one filtered variant.

# Statistical Method:
# - Fisher’s Exact Test: Used to evaluate whether the presence of a variant differs significantly between cases and controls.
# - Odds Ratio: Measures the strength of association between variant presence and case/control status.


#                Carries Variant | Does Not Carry Variant
# -------------------------------|------------------------
# Cases            a            |          b
# Controls         c            |          d

# Where:

#     a = Number of cases carrying at least one variant.
#     b = Number of cases not carrying any variants.
#     c = Number of controls carrying at least one variant.
#     d = Number of controls not carrying any variants.



# Function to perform Fisher's exact test for a single variant type.
perform_fisher_test_per_ind <- function(variant_name, case_count, control_count, total_cases, total_controls) {
  # Define counts for contingency table
  a <- case_count  # Number of cases carrying the variant
  c <- control_count  # Number of controls carrying the variant
  b <- total_cases - a  # Number of cases NOT carrying the variant
  d <- total_controls - c  # Number of controls NOT carrying the variant

  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  colnames(contingency_table) <- c("Carries Variant", "Does Not Carry Variant")
  rownames(contingency_table) <- c("Cases", "Controls")

  fisher_test <- fisher.test(contingency_table)
  odds_ratio <- (a / b) / (c / d)

  return(list(p_value = fisher_test$p.value, odds_ratio = odds_ratio))
}

# Total number of cases and controls
total_cases <- 302
total_controls <- 75

# This function applies the perform_fisher_test function to multiple rows in the result dataframe.
# It tests for associations between the presence of specific variants and the case/control status.

perform_fisher_tests_for_all_rows <- function(result_df) {
  rows <- c(2, 4, 6, 8, 10, 12, 14)  # Rows to be tested
  results <- sapply(rows, function(i) {
    result <- perform_fisher_test_per_ind(
      variant_name = result_df$row_names[i],
      case_count = result_df$Case[i],
      control_count = result_df$Control[i],
      total_cases = total_cases,
      total_controls = total_controls
    )
    return(c(result$p_value, result$odds_ratio))
  })
  return(results)
}

# Apply the function to each result dataframe and compile the results
fisher_results_ind <- list(
  NON_patho_pLI = perform_fisher_tests_for_all_rows(result_df_NON_patho_pLI),
  patho_pLI = perform_fisher_tests_for_all_rows(result_df_patho_pLI),
  patho = perform_fisher_tests_for_all_rows(result_df_patho),
  PTV_pLI = perform_fisher_tests_for_all_rows(result_df_PTV_pLI),
  PTV = perform_fisher_tests_for_all_rows(result_df_PTV)
)

# Create a new dataframe to store the Fisher's test results
fisher_results_ind_df <- do.call(rbind, lapply(names(fisher_results_ind), function(test_name) {
  data.frame(
    Test = test_name,
    Row = c("filtered", "brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants", "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants"),
    fisher_p_value = fisher_results_ind[[test_name]][1, ],
    fisher_odds_ratio = fisher_results_ind[[test_name]][2, ]
  )
}))

# Print the Fisher's test results dataframe
print(fisher_results_ind_df)


#! FISHER 
# Option 1: Case-Control vs. Brain/Not in Brain
# This table tests whether variants in the brain are enriched in cases versus controls.
# If you suspect that variants present in brain-expressed genes (or brain tissues) contribute more to disease, this table makes sense.
# You're asking: Are brain variants more common in cases than controls?

#                 Variant in brain | Variant not in brain 
# -------------------------------|------------------------
# Cases            a            |          b
# Controls         c            |          d



perform_variant_confrontation_enrich <- function(result_df, row) {
    variant_counts <- data.frame(
        Variant = c(result_df[3, "Case"], result_df[row, "Control"]),
        Non_Variant = c(result_df[1, "Case"] - result_df[row, "Case"],
                        result_df[1, "Control"] - result_df[row, "Control"])
    )
    rownames(variant_counts) <- c("Cases", "Controls")

    fisher_result <- fisher.test(variant_counts)
    chi_squared_result <- chisq.test(variant_counts)
    
    return(list(fisher = fisher_result, chi_squared = chi_squared_result))
}

# Rename one of the perform_tests_for_all_rows functions to avoid conflict
perform_variant_confrontation_enrich_tests <- function(result_df) {
    rows <- c(3, 5, 7, 9, 11, 13)  # Rows to be tested
    row_names <- c("brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants", "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")
    
    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        res <- perform_variant_confrontation_enrich(result_df, row)
        list(row_name = row_name, result = res)
    })
    
    names(results) <- row_names
    return(results)
}

# Use the renamed function
result_var_conf_enrich <- list(
    NON_patho_pLI = perform_variant_confrontation_enrich_tests(result_df_NON_patho_pLI),
    patho_pLI = perform_variant_confrontation_enrich_tests(result_df_patho_pLI),
    patho = perform_variant_confrontation_enrich_tests(result_df_patho),
    PTV_pLI = perform_variant_confrontation_enrich_tests(result_df_PTV_pLI),
    PTV = perform_variant_confrontation_enrich_tests(result_df_PTV)
)

result_var_conf_enrich_df <- do.call(rbind, lapply(names(result_var_conf_enrich), function(test_name) {
    do.call(rbind, lapply(result_var_conf_enrich[[test_name]], function(res) {
        data.frame(
            Test = test_name,
            Row = res$row_name,
            fisher_p_value = res$result$fisher$p.value,
            fisher_odds_ratio = res$result$fisher$estimate,
            chi_squared_p_value = res$result$chi_squared$p.value,
            chi_squared_statistic = res$result$chi_squared$statistic
        )
    }))
}))



#!FISHER varaint confrontation normalized 
#NOT GOOD BECAUSE ITS ROUNDING IT UP 
# Contingency table:

#                             |          case                                        |                                 control
# ----------------------------|----------------------------------------------------------------
# variant in schizo                    n var Schizo /totalcase                       |                 n var Schizo/total control 
# variant not in schizo          total variant/totcase - n var Schizo /total         |    total variant/totcontrol -n var Schizo/total control 




perform_variant_confr_normal   <- function(result_df, row) {
    total_cases <- 302
    total_controls <- 75

    # Calculate rate of variants per individual
    case_percentage <- (result_df[row, "Case"] / total_cases) 
    control_percentage <- (result_df[row, "Control"] / total_controls)

    non_case_percentage <-  (result_df[1 , "Case"] / total_cases) - case_percentage
    non_control_percentage <-  (result_df[1  , "Control"] / total_controls) - control_percentage

    # Create contingency table using percentages
    variant_percentages <- data.frame(
        Group = c("Cases", "Controls"),
        Brain_PTV_Percentage = c(case_percentage, control_percentage),
        Non_Brain_PTV_Percentage = c(non_case_percentage, non_control_percentage)
    )

    # Perform Fisher's exact test
    fisher_result <- fisher.test(variant_percentages[, 2:3])

    return(list(fisher = fisher_result))
}

# Rename the other perform_tests_for_all_rows function to avoid conflict
perform_variant_confrontation_normalized_tests <- function(result_df) {
    rows <- c(3, 5, 7, 9, 11, 13)  # Rows to be tested
    row_names <- c("brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants", "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")

    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        res <- perform_variant_confr_normal(result_df, row)
        list(row_name = row_name, result = res)
    })
    
    names(results) <- row_names
    return(results)
}

# Use the renamed function
result_var_conf_norm <- list(
    NON_patho_pLI = perform_variant_confrontation_normalized_tests(result_df_NON_patho_pLI),
    patho_pLI = perform_variant_confrontation_normalized_tests(result_df_patho_pLI),
    patho = perform_variant_confrontation_normalized_tests(result_df_patho),
    PTV_pLI = perform_variant_confrontation_normalized_tests(result_df_PTV_pLI),
    PTV = perform_variant_confrontation_normalized_tests(result_df_PTV)
)

result_var_conf_norm_df <- do.call(rbind, lapply(names(result_var_conf_norm), function(test_name) {
    do.call(rbind, lapply(result_var_conf_norm[[test_name]], function(res) {
        data.frame(
            Test = test_name,
            Row = res$row_name,
            fisher_p_value = res$result$fisher$p.value,
            fisher_odds_ratio = res$result$fisher$estimate
        )
    }))
}))



#! rate comparison 
# Function to calculate case and control rates
rate_confront <- function(result_df, row) {
    total_cases <- 302
    total_controls <- 75

    # Calculate rate of variants per individual
    case_rate <- result_df[row, "Case"] / total_cases
    control_rate <- result_df[row, "Control"] / total_controls

    return(list(case_rate = case_rate, control_rate = control_rate))
}

# Function to apply the test to specified rows
perform_tests_for_all_rows <- function(result_df) {
    rows <- c(1, 3, 5, 7, 9, 11, 13)  # Rows to be tested
    row_names <- c("total", "brain_variants", "brain_filter_ntpm_variants", 
                   "schema_pval_variants", "bipolar_variants", 
                   "sfari_non_syndromic_variants", "schema_or_variants")

    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        rates <- rate_confront(result_df, row)
        data.frame(Row = row_name, Case_Rate = rates$case_rate, Control_Rate = rates$control_rate)
    })

    return(do.call(rbind, results))  # Combine results into a data frame
}

# Run the tests on multiple datasets
rate_confront <- list(
    NON_patho_pLI = perform_tests_for_all_rows(result_df_NON_patho_pLI),
    patho_pLI = perform_tests_for_all_rows(result_df_patho_pLI),
    patho = perform_tests_for_all_rows(result_df_patho),
    PTV_pLI = perform_tests_for_all_rows(result_df_PTV_pLI),
    PTV = perform_tests_for_all_rows(result_df_PTV)
)

# Combine results into a single dataframe with a 'Test' column
rate_confront_df <- do.call(rbind, lapply(names(rate_confront), function(test_name) {
    df <- rate_confront[[test_name]]
    df$Test <- test_name  # Add test name as a column
    return(df)
}))

# Reorder columns for readability
rate_confront_df <- rate_confront_df[, c("Test", "Row", "Case_Rate", "Control_Rate")]

# Print the final table
print(rate_confront_df)




#!POISSON REGRESSION
#!DOESNT HAVE SENSE GIVEN THAT RATE IS BAD 

# Create a dataframe
data <- data.frame(
  group = c(rep("Case", 2), rep("Control", 2)),
  schizophrenia_variant = c(X, total_case_variants - X, Y, total_control_variants - Y),
  individuals = c(302, 302, 75, 75)
)

# Poisson regression
model <- glm(schizophrenia_variant ~ group + offset(log(individuals)), family = poisson, data = data)
summary(model)









# ! Perform binomial test for specific rows
# This script performs binomial tests to determine whether the proportion of individuals with a specific type of variants in cases 
# and controls significantly deviates from an expected distribution.
# The binomial test is used to determine if the proportion of successes in a sample is significantly different from a specified proportion.
# The binomial test is used when analyzing the probability of a certain number of successes in a fixed number of trials.
# A binomial test is particularly useful when you are interested in the proportion of successes (in this case, brain PTV variants) out of a given total (the number of individuals).

# Hypotheses:
# - Null Hypothesis (H0): The proportion of individuals with a specific type of PTV variants in cases and controls is equal to the expected proportion 
#   (based on the overall variant distribution or another baseline).
# - Alternative Hypothesis (H1): The proportion of individuals with a specific type of PTV variants in cases and controls is different from the expected proportion.

# Statistical Method:
# - Binomial Test: Evaluates whether the observed proportion of individuals carrying brain PTV variants differs significantly 
#   from an expected probability.


# Function to perform a binomial test for a specific row in the dataset.
perform_binomial_test <- function(result_df, row) {
    # Extract the number of successes (individuals with the variant)
    success_case <- result_df[row, "Case"]
    success_control <- result_df[row, "Control"]

    # Extract the total number of trials (total individuals considered)
    trials_case <- result_df[2, "Case"]  
    trials_control <- result_df[2, "Control"]

    # Compute expected proportions based on the previous row
    expected_p_case <- result_df[row - 1, "Case"] / result_df[1, "Case"]
    expected_p_control <- result_df[row - 1, "Control"] / result_df[1, "Control"]

    # Debugging output to check values
    print(paste("Row:", row))
    print(paste("Success Case:", success_case))
    print(paste("Trials Case:", trials_case))
    print(paste("Expected P Case:", expected_p_case))
    print(paste("Success Control:", success_control))
    print(paste("Trials Control:", trials_control))
    print(paste("Expected P Control:", expected_p_control))

    # Perform binomial tests if the number of trials is valid
    if (trials_case > 0 && trials_case >= success_case && trials_control > 0 && trials_control >= success_control) {
        binom_case <- binom.test(success_case, trials_case, p = expected_p_case)
        binom_control <- binom.test(success_control, trials_control, p = expected_p_control)
        
        return(list(case = binom_case, control = binom_control))
    } else {
        return(list(case = NA, control = NA))  # Return NA if invalid trials
    }
}

# Apply binomial tests for all specified rows
# It tests if the proportion of individuals with specific variants is significantly different from the expected proportion.

perform_binomial_tests_for_all_rows <- function(result_df) {
    rows <- c(4, 6, 8, 10, 12, 14)  # Rows to be tested
    row_names <- c("brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants", "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")
    
    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        res <- perform_binomial_test(result_df, row)
        list(row_name = row_name, result = res)
    })
    
    names(results) <- row_names
    return(results)
}

# Apply the function to each result dataframe and compile the results

binomial_results <- list(
    NON_patho_pLI = perform_binomial_tests_for_all_rows(result_df_NON_patho_pLI),
    patho_pLI = perform_binomial_tests_for_all_rows(result_df_patho_pLI),
    patho = perform_binomial_tests_for_all_rows(result_df_patho),
    PTV_pLI = perform_binomial_tests_for_all_rows(result_df_PTV_pLI),
    PTV = perform_binomial_tests_for_all_rows(result_df_PTV)
)

# Create a dataframe to store the binomial test results
binomial_results_df <- do.call(rbind, lapply(names(binomial_results), function(test_name) {
    do.call(rbind, lapply(binomial_results[[test_name]], function(res) {
        data.frame(
            Test = test_name,
            Row = res$row_name,
            case_p_value = res$result$case$p.value,
            case_estimate = res$result$case$estimate,
            control_p_value = res$result$control$p.value,
            control_estimate = res$result$control$estimate
        )
    }))
}))

# Print the binomial test results dataframe
print(binomial_results_df)

#! Perform Wilcoxon test only for general number of filtered variants (generated by VEP) 
#! SIMULATING THE proportion pro individual 
# This script performs Wilcoxon tests to compare the distribution of variants per individual between cases and controls.
# ergo if you want to compare the total number of variants carried by cases vs. controls.
# The Wilcoxon test is a non-parametric alternative to the t-test, used when the data do not follow a normal distribution.

# Hypotheses:
# - Null Hypothesis (H0): The distribution of variants per individual in cases and controls is the same.
# - Alternative Hypothesis (H1): The distribution of variants per individual in cases and controls is different.

# Statistical Method:
# - Wilcoxon Rank-Sum Test (Mann-Whitney U Test): Compares two independent groups (cases vs. controls) for differences in distribution.


#  Function to perform Wilcoxon test for the total number of variants per individual in a dataset.
perform_wilcoxon_test <- function(result_df) {
    case_variants <- result_df[1, "Case"]
    control_variants <- result_df[1, "Control"]
    case_individuals <- result_df[2, "Case"]
    control_individuals <- result_df[2, "Control"]
    
    case_ptv_per_individual <- case_variants / case_individuals
    control_ptv_per_individual <- control_variants / control_individuals
    
    # Perform Wilcoxon test
    wilcox_result <- wilcox.test(case_ptv_per_individual, control_ptv_per_individual, alternative = "two.sided")
    
    return(wilcox_result)
}

# Apply Wilcoxon tests for all result dataframes and compile the results
perform_wilcoxon_tests_for_all <- function(result_dfs) {
    wilcoxon_results <- lapply(names(result_dfs), function(test_name) {
        result_df <- result_dfs[[test_name]]
        wilcox_result <- perform_wilcoxon_test(result_df)
        data.frame(
            Test = test_name,
            p_value = wilcox_result$p.value,
            statistic = wilcox_result$statistic
        )
    })
    
    wilcoxon_results_df <- do.call(rbind, wilcoxon_results)
    return(wilcoxon_results_df)
}

# Apply the function to all result dataframes
wilcoxon_results_df <- perform_wilcoxon_tests_for_all(list(
    NON_patho_pLI = result_df_NON_patho_pLI,
    patho_pLI = result_df_patho_pLI,
    patho = result_df_patho,
    PTV_pLI = result_df_PTV_pLI,
    PTV = result_df_PTV
))

# Print the Wilcoxon test results dataframe
print(wilcoxon_results_df)

#! Perform Wilcoxon test for specific type of filtered variants (brain,SCHEMA Etc)
#! SIMULATING THE proportion pro individual 
# This function performs a Wilcoxon test for a given row in the choosen dataframe.
# This script performs Wilcoxon tests to compare the distribution of variants per individual between cases and controls.
# ergo if you want to compare the total number of variants carried by cases vs. controls.
# The Wilcoxon test is a non-parametric alternative to the t-test, used when the data do not follow a normal distribution.

# Hypotheses:
# - Null Hypothesis (H0): The distribution of variants per individual in cases and controls is the same.
# - Alternative Hypothesis (H1): The distribution of variants per individual in cases and controls is different.

# Statistical Method:
# - Wilcoxon Rank-Sum Test (Mann-Whitney U Test): Compares two independent groups (cases vs. controls) for differences in distribution.


perform_wilcoxon_test <- function(result_df, row) {
    case_variants <- result_df[row, "Case"]
    control_variants <- result_df[row, "Control"]
    case_individuals <- result_df[row + 1, "Case"]
    control_individuals <- result_df[row + 1, "Control"]
    
    case_ptv_per_individual <- case_variants / case_individuals
    control_ptv_per_individual <- control_variants / control_individuals
    
    # Ensure there are enough observations for the test
    if (!is.na(case_ptv_per_individual) && !is.na(control_ptv_per_individual) && case_individuals > 1 && control_individuals > 1) {
        # Perform Wilcoxon test
        wilcox_result <- wilcox.test(case_ptv_per_individual, control_ptv_per_individual, alternative = "two.sided")
    } else {
        wilcox_result <- list(p.value = NA, statistic = NA)
    }
    
    return(list(wilcox_result = wilcox_result, case_ptv_per_individual = case_ptv_per_individual, control_ptv_per_individual = control_ptv_per_individual))
}

# Apply Wilcoxon tests for all specified rows
# It tests if the distribution of variants per individual is significantly different between cases and controls for specific variant types.
perform_wilcoxon_tests_for_all_rows <- function(result_df) {
    rows <- c(1, 3, 5, 7, 9, 11, 13)  # Rows to be tested
    row_names <- c("normal_variants", "brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants", "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")
    
    results <- lapply(seq_along(rows), function(i) {
        row <- rows[i]
        row_name <- row_names[i]
        wilcox_result <- perform_wilcoxon_test(result_df, row)
        list(row_name = row_name, result = wilcox_result$wilcox_result, case_ptv_per_individual = wilcox_result$case_ptv_per_individual, control_ptv_per_individual = wilcox_result$control_ptv_per_individual)
    })
    
    names(results) <- row_names
    return(results)
}

# Apply the function to each result dataframe and compile the results

wilcoxon_results <- list(
    NON_patho_pLI = perform_wilcoxon_tests_for_all_rows(result_df_NON_patho_pLI),
    patho_pLI = perform_wilcoxon_tests_for_all_rows(result_df_patho_pLI),
    patho = perform_wilcoxon_tests_for_all_rows(result_df_patho),
    PTV_pLI = perform_wilcoxon_tests_for_all_rows(result_df_PTV_pLI),
    PTV = perform_wilcoxon_tests_for_all_rows(result_df_PTV)
)

# Create a dataframe to store the Wilcoxon test results
wilcoxon_results_df <- do.call(rbind, lapply(names(wilcoxon_results), function(test_name) {
    do.call(rbind, lapply(wilcoxon_results[[test_name]], function(res) {
        data.frame(
            Test = test_name,
            Row = res$row_name,
            p_value = res$result$p.value,
            statistic = res$result$statistic,
            case_ptv_per_individual = res$case_ptv_per_individual,
            control_ptv_per_individual = res$control_ptv_per_individual
        )
    }))
}))

# Print the Wilcoxon test results dataframe
print(wilcoxon_results_df)


