# to test whether cases tend to carry more variants than controls, you need to analyze the counts of variants per individual instead of just binary presence/absence.
 

library(readr)
library(data.table)
# library(tidyverse)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)

#! Load my data   
SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
GWAS_120 <- read.delim("~/GWAS_120.csv")       
BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")
SFARI <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SFARI.csv")

# !Clean GENESETS
# Convert gene names 
# BipEx_Bipolar
convert_col <- gconvert(query = BipEx_Bipolar$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
BipEx_Bipolar <- merge(BipEx_Bipolar, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
BipEx_Bipolar <- BipEx_Bipolar %>% select(Gene, name, everything())

# SCHEMA
convert_col <- gconvert(query = SCHEMA$Gene, organism = "hsapiens", target = "ENSG", mthreshold = Inf, filter_na = TRUE)
SCHEMA <- merge(SCHEMA, convert_col[, c("target", "name")], by.x = "Gene", by.y = "target", all.x = TRUE)
SCHEMA <- SCHEMA %>% select(Gene, name, everything())

# General clean 
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

# Unify all the dataframes 
genes_schema_pval <- unique(SCHEMA_pVal$Gene)
genes_bipolar <- unique(c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene))
genes_sfari_non_syndromic <- unique(SFARI_non_syndromic_lower_3$Gene)
genes_gwas <- unique(GWAS_120$Gene)
genes_schema_or <- unique(SCHEMA_OR$Gene)
genes_sfari_syndromic <- unique(SFARI_syndromic$Gene)


#!LOAD MY DATA

# Specify the folder containing your TSV files
folder_path <- "/home/rachele/vcf_csv"

# List all the TSV files in the folder
tsv_files <- list.files(path = folder_path, pattern = "\\.tsv$", full.names = TRUE)

# Loop through each file and assign it to an object
for (file in tsv_files) {
  # Create a dynamic object name based on the file name (without the extension)
  file_name <- tools::file_path_sans_ext(basename(file))
  
  # Read the TSV file and assign it to a new object
  assign(file_name, read_tsv(file))
}

# Get the names of all dataframes in the environment
dataframes <- ls(pattern = ".*")  # Adjust the pattern if needed to only target your loaded dataframes

#!GENERAL CLEAN
# Loop through each dataframe and modify the SAMPLES column
for (df_name in dataframes) {
  # Use get() to access the dataframe
  df <- get(df_name)
  
  # Check if the SAMPLES column exists
  if ("SAMPLES" %in% names(df)) {
    # Apply gsub to modify the SAMPLES column
    df$SAMPLES <- gsub("_pool", "", df$SAMPLES)
    
    # Assign the modified dataframe back to the same object name
    assign(df_name, df)
  }
}


#! get labels 

# Function to get the labels (case or control) for the outlier values
get_outlier_labels <- function(outlier_value, manifest_df) {
    # Split the outlier values by ';' into a vector
    outlier_values <- strsplit(outlier_value, ",")[[1]]
    
    # Retrieve corresponding labels from the manifest dataframe
    labels <- sapply(outlier_values, function(val) {
        # Find the label in the manifest dataframe (either "case" or "control")
        if (val %in% manifest_df$Sequencing_number) {
            return(manifest_df$new_column[manifest_df$Sequencing_number == val])
        } else {
            return(NA) # Return NA if no match is found
        }
    })
    
    # Concatenate the labels into a single string, separated by commas
    return(paste(labels, collapse = ", "))
}


pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_2 <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$SAMPLES,get_outlier_labels,manifest_correct)
pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$SAMPLES,get_outlier_labels,manifest_correct)
pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$SAMPLES,get_outlier_labels,manifest_correct)
pred_PTV_pLI_VEP_filtered$sample_label_2 <- sapply(pred_PTV_pLI_VEP_filtered$SAMPLES,get_outlier_labels,manifest_correct)
pred_PTV_VEP_filtered$sample_label_2 <- sapply(pred_PTV_VEP_filtered$SAMPLES,get_outlier_labels,manifest_correct)


#! CREATE DATAFRAME FOR STATISTICAL TESTS TO COME 

process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or) {
    
    # Assume expanded_df is already created as per your code
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
    
    # Count the number of variants per individual for each category
    result_df <- expanded_df %>%
        group_by(sample_id, group) %>%
        summarise(
            Number_of_Variants = n(),
            Number_of_Brain_Variants = sum(brain),
            Number_of_Brain_Filter_NTPM_Variants = sum(brain_filter_ntpm),
            Number_of_Schema_Pval_Variants = sum(schema_pval),
            Number_of_Bipolar_Variants = sum(bipolar),
            Number_of_Sfari_Non_Syndromic_Variants = sum(sfari_non_syndromic),
            Number_of_Schema_Or_Variants = sum(schema_or),
            .groups = "drop"
        )
    
    return(result_df)
}

# Create separate data frames for each call to process_df
result_df_NON_patho_pLI <- process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_patho_pLI <- process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_patho <- process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_PTV_pLI <- process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_df_PTV <- process_df(pred_PTV_VEP_filtered, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)

colnames(manifest_correct) <- c("Status", "sample_id", "group")

#! ADD ALL SAMPLES INDEPENDENTLY OF PRESENCE/ABSENCE OF VARIANTS 
# Ensure consistency in column names and perform left join for all dataframes
complete_df_NON_patho_pLI <- manifest_correct %>%
    left_join(result_df_NON_patho_pLI, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_patho_pLI <- manifest_correct %>%
    left_join(result_df_patho_pLI, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_patho <- manifest_correct %>%
    left_join(result_df_patho, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_PTV_pLI <- manifest_correct %>%
    left_join(result_df_PTV_pLI, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_PTV <- manifest_correct %>%
    left_join(result_df_PTV, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))


#! Perform Wilcoxon test only for general number of filtered variants (generated by VEP) 
#! DOES NOT SIMULATE THE proportion pro individual 
# Problem Definition  
# We aim to compare the distribution of **variants per individual** between **cases and controls** 
# using two statistical tests:

# 1. **Welch’s t-test** – a parametric test that compares the means of two independent groups.  
# 2. **Wilcoxon Rank-Sum Test (Mann-Whitney U Test)** – a non-parametric alternative that compares distributions without assuming normality.  

# This approach helps determine whether **cases** carry significantly more variants than **controls**.

# Hypotheses  

## 1. Welch’s t-test (Parametric)
# Null Hypothesis (H₀): The mean number of variants per individual is the same for cases and controls.  
# Alternative Hypothesis (H₁): The mean number of variants per individual is different between cases and controls.  

## 2. Wilcoxon Rank-Sum Test (Non-Parametric)
# Null Hypothesis (H₀): The distribution of variants per individual is the same for cases and controls.  
# Alternative Hypothesis (H₁): The distribution of variants per individual is different between cases and controls.  

# Statistical Methods  

## 1. Welch’s t-test (Parametric Test)
# - Compares the **means** of two independent groups.
# - Does not assume equal variances (Welch’s correction is applied).
# - Sensitive to **outliers** and **non-normal distributions**.
# - Used when data approximately follows a normal distribution.

## 2. Wilcoxon Rank-Sum Test (Non-Parametric Test)
# - Compares the **ranks** of two independent groups instead of means.
# - Does not assume normality.
# - More **robust** to outliers and small sample sizes.
# - Detects differences in both **median** and **distribution shape**.


t_test_results_df_no_lev <- data.frame(
    DataFrame = character(),
    Column = character(),
    t_statistic = numeric(),
    df = numeric(),
    p_value_ttest = numeric(),
    mean_case = numeric(),
    mean_control = numeric(),
    stringsAsFactors = FALSE
)

wilcoxon_results_df_no_lev <- data.frame(
    DataFrame = character(),
    Column = character(),
    W_statistic = numeric(),
    p_value_wilcox = numeric(),
    stringsAsFactors = FALSE
)

# Function to perform tests and store results (without Levene's test for variance equality)
perform_tests_no_lev <- function(df, name) {
    cat("### Performing tests for dataframe:", name, "###\n")
    
    # Filter only 'Case' and 'Control' groups (removing 'UHR_NA')
    df_filtered <- df %>% filter(group %in% c("Case", "Control"))
    
    # Ensure there are exactly two levels in 'group'
    if (length(unique(df_filtered$group)) != 2) {
        cat("Error: Group variable must have exactly 2 levels\n")
        print(unique(df_filtered$group))  # Print the unique levels of 'group'
        return(NULL)  # Skip the test if there are more than two levels
    }
    
    # List of columns to test
    columns_to_test <- c("Number_of_Variants", "Number_of_Brain_Variants", "Number_of_Brain_Filter_NTPM_Variants", "Number_of_Schema_Pval_Variants", "Number_of_Bipolar_Variants", "Number_of_Sfari_Non_Syndromic_Variants", "Number_of_Schema_Or_Variants")
    
    for (column in columns_to_test) {
        # Perform Welch's t-test (doesn't assume equal variance)
        t_test <- t.test(df_filtered[[column]] ~ df_filtered$group, var.equal = FALSE)
        
        # Append t-test results to the results data frame
        t_test_results_df_no_lev <<- rbind(t_test_results_df_no_lev, data.frame(
            DataFrame = name,
            Column = column,
            t_statistic = t_test$statistic,
            df = t_test$parameter,
            p_value_ttest = t_test$p.value,
            mean_case = t_test$estimate[1],
            mean_control = t_test$estimate[2]
        ))
        
        # Wilcoxon Rank-Sum Test (non-parametric alternative)
        wilcox_test <- wilcox.test(df_filtered[[column]] ~ df_filtered$group)
        
        # Append Wilcoxon test results to the results data frame
        wilcoxon_results_df_no_lev <<- rbind(wilcoxon_results_df_no_lev, data.frame(
            DataFrame = name,
            Column = column,
            W_statistic = wilcox_test$statistic,
            p_value_wilcox = wilcox_test$p.value
        ))
        
        cat("\nT-test result for", name, " - ", column, ":\n")
        print(t_test)
        cat("\nWilcoxon Rank-Sum Test result for", name, " - ", column, ":\n")
        print(wilcox_test)
        cat("\n")
    }
}
# Perform tests on all dataframes
perform_tests_no_lev(result_df_NON_patho_pLI, "NON_patho_pLI")
perform_tests_no_lev(result_df_patho_pLI, "patho_pLI")
perform_tests_no_lev(result_df_patho, "patho")
perform_tests_no_lev(result_df_PTV_pLI, "PTV_pLI")
perform_tests_no_lev(result_df_PTV, "PTV")

# Print the results data frames
print(t_test_results_df_no_lev)
print(wilcoxon_results_df_no_lev)
