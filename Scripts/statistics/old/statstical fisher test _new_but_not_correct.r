#statstical fisher test new 
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


# LOAD PERSONAL DATA 
# 
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

    # Function to check the label for the outlier values
    check_outlier_label <- function(outlier_value, manifest_df) {
        # Split the outlier values by ',' and trim whitespace
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

        # If there's more than one label, return the unique ones
        unique_labels <- unique(labels)
        
        # If there's both "case" and "control", return "mixed"; otherwise return the label
        if (length(unique_labels) > 1) {
            return("mixed")
        } else if (length(unique_labels) == 1) {
            return(unique_labels)
        } else {
            return(NA) # If no label found
        }
    }

    # DO the UHR NA CONTROL

        # Function to clean UHR_NA from a dataframe
        clean_uh_na <- function(df) {
        # Step 1: Remove rows where sample_label contains "UHR_NA"
        df <- df %>%
            filter(!grepl("UHR_NA", sample_label))
        
        # Step 2: Remove "UHR_NA" and clean up sample_label_2
        df <- df %>%
            mutate(
            # Remove "UHR_NA" and any surrounding commas/spaces
            sample_label_2 = gsub("\\s*,?\\s*UHR_NA\\s*,?\\s*", ",", sample_label_2),
            # Remove leading/trailing commas and spaces
            sample_label_2 = gsub("^,|,$", "", sample_label_2),
            # Replace multiple commas with a single comma
            sample_label_2 = gsub("\\s*,+\\s*", ",", sample_label_2),
            # Trim leading/trailing spaces
            sample_label_2 = trimws(sample_label_2)
            )
        
        return(df)
        }

        derive_group_label <- function(sample_label_2) {
        # Split the sample_label_2 by commas and trim whitespace
        labels <- trimws(strsplit(sample_label_2, ",")[[1]])
        
        # Determine the group based on the labels
        if (all(labels == "Case")) {
            return("Case")
        } else if (all(labels == "Control")) {
            return("Control")
        } else if (any(labels == "Case") && any(labels == "Control")) {
            return("mixed")
        } else {
            return(NA)  # If no valid labels are found
        }
        }


        # Apply the cleaning function to all dataframes in pred_df_list
        #pred_df_list_cleaned <- lapply(pred_df_list, clean_uh_na)



    # Apply the get_outlier_labels function to the dataframes
    pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_2 <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_PTV_pLI_VEP_filtered$sample_label_2 <- sapply(pred_PTV_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_PTV_VEP_filtered$sample_label_2 <- sapply(pred_PTV_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)

    #other label
    pred_NON_patho_MPC_pLI_VEP_filtered$sample_label <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_PTV_pLI_VEP_filtered$sample_label <- sapply(pred_PTV_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_PTV_VEP_filtered$sample_label <- sapply(pred_PTV_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)

    #apply clean uhna and derive group label 
    pred_NON_patho_MPC_pLI_VEP_filtered <- clean_uh_na(pred_NON_patho_MPC_pLI_VEP_filtered)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered <- clean_uh_na(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered <- clean_uh_na(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered)
    pred_PTV_pLI_VEP_filtered <- clean_uh_na(pred_PTV_pLI_VEP_filtered)
    pred_PTV_VEP_filtered <- clean_uh_na(pred_PTV_VEP_filtered)

    pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_3 <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_2,derive_group_label)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_3 <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_2,derive_group_label)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_3 <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_2,derive_group_label)
    pred_PTV_pLI_VEP_filtered$sample_label_3 <- sapply(pred_PTV_pLI_VEP_filtered$sample_label_2,derive_group_label)
    pred_PTV_VEP_filtered$sample_label_3 <- sapply(pred_PTV_VEP_filtered$sample_label_2,derive_group_label)



# CREATE DATAFRAME TO BE USED FOR STATISTICS 
# Define a function to process each dataframe and return the result
# This function processes the data to separate rows, clean up sample IDs, and filter based on specific gene sets
process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or) {
  
  expanded_df <- df %>%
    separate_rows(SAMPLES, sep = ",") %>%  # Split SAMPLES only
    mutate(
      sample_id = trimws(SAMPLES),
      group = sample_label_3  # Use pre-computed labels
    ) %>%
    filter(group %in% c("Case", "Control")) %>%  # Remove mixed/NA
    distinct(SYMBOL, POS, sample_id, .keep_all = TRUE)  # Unique variant/sample
  
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
#!statitics


#!  VARIANTS CONFRONATION
# fisher variants 
#Question: "Do specific gene sets (e.g., brain genes) have a different proportion of variants in cases vs. controls?"

    perform_variant_confrontation <- function(result_df, row) {
    # Input validation
    if (any(result_df[row, c("Case", "Control")] < 0) || 
        result_df[1, "Case"] < result_df[row, "Case"] || 
        result_df[1, "Control"] < result_df[row, "Control"]) {
        stop("Invalid counts: Negative values or subgroup exceeds total group")
    }

    # Build contingency table
    variant_counts <- matrix(
        c(
        result_df[row, "Case"], 
        result_df[1, "Case"] - result_df[row, "Case"],
        result_df[row, "Control"], 
        result_df[1, "Control"] - result_df[row, "Control"]
        ),
        nrow = 2,
        dimnames = list(
        Group = c("Cases", "Controls"),
        Variant = c("In Group", "Not In Group")
        )
    )

    # Auto-select test based on expected counts
    expected <- chisq.test(variant_counts)$expected
    if (any(expected < 5)) {
        test <- fisher.test(variant_counts)
        method <- "Fisher"
    } else {
        test <- chisq.test(variant_counts)
        method <- "ChiSq"
    }

    return(list(
        contingency_table = variant_counts,
        p_value = test$p.value,
        odds_ratio = ifelse(method == "Fisher", test$estimate, (variant_counts[1,1]/variant_counts[1,2])/(variant_counts[2,1]/variant_counts[2,2])),
        method = method
    ))
    }

    perform_variant_confrontation_tests <- function(result_df, test_name) {
    rows <- c(3, 5, 7, 9, 11, 13)
    row_names <- c("brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants",
                    "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")

    results <- lapply(seq_along(rows), function(i) {
        res <- perform_variant_confrontation(result_df, rows[i])
        data.frame(
        Test = test_name,  # Explicit test name
        Row = row_names[i],
        Cases_In = res$contingency_table[1,1],
        Cases_Out = res$contingency_table[1,2],
        Controls_In = res$contingency_table[2,1],
        Controls_Out = res$contingency_table[2,2],
        OR = res$odds_ratio,
        p_value = res$p_value,
        Method = res$method
        )
    })
    
    do.call(rbind, results)
    }

    results <- list(
    NON_patho_pLI = perform_variant_confrontation_tests(result_df_NON_patho_pLI, "NON_patho_pLI"),
    patho_pLI = perform_variant_confrontation_tests(result_df_patho_pLI, "patho_pLI"),
    patho = perform_variant_confrontation_tests(result_df_patho, "patho"),
    PTV_pLI = perform_variant_confrontation_tests(result_df_PTV_pLI, "PTV_pLI"),
    PTV = perform_variant_confrontation_tests(result_df_PTV, "PTV")
    )

    final_results <- do.call(rbind, results)
    # After creating final_results
    final_results$p_adj_global <- p.adjust(final_results$p_value, method = "fdr")
    final_results <- final_results %>%
    group_by(Test) %>%  # Group by dataset (e.g., NON_patho_pLI)
    mutate(p_adj_per_dataset = p.adjust(p_value, method = "fdr")) %>%
    ungroup()



#! Perform Fisher's exact test for individuals carries at least 1 varaint
#Question:
"Are individuals with ≥1 variant in a gene set more common in cases vs. controls?"
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

# fisher individual 
    perform_fisher_test_per_ind <- function(case_count, control_count, total_cases, total_controls) {
    # Input validation
    if (case_count > total_cases || control_count > total_controls) {
        stop("Carrier counts exceed group totals")
    }
    
    # Build contingency table with fail-safe
    a <- case_count
    c <- control_count
    b <- total_cases - a
    d <- total_controls - c
    
    contingency_table <- matrix(
        c(a, b, c, d), 
        nrow = 2,
        dimnames = list(
        Group = c("Cases", "Controls"),
        Carrier = c("Yes", "No")
        )
    )

    # Handle zero margins
    if (any(contingency_table < 0)) {
        return(list(
        error = "Negative counts in contingency table",
        table = contingency_table
        ))
    }

    test <- fisher.test(contingency_table)
    
    return(list(
        p_value = test$p.value,
        odds_ratio = test$estimate,
        ci_low = test$conf.int[1],
        ci_high = test$conf.int[2],
        table = contingency_table
    ))
    }

    perform_fisher_tests_for_all_rows <- function(result_df) {
    rows <- c(2, 4, 6, 8, 10, 12, 14)
    row_names <- c("filtered", "brain_variants", "brain_filter_ntpm_variants",
                    "schema_pval_variants", "bipolar_variants", 
                    "sfari_non_syndromic_variants", "schema_or_variants")

    results <- lapply(seq_along(rows), function(i) {
        res <- perform_fisher_test_per_ind(
        result_df$Case[rows[i]],
        result_df$Control[rows[i]],
        sum(result_df$Case[c(1, rows[i])]),  # Total cases
        sum(result_df$Control[c(1, rows[i])]) # Total controls
        )
        
        data.frame(
        Row = row_names[i],
        Case_Carriers = res$table[1,1],
        Case_NonCarriers = res$table[1,2],
        Control_Carriers = res$table[2,1],
        Control_NonCarriers = res$table[2,2],
        OR = round(res$odds_ratio, 2),
        CI = paste0("[", round(res$ci_low,2), ", ", round(res$ci_high,2), "]"),
        p_value = format.pval(res$p_value, eps = 0.0001)
        )
    })
    
    do.call(rbind, results)
    }


    # Apply the function to each result dataframe and compile the results
    fisher_results_ind <- list(
    NON_patho_pLI = perform_fisher_tests_for_all_rows(result_df_NON_patho_pLI),
    patho_pLI = perform_fisher_tests_for_all_rows(result_df_patho_pLI),
    patho = perform_fisher_tests_for_all_rows(result_df_patho),
    PTV_pLI = perform_fisher_tests_for_all_rows(result_df_PTV_pLI),
    PTV = perform_fisher_tests_for_all_rows(result_df_PTV)
    )

    fisher_results_ind_df <- do.call(rbind, lapply(names(fisher_results_ind), function(test_name) {
    df <- fisher_results_ind[[test_name]]
    df$Test <- test_name  # Add Test column
    df
    }))

    # Global FDR across all tests
    fisher_results_ind_df$p_adj_global <- p.adjust(fisher_results_ind_df$p_value, method = "fdr")

    # Per-Test FDR
    fisher_results_ind_df <- fisher_results_ind_df %>%
    group_by(Test) %>%
    mutate(p_adj_per_test = p.adjust(p_value, method = "fdr")) %>%
    ungroup()


#  rate comparison 
#Question: "What is the proportion of individuals carrying variants in a gene set (cases vs. controls)?"
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
        rows <- c(1, 4, 7, 10, 13, 16, 19)  # Rows to be tested
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
    

# fisher case contro 
# Question:"Does a dataset (e.g., NON_patho_pLI) contain more variants originating from cases than expected by chance?"
# "Do these variant categories contain significantly more case-associated variants than expected by chance?"
# DO the UHR NA CONTROL

    # Function to clean UHR_NA from a dataframe
    clean_uh_na <- function(df) {
    # Step 1: Remove rows where sample_label contains "UHR_NA"
    df <- df %>%
        filter(!grepl("UHR_NA", sample_label))
    
    # Step 2: Remove "UHR_NA" and clean up sample_label_2
    df <- df %>%
        mutate(
        # Remove "UHR_NA" and any surrounding commas/spaces
        sample_label_2 = gsub("\\s*,?\\s*UHR_NA\\s*,?\\s*", ",", sample_label_2),
        # Remove leading/trailing commas and spaces
        sample_label_2 = gsub("^,|,$", "", sample_label_2),
        # Replace multiple commas with a single comma
        sample_label_2 = gsub("\\s*,+\\s*", ",", sample_label_2),
        # Trim leading/trailing spaces
        sample_label_2 = trimws(sample_label_2)
        )
    
    return(df)
    }

    derive_group_label <- function(sample_label_2) {
    # Split the sample_label_2 by commas and trim whitespace
    labels <- trimws(strsplit(sample_label_2, ",")[[1]])
    
    # Determine the group based on the labels
    if (all(labels == "Case")) {
        return("Case")
    } else if (all(labels == "Control")) {
        return("Control")
    } else if (any(labels == "Case") && any(labels == "Control")) {
        return("mixed")
    } else {
        return(NA)  # If no valid labels are found
    }
    }


    # Apply the cleaning function to all dataframes in pred_df_list
    pred_df_list_cleaned <- lapply(pred_df_list, clean_uh_na)

# without mixed 
    # Create a list of your dataframes
    pred_df_list <- list(
    NON_patho_pLI = pred_NON_patho_MPC_pLI_VEP_filtered,
    patho_pLI = pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered,
    patho = pred_patho_MPC_ALPHAMISSENSE_VEP_filtered,
    PTV_pLI = pred_PTV_pLI_VEP_filtered,
    PTV = pred_PTV_VEP_filtered
    )

    # Function to create contingency table
    create_contingency_table <- function(df) {
    # Ensure variants are unique based on CHROM, POS, REF, ALT
    df <- df %>%
        distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)
    
    # Get group labels for each variant
    df$group <- sapply(df$sample_label_2, derive_group_label)
    
    # Calculate counts
    total_case = sum(df$group %in% c("Case"), na.rm = TRUE)
    total_control = sum(df$group %in% c("Control"), na.rm = TRUE)
    total_variants = nrow(df)
    
    # Create contingency table
    other_case = total_variants - total_case
    other_control = total_variants - total_control
    
    contingency_table <- matrix(
    c(total_case, other_case,  # Row 1: Case
      total_control, other_control),  # Row 2: Control
    nrow = 2,
    byrow = TRUE,  # Fill the matrix row-wise
    dimnames = list(
      Group = c("Case", "Control"),
      Variant = c("In Group", "Not In Group")
    )
  )
    # Return counts along with the contingency table
    return(list(
        contingency_table = contingency_table,
        total_variants = total_variants,
        total_case = total_case,
        total_control = total_control
    ))
    }

    # Generate contingency tables and counts for all dataframes
    contingency_results <- lapply(pred_df_list, create_contingency_table)

    # Function to perform Fisher's exact test on a contingency table
    perform_fisher_test <- function(contingency_table) {
    # Perform Fisher's exact test
    test_result <- fisher.test(contingency_table)
    
    # Extract results
    list(
        p_value = test_result$p.value,
        odds_ratio = test_result$estimate,
        conf_int = test_result$conf.int
    )
    }

    # Apply Fisher's test to all contingency tables
    fisher_results <- lapply(contingency_results, function(res) {
    list(
        test_result = perform_fisher_test(res$contingency_table),
        total_variants = res$total_variants,
        total_case = res$total_case,
        total_control = res$total_control
    )
    })

    # Combine results into a dataframe
    fisher_results_df <- do.call(rbind, lapply(names(fisher_results), function(test_name) {
    res <- fisher_results[[test_name]]
    data.frame(
        Test = test_name,
        Total_Variants = res$total_variants,
        Case_Variants = res$total_case,
        Control_Variants = res$total_control,
        OR = res$test_result$odds_ratio,
        CI_Low = res$test_result$conf_int[1],
        CI_High = res$test_result$conf_int[2],
        p_value = res$test_result$p_value
    )
    }))

    # Add FDR correction for multiple testing
    fisher_results_df$p_adj <- p.adjust(fisher_results_df$p_value, method = "fdr")

    # Print the results
    print(fisher_results_df)

# wRong
# with mixed

    create_contingency_table <- function(df) {
    # Ensure variants are unique based on CHROM, POS, REF, ALT
    df <- df %>%
        distinct(CHROM, POS, REF, ALT, .keep_all = TRUE)
    
    # Get group labels for each variant
    df$group <- sapply(df$SAMPLES, check_outlier_label, manifest_df = manifest_correct)
    
    # Calculate counts
    total_case = sum(df$group %in% c("Case", "mixed"), na.rm = TRUE)  # Include mixed
    total_control = sum(df$group %in% c("Control", "mixed"), na.rm = TRUE)  # Include mixed
    total_variants = nrow(df)
    
    # Create contingency table
    other_case = total_variants - total_case
    other_control = total_variants - total_control
    
    contingency_table <- matrix(
        c(total_case, other_case,
        total_control, other_control),
        nrow = 2,
        dimnames = list(
        Group = c("Case", "Control"),
        Variant = c("In Group", "Other")
        )
    )
    
    return(contingency_table)
    }

    # Generate contingency tables and counts for all dataframes
    contingency_results <- lapply(pred_df_list, create_contingency_table)

    # Function to perform Fisher's exact test on a contingency table
    perform_fisher_test <- function(contingency_table) {
    # Perform Fisher's exact test
    test_result <- fisher.test(contingency_table)
    
    # Extract results
    list(
        p_value = test_result$p.value,
        odds_ratio = test_result$estimate,
        conf_int = test_result$conf.int
    )
    }

    # Apply Fisher's test to all contingency tables
    fisher_results_mixed <- lapply(contingency_results, function(res) {
    list(
        test_result = perform_fisher_test(res$contingency_table),
        total_variants = res$total_variants,
        total_case = res$total_case,
        total_control = res$total_control
    )
    })

    # Combine results into a dataframe
    fisher_results_mixed_df <- do.call(rbind, lapply(names(fisher_results_mixed), function(test_name) {
    res <- fisher_results_mixed[[test_name]]
    data.frame(
        Test = test_name,
        Total_Variants = res$total_variants,
        Case_Variants = res$total_case,
        Control_Variants = res$total_control,
        OR = res$test_result$odds_ratio,
        CI_Low = res$test_result$conf_int[1],
        CI_High = res$test_result$conf_int[2],
        p_value = res$test_result$p_value
    )
    }))

    # Add FDR correction for multiple testing
    fisher_results_mixed_df$p_adj <- p.adjust(fisher_results_mixed_df$p_value, method = "fdr")

    # Print the results
    print(fisher_results_mixed_df)