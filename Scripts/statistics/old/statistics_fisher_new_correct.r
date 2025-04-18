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
library(ggplot2)
# Load my data from CSV files
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
    # Specify the folder containing your TSV files
    folder_path <- "/home/rachele/SNVs/results/euro"

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

# ADD LABELS 
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

     df <- df %>%
    mutate(
      # First remove all UHR_NA instances
      sample_label_2 = gsub("UHR_NA", "", sample_label_2),
      # Then clean up resulting comma patterns
      sample_label_2 = gsub(",+", ",", sample_label_2),  # Replace multiple commas with single
      sample_label_2 = gsub("^,|,$", "", sample_label_2),  # Remove leading/trailing commas
      sample_label_2 = trimws(sample_label_2),  # Trim whitespace
      # Handle cases where we're left with just a comma
      sample_label_2 = ifelse(sample_label_2 == ",", NA, sample_label_2)
    )
      # Additional check for empty strings
  df$sample_label_2[df$sample_label_2 == ""] <- NA    
    
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




# Apply the get_outlier_labels function to the dataframes
    pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_2 <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_NON_patho_MPC_VEP_filtered$sample_label_2 <- sapply(pred_NON_patho_MPC_VEP_filtered$SAMPLES, get_outlier_labels,manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_2 <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_PTV_pLI_VEP_filtered$sample_label_2 <- sapply(pred_PTV_pLI_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)
    pred_PTV_VEP_filtered$sample_label_2 <- sapply(pred_PTV_VEP_filtered$SAMPLES, get_outlier_labels, manifest_correct)

    #other label
    pred_NON_patho_MPC_pLI_VEP_filtered$sample_label <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_NON_patho_MPC_VEP_filtered$sample_label <- sapply(pred_NON_patho_MPC_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_PTV_pLI_VEP_filtered$sample_label <- sapply(pred_PTV_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
    pred_PTV_VEP_filtered$sample_label <- sapply(pred_PTV_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)

    #apply clean uhna and derive group label 
    pred_NON_patho_MPC_pLI_VEP_filtered <- clean_uh_na(pred_NON_patho_MPC_pLI_VEP_filtered)
    pred_NON_patho_MPC_VEP_filtered <- clean_uh_na(pred_NON_patho_MPC_VEP_filtered)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered <- clean_uh_na(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered <- clean_uh_na(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered)
    pred_PTV_pLI_VEP_filtered <- clean_uh_na(pred_PTV_pLI_VEP_filtered)
    pred_PTV_VEP_filtered <- clean_uh_na(pred_PTV_VEP_filtered)

    pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_3 <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_2,derive_group_label)
    pred_NON_patho_MPC_VEP_filtered$sample_label_3 <- sapply(pred_NON_patho_MPC_VEP_filtered$sample_label_2,derive_group_label)
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_3 <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_2,derive_group_label)
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_3 <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_2,derive_group_label)
    pred_PTV_pLI_VEP_filtered$sample_label_3 <- sapply(pred_PTV_pLI_VEP_filtered$sample_label_2,derive_group_label)
    pred_PTV_VEP_filtered$sample_label_3 <- sapply(pred_PTV_VEP_filtered$sample_label_2,derive_group_label)

    pred_NON_patho_MPC_pLI_VEP_filtered <-pred_NON_patho_MPC_pLI_VEP_filtered %>% mutate(count = sapply(strsplit(SAMPLES, ","), length))
    pred_NON_patho_MPC_VEP_filtered <-pred_NON_patho_MPC_VEP_filtered%>% mutate(count = sapply(strsplit(SAMPLES, ","), length))
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered <-pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered %>% mutate(count = sapply(strsplit(SAMPLES, ","), length))
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered <-pred_patho_MPC_ALPHAMISSENSE_VEP_filtered%>% mutate(count = sapply(strsplit(SAMPLES, ","), length))
    pred_PTV_pLI_VEP_filtered <-pred_PTV_pLI_VEP_filtered%>% mutate(count = sapply(strsplit(SAMPLES, ","), length))
    pred_PTV_VEP_filtered <-pred_PTV_VEP_filtered%>% mutate(count = sapply(strsplit(SAMPLES, ","), length))






# manifest for UHR_NA
    UHR_NA_SAMPLE_IDS <- read.csv("~/UHR_NA_SAMPLE_IDS.csv", sep="")
    #example S36827

# CREATE DATAFRAME TO BE USED FOR STATISTICS 
    # Define a function to process each dataframe and return the result
    # This function processes the data to separate rows, clean up sample IDs, and filter based on specific gene sets
    process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE) {

        if(private){
        df <- df %>% filter(count == 1)  # Only filter if private=TRUE
    }


    #dataframe to count individuals 
    expanded_df <- df %>%
        separate_rows(SAMPLES, sep = ",") %>%  # Split SAMPLES only
        mutate(
        sample_id = trimws(SAMPLES),
        group = sample_label_3  # Use pre-computed labels
        ) %>%
        filter(group %in% c("Case", "Control")) %>%  # Remove mixed/NA
        filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%  # Remove UHR_NA
        distinct(SYMBOL, CHROM,REF,ALT, POS, sample_id, .keep_all = TRUE)  # Unique variant/sample

    #dataframe for count variants 
    variant_df <- df %>% 
    mutate(group = sample_label_3) %>%  # Use pre-computed labels
    filter(group %in% c("Case", "Control")) %>%  # Remove mixed/NA
    distinct(SYMBOL, CHROM,REF,ALT, POS , .keep_all = TRUE)  # Unique variant/sample 
    
    # Filter based on the genes in the brain, brain filter_ntpm, and additional gene sets
    expanded_df <- expanded_df %>%
        mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
            brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
            schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
            bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
            sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
            schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))

    variant_df <- variant_df %>%
            mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
            brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
            schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
            bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
            sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
            schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))
    
    num_non_pathogenic_pLi_variants <- variant_df %>%
        group_by(group) %>%
        summarize(num_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_non_pathogenic_pLi <- expanded_df %>%
        group_by(group) %>%
        summarize(unique_individuals = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_non_pathogenic_pLi <- expanded_df %>%
        group_by(group) %>%  
        summarize(num_individuals = n(), .groups = "drop")
    
    num_brain_variants <- variant_df %>%
        filter(brain == 1) %>%
        group_by(group) %>%
        summarize(num_brain_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_brain_variants <- expanded_df %>%
        filter(brain == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_brain = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_brain_variants <- expanded_df %>%
        filter(brain == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_brain = n(), .groups = "drop")  
    
    num_brain_filter_ntpm_variants <- variant_df %>%
        filter(brain_filter_ntpm == 1) %>%
        group_by(group) %>%
        summarize(num_brain_filter_ntpm_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_brain_filter_ntpm_variants <- expanded_df %>%
        filter(brain_filter_ntpm == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_brain_filter_ntpm = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_brain_filter_ntpm_variants <- expanded_df %>%
    filter(brain_filter_ntpm == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_brain_filter_ntpm = n(), .groups = "drop")
    
    num_schema_pval_variants <- variant_df %>%
        filter(schema_pval == 1) %>%
        group_by(group) %>%
        summarize(num_schema_pval_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_schema_pval_variants <- expanded_df %>%
        filter(schema_pval == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_schema_pval = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_schema_pval_variants <- expanded_df %>%
        filter(schema_pval == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_schema_pval = n(), .groups = "drop")  
    
    num_bipolar_variants <- variant_df %>%
        filter(bipolar == 1) %>%
        group_by(group) %>%
        summarize(num_bipolar_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_bipolar_variants <- expanded_df %>%
        filter(bipolar == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_bipolar = n_distinct(sample_id), .groups = "drop")
    
    num_individuals_with_bipolar_variants <- expanded_df %>%
        filter(bipolar == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_bipolar = n(), .groups = "drop")
    
    num_sfari_non_syndromic_variants <- variant_df %>%
        filter(sfari_non_syndromic == 1) %>%
        group_by(group) %>%
        summarize(num_sfari_non_syndromic_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_sfari_non_syndromic_variants <- expanded_df %>%
        filter(sfari_non_syndromic == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_sfari_non_syndromic = n_distinct(sample_id), .groups = "drop")

    num_individuals_with_sfari_non_syndromic_variants <- expanded_df %>%
        filter(sfari_non_syndromic == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_sfari_non_syndromic = n(), .groups = "drop")  
    
    num_schema_or_variants <- variant_df %>%
        filter(schema_or == 1) %>%
        group_by(group) %>%
        summarize(num_schema_or_variants = n(), .groups = "drop")
    
    num_individuals_unique_with_schema_or_variants <- expanded_df %>%
        filter(schema_or == 1) %>%
        group_by(group) %>%
        summarize(unique_individuals_schema_or = n_distinct(sample_id), .groups = "drop")
    
    num_individuals_with_schema_or_variants <- expanded_df %>%
        filter(schema_or == 1) %>%
        group_by(group) %>%
        summarize(num_individuals_schema_or = n(), .groups = "drop")
    
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
    groups <- c("Case", "Control")
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

    num_individuals_unique_with_non_pathogenic_pLi <- num_individuals_unique_with_non_pathogenic_pLi %>% complete(group = groups, fill = list(unique_individuals = 0))
    num_individuals_unique_with_brain_variants <- num_individuals_unique_with_brain_variants %>% complete(group = groups, fill = list(unique_individuals_brain = 0))
    num_individuals_unique_with_brain_filter_ntpm_variants <- num_individuals_unique_with_brain_filter_ntpm_variants %>% complete(group = groups, fill = list(unique_individuals_brain_filter_ntpm = 0))
    num_individuals_unique_with_schema_pval_variants <- num_individuals_unique_with_schema_pval_variants %>% complete(group = groups, fill = list(unique_individuals_schema_pval = 0))
    num_individuals_unique_with_bipolar_variants <- num_individuals_unique_with_bipolar_variants %>% complete(group = groups, fill = list(unique_individuals_bipolar = 0))
    num_individuals_unique_with_sfari_non_syndromic_variants <- num_individuals_unique_with_sfari_non_syndromic_variants %>% complete(group = groups, fill = list(unique_individuals_sfari_non_syndromic = 0))
    num_individuals_unique_with_schema_or_variants <- num_individuals_unique_with_schema_or_variants %>% complete(group = groups, fill = list(unique_individuals_schema_or = 0))


    # Define the row names based on the dataframe's name
    row_names <- c(
        paste("Number of", name ,"variants"),
        paste("Number of individuals with",name ,"variants"),
        paste("Number of individuals unique with",name ,"variants"),
        paste("Number of brain", name ,"variants"),
        paste("Number of individuals with brain",name ,"variants"),
        paste("Number of individuals unique with brain",name ,"variants"),
        paste("Number of brain filter_ntpm", name ,"variants"),
        paste("Number of individuals with brain filter_ntpm",name ,"variants"),
        paste("Number of individuals unique with brain filter_ntpm",name ,"variants"),
        paste("Number of schema_pval", name ,"variants"),
        paste("Number of individuals with schema_pval",name ,"variants"),
        paste("Number of individuals unique with schema_pval",name ,"variants"),
        paste("Number of bipolar", name ,"variants"),
        paste("Number of individuals with bipolar",name ,"variants"),
        paste("Number of individuals unique with bipolar",name ,"variants"),
        paste("Number of sfari_non_syndromic", name ,"variants"),
        paste("Number of individuals with sfari_non_syndromic",name ,"variants"),
        paste("Number of individuals unique with sfari_non_syndromic",name ,"variants"),
        paste("Number of schema_or", name ,"variants"),
        paste("Number of individuals with schema_or",name ,"variants"),
        paste("Number of individuals unique with schema_or",name ,"variants"),
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
            num_individuals_with_non_pathogenic_pLi$num_individuals[num_individuals_with_non_pathogenic_pLi$group == "Case"],  # Corrected column
            num_individuals_unique_with_non_pathogenic_pLi$unique_individuals[num_individuals_unique_with_non_pathogenic_pLi$group == "Case"],
            num_brain_variants$num_brain_variants[num_brain_variants$group == "Case"],
            num_individuals_with_brain_variants$num_individuals_brain[num_individuals_with_brain_variants$group == "Case"],  # Corrected column
            num_individuals_unique_with_brain_variants$unique_individuals_brain[num_individuals_unique_with_brain_variants$group == "Case"],
            num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "Case"],
            num_individuals_with_brain_filter_ntpm_variants$num_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "Case"],  # Corrected column
            num_individuals_unique_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_unique_with_brain_filter_ntpm_variants$group == "Case"],
            num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "Case"],
            num_individuals_with_schema_pval_variants$num_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "Case"],  # Corrected column
            num_individuals_unique_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_unique_with_schema_pval_variants$group == "Case"],
            num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "Case"],
            num_individuals_with_bipolar_variants$num_individuals_bipolar[num_individuals_with_bipolar_variants$group == "Case"],  # Corrected column
            num_individuals_unique_with_bipolar_variants$unique_individuals_bipolar[num_individuals_unique_with_bipolar_variants$group == "Case"],
            num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "Case"],
            num_individuals_with_sfari_non_syndromic_variants$num_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "Case"],  # Corrected column
            num_individuals_unique_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_unique_with_sfari_non_syndromic_variants$group == "Case"],
            num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "Case"],
            num_individuals_with_schema_or_variants$num_individuals_schema_or[num_individuals_with_schema_or_variants$group == "Case"],  # Corrected column
            num_individuals_unique_with_schema_or_variants$unique_individuals_schema_or[num_individuals_unique_with_schema_or_variants$group == "Case"],
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
            num_individuals_with_non_pathogenic_pLi$num_individuals[num_individuals_with_non_pathogenic_pLi$group == "Control"],  # Corrected column
            num_individuals_unique_with_non_pathogenic_pLi$unique_individuals[num_individuals_unique_with_non_pathogenic_pLi$group == "Control"],
            num_brain_variants$num_brain_variants[num_brain_variants$group == "Control"],
            num_individuals_with_brain_variants$num_individuals_brain[num_individuals_with_brain_variants$group == "Control"],  # Corrected column
            num_individuals_unique_with_brain_variants$unique_individuals_brain[num_individuals_unique_with_brain_variants$group == "Control"],
            num_brain_filter_ntpm_variants$num_brain_filter_ntpm_variants[num_brain_filter_ntpm_variants$group == "Control"],
            num_individuals_with_brain_filter_ntpm_variants$num_individuals_brain_filter_ntpm[num_individuals_with_brain_filter_ntpm_variants$group == "Control"],  # Corrected column
            num_individuals_unique_with_brain_filter_ntpm_variants$unique_individuals_brain_filter_ntpm[num_individuals_unique_with_brain_filter_ntpm_variants$group == "Control"],
            num_schema_pval_variants$num_schema_pval_variants[num_schema_pval_variants$group == "Control"],
            num_individuals_with_schema_pval_variants$num_individuals_schema_pval[num_individuals_with_schema_pval_variants$group == "Control"],  # Corrected column
            num_individuals_unique_with_schema_pval_variants$unique_individuals_schema_pval[num_individuals_unique_with_schema_pval_variants$group == "Control"],
            num_bipolar_variants$num_bipolar_variants[num_bipolar_variants$group == "Control"],
            num_individuals_with_bipolar_variants$num_individuals_bipolar[num_individuals_with_bipolar_variants$group == "Control"],  # Corrected column
            num_individuals_unique_with_bipolar_variants$unique_individuals_bipolar[num_individuals_unique_with_bipolar_variants$group == "Control"],
            num_sfari_non_syndromic_variants$num_sfari_non_syndromic_variants[num_sfari_non_syndromic_variants$group == "Control"],
            num_individuals_with_sfari_non_syndromic_variants$num_individuals_sfari_non_syndromic[num_individuals_with_sfari_non_syndromic_variants$group == "Control"],  # Corrected column
            num_individuals_unique_with_sfari_non_syndromic_variants$unique_individuals_sfari_non_syndromic[num_individuals_unique_with_sfari_non_syndromic_variants$group == "Control"],
            num_schema_or_variants$num_schema_or_variants[num_schema_or_variants$group == "Control"],
            num_individuals_with_schema_or_variants$num_individuals_schema_or[num_individuals_with_schema_or_variants$group == "Control"],  # Corrected column
            num_individuals_unique_with_schema_or_variants$unique_individuals_schema_or[num_individuals_unique_with_schema_or_variants$group == "Control"],
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

    # Process each dataframe
    result_df_list <- list(
    process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or),
    process_df(pred_NON_patho_MPC_VEP_filtered, "NON_patho",brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or ),
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
    result_df_NON_patho <- process_df(pred_NON_patho_MPC_VEP_filtered, "NON_patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or )
    result_df_patho_pLI <- process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
    result_df_patho <- process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
    result_df_PTV_pLI <- process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
    result_df_PTV <- process_df(pred_PTV_VEP_filtered, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)

#when there is an NA put 0 
final_result_df[is.na(final_result_df)] <-0
result_df_NON_patho_pLI[is.na(result_df_NON_patho_pLI)] <- 0
result_df_NON_patho[is.na(result_df_NON_patho)] <- 0
result_df_patho_pLI[is.na(result_df_patho_pLI)] <- 0
result_df_patho[is.na(result_df_patho)] <- 0
result_df_PTV_pLI[is.na(result_df_PTV_pLI)] <- 0
result_df_PTV[is.na(result_df_PTV)] <- 0


# fisher variants
##Question: "Do specific gene sets (e.g., brain genes) have a different proportion of variants in cases vs. controls?"
    perform_variant_confrontation <- function(result_df, row) {
        # Input validation
        if (any(result_df[row, c("Case", "Control")] < 0) || 
            result_df[1, "Case"] < result_df[row, "Case"] || 
            result_df[1, "Control"] < result_df[row, "Control"]) {
            stop("Invalid counts: Negative values or subgroup exceeds total group")
        }

        # Calculate counts
        total_case <- result_df[1, "Case"]  # Total cases (from the first row)
        total_control <- result_df[1, "Control"]  # Total controls (from the first row)
        
        case_with_variant <- result_df[row, "Case"]  # Cases with the variant
        control_with_variant <- result_df[row, "Control"]  # Controls with the variant
        
        case_without_variant <- total_case - case_with_variant  # Cases without the variant
        control_without_variant <- total_control - control_with_variant  # Controls without the variant

        # Create contingency table
        variant_counts <- matrix(
            c(
                case_with_variant, case_without_variant,  # Row 1: Cases
                control_with_variant, control_without_variant  # Row 2: Controls
            ),
            nrow = 2,
            byrow = TRUE,  # Fill the matrix row-wise
            dimnames = list(
                Group = c("Case", "Control"),
                Variant = c("In Group", "Not In Group")
            )
        )

        # Auto-select test based on expected counts
        expected <- chisq.test(variant_counts)$expected
        if (any(expected < 5)) {
            test <- fisher.test(variant_counts)
            method <- "Fisher"
            odds_ratio <- test$estimate
        } else {
            test <- chisq.test(variant_counts)
            method <- "ChiSq"
            odds_ratio <- (variant_counts[1, 1] / variant_counts[1, 2]) / (variant_counts[2, 1] / variant_counts[2, 2])
        }

        return(list(
            contingency_table = variant_counts,
            p_value = test$p.value,
            odds_ratio = odds_ratio,
            method = method,
            total_case = total_case,
            total_control = total_control,
            case_with_variant = case_with_variant,
            control_with_variant = control_with_variant,
            case_without_variant = case_without_variant,
            control_without_variant = control_without_variant
        ))
    }

    perform_variant_confrontation_tests <- function(result_df, test_name) {
        # Define rows to test and their corresponding names
        rows <- c(4, 7, 10, 13, 16, 19)
        row_names <- c("brain_variants", "brain_filter_ntpm_variants", "schema_pval_variants",
                    "bipolar_variants", "sfari_non_syndromic_variants", "schema_or_variants")

        # Perform tests for each row
        results <- lapply(seq_along(rows), function(i) {
            res <- perform_variant_confrontation(result_df, rows[i])
            data.frame(
                Test = test_name,  # Explicit test name
                Row = row_names[i],
                Cases_In = res$contingency_table[1, 1],
                Cases_Out = res$contingency_table[1, 2],
                Controls_In = res$contingency_table[2, 1],
                Controls_Out = res$contingency_table[2, 2],
                OR = res$odds_ratio,
                p_value = res$p_value,
                Method = res$method,
                Total_Cases = res$total_case,
                Total_Controls = res$total_control,
                Cases_With_Variant = res$case_with_variant,
                Controls_With_Variant = res$control_with_variant,
                Cases_Without_Variant = res$case_without_variant,
                Controls_Without_Variant = res$control_without_variant
            )
        })
        
        # Combine results into a single dataframe
        do.call(rbind, results)
    }

    # Perform tests for all datasets
    results <- list(
        NON_patho_pLI = perform_variant_confrontation_tests(result_df_NON_patho_pLI, "NON_patho_pLI"),
        NON_patho = perform_variant_confrontation_tests(result_df_NON_patho,"NON_patho"),
        patho_pLI = perform_variant_confrontation_tests(result_df_patho_pLI, "patho_pLI"),
        patho = perform_variant_confrontation_tests(result_df_patho, "patho"),
        PTV_pLI = perform_variant_confrontation_tests(result_df_PTV_pLI, "PTV_pLI"),
        PTV = perform_variant_confrontation_tests(result_df_PTV, "PTV")
    )

    # Combine all results into a single dataframe
    final_results <- do.call(rbind, results)

    # Apply multiple testing corrections
    final_results$p_adj_global <- p.adjust(final_results$p_value, method = "fdr")
    final_results <- final_results %>%
        group_by(Test) %>%  # Group by dataset (e.g., NON_patho_pLI)
        mutate(p_adj_per_test = p.adjust(p_value, method = "fdr")) %>%
        ungroup()

    # View final results
    print(final_results)

# fisher individuals 
#Question:"Are individuals with ≥1 variant in a gene set more common in cases vs. controls?"
    # Function to perform Fisher's exact test for a single variant type
    perform_fisher_test_per_ind <- function(case_count, control_count, total_cases, total_controls) {
        # Input validation
        if (case_count > total_cases || control_count > total_controls) {
            stop("Carrier counts exceed group totals")
        }
        
        # Build contingency table
        a <- case_count  # Cases with the variant
        c <- control_count  # Controls with the variant
        b <- total_cases - a  # Cases without the variant
        d <- total_controls - c  # Controls without the variant
        
        contingency_table <- matrix(
            c(a, b, c, d), 
            nrow = 2,
            byrow = TRUE,  # Fill the matrix row-wise
            dimnames = list(
                Group = c("Case", "Control"),
                Variant = c("Yes", "No")
            )
        )

        # Handle zero margins
        if (any(contingency_table < 0)) {
            return(list(
                error = "Negative counts in contingency table",
                table = contingency_table
            ))
        }

        # Perform Fisher's exact test
        test <- fisher.test(contingency_table)
        
        return(list(
            p_value = test$p.value,
            odds_ratio = test$estimate,
            ci_low = test$conf.int[1],
            ci_high = test$conf.int[2],
            table = contingency_table
        ))
    }

    # Function to perform Fisher's tests for all rows in the result dataframe
    perform_fisher_tests_for_all_rows <- function(result_df) {
        # Define rows to test and their corresponding names
        rows <- c(3, 6, 9, 12, 15, 18, 21)
        row_names <- c("filtered", "brain_variants", "brain_filter_ntpm_variants",
                    "schema_pval_variants", "bipolar_variants", 
                    "sfari_non_syndromic_variants", "schema_or_variants")

        # Perform tests for each row
        results <- lapply(seq_along(rows), function(i) {
            res <- perform_fisher_test_per_ind(
                result_df$Case[rows[i]],
                result_df$Control[rows[i]],
                302,  # Total cases (fixed value)
                75    # Total controls (fixed value)
            )
            
            data.frame(
                Row = row_names[i],
                Case_Carriers = res$table[1, 1],
                Case_NonCarriers = res$table[1, 2],
                Control_Carriers = res$table[2, 1],
                Control_NonCarriers = res$table[2, 2],
                OR = round(res$odds_ratio, 2),
                CI = paste0("[", round(res$ci_low, 2), ", ", round(res$ci_high, 2), "]"),
                p_value = format.pval(res$p_value, eps = 0.0001)
            )
        })
        
        do.call(rbind, results)
    }

    # Apply the function to each result dataframe and compile the results
    fisher_results_ind <- list(
        NON_patho_pLI = perform_fisher_tests_for_all_rows(result_df_NON_patho_pLI),
        NON_patho = perform_fisher_tests_for_all_rows(result_df_NON_patho),
        patho_pLI = perform_fisher_tests_for_all_rows(result_df_patho_pLI),
        patho = perform_fisher_tests_for_all_rows(result_df_patho),
        PTV_pLI = perform_fisher_tests_for_all_rows(result_df_PTV_pLI),
        PTV = perform_fisher_tests_for_all_rows(result_df_PTV)
    )

    # Combine all results into a single dataframe
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

    # View final results
    print(fisher_results_ind_df)

# fisher case control 
#question = is the number of variance between case and control signifcant diffrent 


    # Create a list of your dataframes
    pred_df_list <- list(
    NON_patho_pLI = pred_NON_patho_MPC_pLI_VEP_filtered,
    NON_patho = pred_NON_patho_MPC_VEP_filtered,
    patho_pLI = pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered,
    patho = pred_patho_MPC_ALPHAMISSENSE_VEP_filtered,
    PTV_pLI = pred_PTV_pLI_VEP_filtered,
    PTV = pred_PTV_VEP_filtered
    )

    create_contingency_table <- function(df) {
    # Split variants by sample and group
    expanded_df <- df %>%
        separate_rows(SAMPLES, sep = ",") %>%
        mutate(
        sample_id = trimws(SAMPLES),
        group = sapply(sample_label_2, derive_group_label)
        ) %>%
        filter(group %in% c("Case", "Control")) %>%
        filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%
        distinct(SYMBOL,CHROM, POS,ID, REF, ALT, group, .keep_all = TRUE)
    
    # Get counts
    total_case <- sum(expanded_df$group == "Case")
    total_control <- sum(expanded_df$group == "Control")
    total_variants <- nrow(expanded_df)
    
    # Build correct contingency table
    contingency_table <- matrix(
        c(total_case, total_control,
        total_variants - total_case, total_variants - total_control),
        nrow = 2,
        dimnames = list(
        Group = c("In Group", "Not In Group"),
        Status = c("Case", "Control")
        )
    )
    
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
        NON_patho = perform_tests_for_all_rows(result_df_NON_patho),
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



# visualization 
manhattan_plot_V <- ggplot(final_results, aes(x = Row, y = -log10(p_adj_global))) +
  geom_point(aes(color = Test, size = OR), alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  facet_grid(Test ~ ., scales = "free_y") +
  scale_color_viridis_d(option = "plasma") +
  labs(title = "fisher_chiSquare_Across Tests",
       x = "Gene Set Category",
       y = "-log10(FDR-adjusted p-value)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.y = element_text(angle = 0))

manhattan_plot_CC <- ggplot(fisher_results_df, aes(x = Test, y = -log10(p_adj))) +
    geom_point(aes(color = Test, size = OR), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    scale_color_viridis_d(option = "plasma") +
    labs(title = "fisher_chiSquare_Across Tests",
         x = "Gene Set Category",
         y = "-log10(FDR-adjusted p-value)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(angle = 0))


manhattan_plot_ind <- ggplot(fisher_results_ind_df, aes(x = Row, y = -log10(p_adj_global))) +
    geom_point(aes(color = Test, size = OR), alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    facet_grid(Test ~ ., scales = "free_y") +
    scale_color_viridis_d(option = "plasma") +
    labs(title = "fisher_chiSquare_Across Tests",
         x = "Gene Set Category",
         y = "-log10(FDR-adjusted p-value)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(angle = 0))


volcano_plot_ind <- fisher_results_ind_df %>%
  mutate(log_OR = log2(OR),
         log_p = -log10(p_adj_global)) %>%
  ggplot(aes(x = log_OR, y = log_p)) +
  geom_point(aes(color = Test, shape = p_adj_global < 0.05), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), color = "red") +
  ggrepel::geom_text_repel(
    data = . %>% filter(p_adj_global < 0.05),
    aes(label = Row), max.overlaps = 20
  ) +
  scale_shape_manual(values = c(1, 16)) +
  labs(title = "Variant Carrier Association",
       x = "log2(Odds Ratio)",
       y = "-log10(FDR p-value)")

volcano_plot_CC <- fisher_results_df %>%
    mutate(log_OR = log2(OR),
           log_p = -log10(p_adj)) %>%
    ggplot(aes(x = log_OR, y = log_p)) +
    geom_point(aes(color = Test, shape = p_adj < 0.05), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "red") +
    ggrepel::geom_text_repel(
        data = . %>% filter(p_adj < 0.05),
        aes(label = Test), max.overlaps = 20
    ) +
    scale_shape_manual(values = c(1, 16)) +
    labs(title = "Variant Carrier Association",
         x = "log2(Odds Ratio)",
         y = "-log10(FDR p-value)")

volcano_plot_V <- final_results %>%
    mutate(log_OR = log2(OR),
           log_p = -log10(p_adj_global)) %>%
    ggplot(aes(x = log_OR, y = log_p)) +
    geom_point(aes(color = Test, shape = p_adj_global < 0.05), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "red") +
    ggrepel::geom_text_repel(
        data = . %>% filter(p_adj_global < 0.05),
        aes(label = Row), max.overlaps = 20
    ) +
    scale_shape_manual(values = c(1, 16)) +
    labs(title = "Variant Carrier Association",
         x = "log2(Odds Ratio)",
         y = "-log10(FDR p-value)")


dot_plot_V <- final_results %>%
  ggplot(aes(x = OR, y = Row)) +
  geom_point(aes(size = Cases_With_Variant, 
                 color = -log10(p_adj_global))) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_color_viridis_c(option = "magma") +
  facet_grid(Test ~ ., scales = "free", space = "free") +
  labs(title = "Variant Enrichment Patterns",
       x = "Odds Ratio (Case vs Control)",
       y = "Gene Set Category") +
  theme_minimal()       


  dot_plot_ind <- fisher_results_ind_df %>%
    ggplot(aes(x = OR, y = Row)) +
    geom_point(aes(size = Case_Carriers, 
                   color = -log10(p_adj_global))) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_color_viridis_c(option = "magma") +
    facet_grid(Test ~ ., scales = "free", space = "free") +
    labs(title = "Variant Enrichment Patterns",
         x = "Odds Ratio (Case vs Control)",
         y = "Gene Set Category") +
    theme_minimal()      


setwd("/home/rachele/SNVs/results/stats/euro")


#save plots
ggsave("Vulcano_CaseControl_association.png", plot = volcano_plot_CC, width = 10 , height = 8 )
ggsave("Vulcano_Individuals_association_ind.png", plot = volcano_plot_ind, width = 10 , heigh =8 )
ggsave("Vulcano_Variants_association.png",plot=volcano_plot_V,width=10,heigh=8)
ggsave("Manhattan_CaseControl.png", plot= manhattan_plot_CC,width=10 , heigh =8  )
ggsave("Manhattan_Individuals_ind.png", plot= manhattan_plot_ind,width=10 , heigh =8 )
ggsave("Manhattan_Variants.png",plot= manhattan_plot_V,width =10 , heigh =8 )
ggsave("Dot_Individuals_association_ind.png", plot = dot_plot_ind, width = 10 , heigh =8 )
ggsave("Dot_Variants_association.png", plot = dot_plot_V, width = 10 , heigh =8 )

#save tabels 
write.csv(final_results, "Fisher_Variants.csv")
write.csv(fisher_results_ind_df, "Fisher_Individuals.csv")
write.csv(fisher_results_df, "Fisher_CaseControl.csv")
write.csv(rate_confront_df,"rate_confront.csv")




