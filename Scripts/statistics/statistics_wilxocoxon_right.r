# statistics wilxocoxon right 


library(readr)
library(data.table)
# library(tidyverse)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)



#  Load my datasets 
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

# LOAD MINE     
# LOAD PERSONAL DATA 
    # Specify the folder containing your TSV files
    folder_path <- "/home/rachele/SNVs/results/global"

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


    UHR_NA_SAMPLE_IDS <- read.csv("~/UHR_NA_SAMPLE_IDS.csv", sep="")
    #example S36827

colnames(manifest_correct) <- c("Status", "sample_id", "group")


# stuff
process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE) {

    if(private){
        df <- df %>% filter(count == 1)  # Only filter if private=TRUE
    }
    
    # Assume expanded_df is already created as per your code
   expanded_df <- df %>%
    separate_rows(SAMPLES, sep = ",") %>%
    mutate(sample_id = trimws(SAMPLES)) %>%
    left_join(
        manifest_correct %>% distinct(sample_id, group), 
        by = "sample_id"
    ) %>%
    filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%
    distinct(SYMBOL, CHROM, REF, ALT, POS, sample_id, .keep_all = TRUE)


    expanded_df_pure <- df %>% 
        separate_rows(SAMPLES, sep = ",") %>%  # Split SAMPLES only
        mutate(
            sample_id = trimws(SAMPLES),
            group = sample_label_3  # Use pre-computed labels
        ) %>% filter(group %in% c("Case", "Control")) %>% 
        filter(!sample_id %in% UHR_NA_SAMPLE_IDS$V1) %>%  # Remove UHR_NA
        distinct(SYMBOL, CHROM,REF,ALT, POS, sample_id, .keep_all = TRUE)  # Unique variant/sample
        
   
    # Filter based on the genes in the brain, brain filter_ntpm, and additional gene sets
    expanded_df <- expanded_df %>%
        mutate(brain = ifelse(SYMBOL %in% brain_gene_consensus_filtered_consensus_no_pitular$Gene.name, 1, 0),
               brain_filter_ntpm = ifelse(SYMBOL %in% brain_gene_consensus_ntm_consensus_no_pitular$Gene.name, 1, 0),
               schema_pval = ifelse(SYMBOL %in% genes_schema_pval, 1, 0),
               bipolar = ifelse(SYMBOL %in% genes_bipolar, 1, 0),
               sfari_non_syndromic = ifelse(SYMBOL %in% genes_sfari_non_syndromic, 1, 0),
               schema_or = ifelse(SYMBOL %in% genes_schema_or, 1, 0))  

    expanded_df_pure <- expanded_df_pure %>%
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

    result_df_pure <- expanded_df_pure %>%  
    group_by(sample_id, group) %>%
        summarise(
            Number_of_Variants_pure = n(),
            Number_of_Brain_Variants_pure = sum(brain),
            Number_of_Brain_Filter_NTPM_Variants_pure = sum(brain_filter_ntpm),
            Number_of_Schema_Pval_Variants_pure = sum(schema_pval),
            Number_of_Bipolar_Variants_pure = sum(bipolar),
            Number_of_Sfari_Non_Syndromic_Variants_pure = sum(sfari_non_syndromic),
            Number_of_Schema_Or_Variants_pure = sum(schema_or),
            .groups = "drop"
        )   
    merged_results <- result_df %>%
    full_join( result_df_pure,
    by = "sample_id"  ) %>% mutate(
    group = coalesce(group.x, group.y)  # Combine group.x and group.y into "group"
    ) %>%select(-group.x, -group.y) %>%       # Remove the original group columns
    distinct()                            # Remove duplicates

    return(merged_results)
}

# Create separate data frames for each call to process_df
    result_df_NON_patho_pLI <- process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_NON_patho <- process_df(pred_NON_patho_MPC_VEP_filtered, "NON_patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE )
    result_df_patho_pLI <- process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_patho <- process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_PTV_pLI <- process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)
    result_df_PTV <- process_df(pred_PTV_VEP_filtered, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=FALSE)


    result_df_NON_patho_pLI <- process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_NON_patho <- process_df(pred_NON_patho_MPC_VEP_filtered, "NON_patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE )
    result_df_patho_pLI <- process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_patho <- process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_PTV_pLI <- process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)
    result_df_PTV <- process_df(pred_PTV_VEP_filtered, "PTV", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or,private=TRUE)


manifest_correct <- manifest_correct %>% filter(!grepl("UHR_NA", group))


# ADD ALL SAMPLES INDEPENDENTLY OF PRESENCE/ABSENCE OF VARIANTS 
# Ensure consistency in column names and perform left join for all dataframes
complete_df_NON_patho_pLI <- manifest_correct %>%
    left_join(result_df_NON_patho_pLI, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_df_NON_patho <- manifest_correct %>%
    left_join(result_df_NON_patho, by = c("sample_id", "group")) %>%
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




# statistical test 
# without normalize

    analyze_metrics <- function(df) {
    # 1. Test each variant count metric
    columns_to_test <- grep("^Number_of", names(df), value = TRUE)
    
    results <- map_dfr(columns_to_test, function(col) {
        # Initialize result row
        result_row <- tibble(
        metric = col,
        case_zero = FALSE,
        control_zero = FALSE,
        shapiro_case_p = NA_real_,
        shapiro_control_p = NA_real_,
        primary_test = NA_character_,
        primary_p = NA_real_,
        secondary_test = NA_character_,
        secondary_p = NA_real_,
        note = NA_character_
        )
        
        # Split data
        case_data <- df[[col]][df$group == "Case"]
        control_data <- df[[col]][df$group == "Control"]
        
        # Check for all-zero cases
        result_row$case_zero <- all(case_data == 0)
        result_row$control_zero <- all(control_data == 0)
        
        if (result_row$case_zero || result_row$control_zero) {
        result_row$note <- case_when(
            result_row$case_zero && result_row$control_zero ~ "All zeros in both groups",
            result_row$case_zero ~ "All zeros in Case group",
            TRUE ~ "All zeros in Control group"
        )
        return(result_row)
        }
        
        # Handle constant values for Shapiro
        safe_shapiro <- function(x) {
        if (length(unique(x)) < 3) return(list(p.value = NA))
        tryCatch(shapiro.test(x), error = function(e) list(p.value = NA))
        }
        
        # Normality checks
        shapiro_case <- safe_shapiro(case_data)
        shapiro_control <- safe_shapiro(control_data)
        result_row$shapiro_case_p <- shapiro_case$p.value
        result_row$shapiro_control_p <- shapiro_control$p.value
        
        # Determine test strategy
        if (!any(is.na(c(shapiro_case$p.value, shapiro_control$p.value))) &&
        shapiro_case$p.value > 0.05 && 
        shapiro_control$p.value > 0.05) {
        primary <- "t-test"
        secondary <- "wilcox"
        } else {
        primary <- "wilcox"
        secondary <- "t-test"
        }
        
        # Perform tests
        tryCatch({
        t_res <- t.test(df[[col]] ~ df$group) %>% tidy()
        w_res <- wilcox.test(df[[col]] ~ df$group) %>% tidy()
        
        result_row$primary_test <- primary
        result_row$primary_p <- if(primary == "t-test") t_res$p.value else w_res$p.value
        result_row$secondary_test <- secondary
        result_row$secondary_p <- if(secondary == "t-test") t_res$p.value else w_res$p.value
        }, error = function(e) {
        result_row$note <- paste("Test error:", e$message)
        })



        result_row
        })
           # Add FDR-adjusted p-values
         results %>%
         mutate(
            primary_p_adj = p.adjust(primary_p, method = "fdr"),
            secondary_p_adj = p.adjust(secondary_p, method = "fdr")
           ) %>%
          select(
          metric, 
           primary_test, primary_p, primary_p_adj,
             secondary_test, secondary_p, secondary_p_adj,
             everything()
    )
    }
        
    # apply 
    results_pred_NON_patho_pLI <- analyze_metrics(complete_df_NON_patho_pLI) %>%  select(-case_zero, -control_zero)
    results_pred_NON_patho <- analyze_metrics(complete_df_NON_patho) %>%  select(-case_zero, -control_zero)
    results_pred_patho_pLI <- analyze_metrics(complete_df_patho_pLI) %>%  select(-case_zero, -control_zero)
    results_pred_patho <- analyze_metrics(complete_df_patho) %>%  select(-case_zero, -control_zero)
    results_pred_PTV_pLI <- analyze_metrics(complete_df_PTV_pLI) %>%  select(-case_zero, -control_zero)
    results_pred_PTV <- analyze_metrics(complete_df_PTV) %>%  select(-case_zero, -control_zero)






# VISUALIZATION
setwd("/home/rachele/SNVs/results/stats/global/private")
# Save all the plots indiscriminate from significance  
        plot <- function(df, results_df, title_suffix) {
                # 2. Prepare annotation data from test results
                annotation_data <- results_df %>%
                    filter(str_detect(metric, "^Number_of")) %>%
                    mutate(
                        metric_clean = gsub("Number_of_", "", metric),
                        label = paste0(primary_test, "\nFDR p = ", format.pval(primary_p_adj, digits = 2)),
                        group = NA  # Add dummy group variable
                    )
                
                # 3. Reshape for plotting
                plot_data <- df %>%
                    select(sample_id, group, starts_with("Number_of")) %>%
                    pivot_longer(
                        cols = -c(sample_id, group),
                        names_to = "metric",
                        values_to = "value"
                    ) %>%
                    mutate(
                        metric_clean = gsub("Number_of_", "", metric)
                    )
                
                # 4. Create the plot
                ggplot(plot_data, aes(x = group, y = value, fill = group)) +
                    geom_boxplot(outlier.shape = NA) +
                    geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
                    facet_wrap(~ metric_clean, scales = "free_y", ncol = 3) +
                    geom_text(
                        data = annotation_data,
                        aes(x = 1.5, y = Inf, label = label),
                        inherit.aes = FALSE,  # Don't inherit aesthetics from main plot
                        vjust = 1.5, size = 3, color = "black"
                    ) +
                    labs(
                        title = paste("Case vs. Control:", title_suffix),
                        x = "Group",
                        y = "Variants per individual (raw)",
                        fill = "Group"
                    ) +
                    theme_bw() +
                    theme(
                        strip.background = element_blank(),
                        strip.text = element_text(face = "bold")
                    )
            }

            plot_NON_patho_pLI <- plot(
            complete_df_NON_patho_pLI, 
            results_pred_NON_patho_pLI,
            "Non-pathogenic (pLI filtered)"
            )

            plot_NON_patho <- plot(
            complete_df_NON_patho, 
            results_pred_NON_patho,
            "Non-pathogenic"
            )

            plot_patho_pLI <- plot(
            complete_df_patho_pLI, 
            results_pred_patho_pLI,
            "Pathogenic (pLI filtered)"
            )

            plot_patho <- plot(
            complete_df_patho, 
            results_pred_patho,
            "Pathogenic"
            )

            plot_PTV_pLI <- plot(
            complete_df_PTV_pLI, 
            results_pred_PTV_pLI,
            "PTV (pLI filtered)"
            )

            plot_PTV <- plot(
            complete_df_PTV, 
            results_pred_PTV,
            "PTV"
            )

            # Display one plot as example
            plot_NON_patho

            # To save all plots:
            plots <- list(
            plot_NON_patho_pLI, plot_NON_patho,
            plot_patho_pLI, plot_patho,
            plot_PTV_pLI, plot_PTV
            )


            walk2(plots, c("NON_patho_pLI", "NON_patho", "patho_pLI", "patho", "PTV_pLI", "PTV"), 
                ~ ggsave(paste0("variant_comparison_", .y, ".png"), .x, width = 10, height = 8))


        #save tabels 
        combined_results_non_norm <- bind_rows(
        "NON_patho_pLI" = results_pred_NON_patho_pLI,
        "NON_patho" = results_pred_NON_patho,
        "patho_pLI" = results_pred_patho_pLI,
        "patho" = results_pred_patho,
        "PTV_pLI" = results_pred_PTV_pLI,
        "PTV" = results_pred_PTV,
        .id = "analysis_type"  # Adds a column to identify the source
        )

# ONLY  show and save significant results
    plot_significant <- function(df, results_df, title_suffix) {
    # Filter for significant metrics (Wilcoxon FDR p < 0.05)
    sig_metrics <- results_df %>%
        filter(
        str_detect(metric, "^Number_of"),
        primary_p_adj < 0.05, 
        primary_test == "wilcox"  # Focus on Wilcoxon results
        ) %>%
        pull(metric)
    
    # Return NULL if no significant results
    if (length(sig_metrics) == 0) {
        message("No significant results for: ", title_suffix)
        return(NULL)
    }
    
    # Prepare annotation data (only for significant metrics)
    annotation_data <- results_df %>%
        filter(metric %in% sig_metrics) %>%
        mutate(
        metric_clean = gsub("Number_of_", "", metric),
        label = paste0("Wilcoxon\nFDR p = ", format.pval(primary_p_adj, digits = 2)),
        group = NA  # Dummy group for positioning
        )
    
    # Reshape data for plotting
    plot_data <- df %>%
        select(sample_id, group, all_of(sig_metrics)) %>%
        pivot_longer(
        cols = -c(sample_id, group),
        names_to = "metric",
        values_to = "value"
        ) %>%
        mutate(
        metric_clean = gsub("Number_of_", "", metric)
        )
    
    # Create the plot
    ggplot(plot_data, aes(x = group, y = value, fill = group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
        facet_wrap(~ metric_clean, scales = "free_y") +
        geom_text(
        data = annotation_data,
        aes(x = 1.5, y = Inf, label = label),
        vjust = 1.5, size = 3, color = "black"
        ) +
        labs(
        title = paste("Significant results:", title_suffix),
        x = "Group",
        y = "Variants per individual",
        fill = "Group"
        ) +
        theme_bw() +
        theme(
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")
        )
    }

    # Generate and save only significant plots
    sig_plots <- list(
    plot_significant(complete_df_NON_patho_pLI, results_pred_NON_patho_pLI, "Non-pathogenic (pLI)"),
    plot_significant(complete_df_NON_patho, results_pred_NON_patho, "Non-pathogenic"),
    plot_significant(complete_df_patho_pLI, results_pred_patho_pLI, "Pathogenic (pLI)"),
    plot_significant(complete_df_patho, results_pred_patho, "Pathogenic"),
    plot_significant(complete_df_PTV_pLI, results_pred_PTV_pLI, "PTV (pLI)"),
    plot_significant(complete_df_PTV, results_pred_PTV, "PTV")
    ) %>% 
    compact()  # Remove NULL entries (non-significant results)

    # Save significant plots
    walk2(
    sig_plots, 
    c("NON_patho_pLI", "NON_patho", "patho_pLI", "patho", "PTV_pLI", "PTV")[1:length(sig_plots)],
    ~ ggsave(
        paste0("SIGNIFICANT_variant_comparison_", .y, ".png"), 
        .x, 
        width = 10, 
        height = 6
    )
    )


# Save the combined results to a CSV file
write.csv(combined_results_non_norm, "combined_results_non_norm.csv", row.names = FALSE)
















































#posssibilty to do more 


# Load the combined results (normalized or non-normalized)
combined_results <- combined_results_norm  # or combined_results_non_norm

# Filter for the metrics we care about (Number_of_Variants and Number_of_Variants_pure)
p_value_data <- combined_results %>%
    filter(metric %in% c("Number_of_Variants", "Number_of_Variants_pure")) %>%
    select(analysis_type, metric, primary_p_adj) %>%
    mutate(
        metric_clean = case_when(
            metric == "Number_of_Variants" ~ "All Samples",
            metric == "Number_of_Variants_pure" ~ "Pure Samples"
        )
    )


    # List of all complete dataframes
complete_dfs <- list(
  "NON_patho_pLI" = complete_df_NON_patho_pLI,
  "NON_patho" = complete_df_NON_patho,
  "patho_pLI" = complete_df_patho_pLI,
  "patho" = complete_df_patho,
  "PTV_pLI" = complete_df_PTV_pLI,
  "PTV" = complete_df_PTV
)

# Extract and combine the relevant columns
combined_variants <- map_dfr(complete_dfs, ~ {
  .x %>%
    select(sample_id, group, Number_of_Variants, Number_of_Variants_pure) %>%
    pivot_longer(
      cols = c(Number_of_Variants, Number_of_Variants_pure),
      names_to = "variant_type",
      values_to = "count"
    ) %>%
    mutate(
      variant_type = case_when(
        variant_type == "Number_of_Variants" ~ "All Samples",
        variant_type == "Number_of_Variants_pure" ~ "Pure Samples"
      )
    )
}, .id = "analysis_type")




# Merge with the plotting data
plot_data <- combined_variants %>%
    left_join(p_value_data, by = c("analysis_type", "variant_type" = "metric_clean"))


# 3. Create the plot
variant_comparison_plot <- ggplot(plot_data, aes(x = group, y = count, fill = group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 1) +
  facet_grid(variant_type ~ analysis_type, scales = "free_y") +
  
  # Add p-value annotations (using p_value_data)
  geom_text(
    data = p_value_data,
    aes(x = 1.5, y = Inf, label = paste("FDR p =", format.pval(primary_p_adj, digits = 2))),
    inherit.aes = FALSE,  # Ignore aesthetics from main plot
    vjust = 1.5, size = 3, color = "black"
  ) +
  
  labs(
    title = "Case vs. Control: Variant Counts (All vs. Pure Samples)",
    x = "Group",
    y = "Number of Variants (normalized per individual)",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Display the plot
print(variant_comparison_plot)

