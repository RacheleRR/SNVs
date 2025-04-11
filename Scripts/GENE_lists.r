#GET GENE NAMES ASSOCIATED TO SCHIZO ETC 


# Load necessary libraries for data manipulation and analysis
library(readr)
library(data.table)
library(dplyr)
library(gprofiler2)
library(httr)
library(readxl)
library(tidyr)
library(stringr)
library(openxlsx)
# Load my data from CSV files
  SCHEMA <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SCHEMA.csv")
  GWAS_120 <- read.delim("~/GWAS_120.csv")       
  BipEx_Bipolar <- read_csv("/home/rachele/Documents/old/geneset_confrontation/Bipolar_Disorder.csv")
  SFARI <- read_csv("/home/rachele/Documents/old/geneset_confrontation/SFARI.csv")

# Clean GENESETS
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


# Apply 


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

    pred_NON_patho_MPC_pLI_VEP_filtered <-pred_NON_patho_MPC_pLI_VEP_filtered %>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    pred_NON_patho_MPC_VEP_filtered <-pred_NON_patho_MPC_VEP_filtered%>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered <-pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered %>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    pred_patho_MPC_ALPHAMISSENSE_VEP_filtered <-pred_patho_MPC_ALPHAMISSENSE_VEP_filtered%>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    pred_PTV_pLI_VEP_filtered <-pred_PTV_pLI_VEP_filtered%>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    pred_PTV_VEP_filtered <-pred_PTV_VEP_filtered%>% mutate(count = sapply(strsplit(sample_label_2, ","), length))
    







# Function to split dataframe into case and control
split_case_control <- function(df) {
case_df <- df[df$sample_label_3 == "Case", ]
control_df <- df[df$sample_label_3 == "Control", ]
case_df <- case_df %>% distinct(SYMBOL, CHROM,REF,ALT, POS, ID, .keep_all = TRUE) 
control_df <- control_df %>% distinct(SYMBOL, CHROM,REF,ALT, POS,ID,.keep_all=TRUE)
return(list(case = case_df, control = control_df))
}


# Applying the function to each dataset
split_NON_patho <- split_case_control(pred_NON_patho_MPC_VEP_filtered)
split_NON_patho_pLI <- split_case_control(pred_NON_patho_MPC_pLI_VEP_filtered)
split_patho_MPC_ALPHAMISSENSE_pLI <- split_case_control(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered)
split_patho_MPC_ALPHAMISSENSE <- split_case_control(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered)
split_PTV_pLI <- split_case_control(pred_PTV_pLI_VEP_filtered)
split_PTV <- split_case_control(pred_PTV_VEP_filtered)
# Accessing the case and control dataframes
case_NON_patho <- split_NON_patho$case
control_NON_patho <- split_NON_patho$control
case_NON_patho_pLI <- split_NON_patho_pLI$case
control_NON_patho_pLI <- split_NON_patho_pLI$control
case_patho_MPC_ALPHAMISSENSE_pLI <- split_patho_MPC_ALPHAMISSENSE_pLI$case
control_patho_MPC_ALPHAMISSENSE_pLI <- split_patho_MPC_ALPHAMISSENSE_pLI$control
case_patho_MPC_ALPHAMISSENSE <- split_patho_MPC_ALPHAMISSENSE$case
control_patho_MPC_ALPHAMISSENSE <- split_patho_MPC_ALPHAMISSENSE$control
case_PTV_pLI <- split_PTV_pLI$case
control_PTV_pLI <- split_PTV_pLI$control
case_PTV <- split_PTV$case
control_PTV <- split_PTV$control



# set working directory 
 setwd("/home/rachele/SNVs/results/stats/global/non_private/genes")
# Function to filter unique genes and save as TSV
save_unique_genes <- function(df, filename) {
unique_df <- df[!duplicated(df$SYMBOL), ]  # Remove duplicates based on GENE column
write.table(unique_df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}

# Apply to case and control datasets
save_unique_genes(case_NON_patho, "case_NON_patho.tsv")
save_unique_genes(control_NON_patho, "control_NON_patho.tsv")
save_unique_genes(case_NON_patho_pLI, "case_NON_patho_pLI.tsv")
save_unique_genes(control_NON_patho_pLI, "control_NON_patho_pLI.tsv")
save_unique_genes(case_patho_MPC_ALPHAMISSENSE_pLI, "case_patho_MPC_ALPHAMISSENSE_pLI.tsv")
save_unique_genes(control_patho_MPC_ALPHAMISSENSE_pLI, "control_patho_MPC_ALPHAMISSENSE_pLI.tsv")
save_unique_genes(case_patho_MPC_ALPHAMISSENSE, "case_patho_MPC_ALPHAMISSENSE.tsv")
save_unique_genes(control_patho_MPC_ALPHAMISSENSE, "control_patho_MPC_ALPHAMISSENSE.tsv")
save_unique_genes(case_PTV_pLI, "case_PTV_pLI.tsv")
save_unique_genes(control_PTV_pLI, "control_PTV_pLI.tsv")
save_unique_genes(case_PTV, "case_PTV.tsv")
save_unique_genes(control_PTV, "control_PTV.tsv")


# GENES IN SCHIZO ETC 
# Define the unique gene lists
genes_schema_pval <- unique(SCHEMA_pVal$Gene)
genes_bipolar <- unique(c(BipEx_Bipolar_p_val_PTV$Gene, BipEx_Bipolar_p_val_Missense$Gene))
genes_sfari_non_syndromic <- unique(SFARI_non_syndromic_lower_3$Gene)
genes_gwas <- unique(GWAS_120$Gene)
genes_schema_or <- unique(SCHEMA_OR$Gene)
genes_sfari_syndromic <- unique(SFARI_syndromic$Gene)
genes_brain_ntm <- unique(brain_gene_consensus_ntm_consensus_no_pitular$Gene.name)
genes_brain_filt <- unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene.name)

# Function to filter genes that are present in both the case dataset and a given gene list
filter_genes <- function(case_df, gene_list, type) {
  filtered_df <- case_df %>% filter(SYMBOL %in% gene_list) %>% mutate(type = type)
  return(filtered_df)
}

# Filter the case datasets based on the gene lists
filtered_case_patho_MPC_ALPHAMISSENSE_schema_pval <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_schema_pval, "schema_pval")
filtered_case_patho_MPC_ALPHAMISSENSE_bipolar <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_bipolar, "bipolar")
filtered_case_patho_MPC_ALPHAMISSENSE_sfari_non_syndromic <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_patho_MPC_ALPHAMISSENSE_gwas <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_gwas, "gwas")
filtered_case_patho_MPC_ALPHAMISSENSE_schema_or <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_schema_or, "schema_or")
filtered_case_patho_MPC_ALPHAMISSENSE_sfari_syndromic <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_patho_MPC_ALPHAMISSENSE_brain_ntm <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_brain_ntm, "brain_ntm")
filtered_case_patho_MPC_ALPHAMISSENSE_brain_filt <- filter_genes(case_patho_MPC_ALPHAMISSENSE, genes_brain_filt, "brain_filt")

filtered_case_PTV_schema_pval <- filter_genes(case_PTV, genes_schema_pval, "schema_pval")
filtered_case_PTV_bipolar <- filter_genes(case_PTV, genes_bipolar, "bipolar")
filtered_case_PTV_sfari_non_syndromic <- filter_genes(case_PTV, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_PTV_gwas <- filter_genes(case_PTV, genes_gwas, "gwas")
filtered_case_PTV_schema_or <- filter_genes(case_PTV, genes_schema_or, "schema_or")
filtered_case_PTV_sfari_syndromic <- filter_genes(case_PTV, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_PTV_brain_ntm <- filter_genes(case_PTV, genes_brain_ntm, "brain_ntm")
filtered_case_PTV_brain_filt <- filter_genes(case_PTV, genes_brain_filt, "brain_filt")

filtered_case_PTV_pLI_schema_pval <- filter_genes(case_PTV_pLI, genes_schema_pval, "schema_pval")
filtered_case_PTV_pLI_bipolar <- filter_genes(case_PTV_pLI, genes_bipolar, "bipolar")
filtered_case_PTV_pLI_sfari_non_syndromic <- filter_genes(case_PTV_pLI, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_PTV_pLI_gwas <- filter_genes(case_PTV_pLI, genes_gwas, "gwas")
filtered_case_PTV_pLI_schema_or <- filter_genes(case_PTV_pLI, genes_schema_or, "schema_or")
filtered_case_PTV_pLI_sfari_syndromic <- filter_genes(case_PTV_pLI, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_PTV_pLI_brain_ntm <- filter_genes(case_PTV_pLI, genes_brain_ntm, "brain_ntm")
filtered_case_PTV_pLI_brain_filt <- filter_genes(case_PTV_pLI, genes_brain_filt, "brain_filt")

filtered_case_patho_MPC_ALPHAMISSENSE_pLI_schema_pval <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_schema_pval, "schema_pval")
filtered_case_patho_MPC_ALPHAMISSENSE_pLI_bipolar <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_bipolar, "bipolar")
filtered_case_patho_MPC_ALPHAMISSENSE_pLI_sfari_non_syndromic <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_sfari_non_syndromic, "sfari_non_syndromic")
filtered_case_patho_MPC_ALPHAMISSENSE_pLI_gwas <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_gwas, "gwas")
filtered_case_patho_MPC_ALPHAMISSENSE_pLI_schema_or <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_schema_or, "schema_or")
filtered_case_patho_MPC_ALPHAMISSENSE_pLI_sfari_syndromic <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_sfari_syndromic, "sfari_syndromic")
filtered_case_patho_MPC_ALPHAMISSENSE_pLI_brain_ntm <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_brain_ntm, "brain_ntm")
filtered_case_patho_MPC_ALPHAMISSENSE_pLI_brain_filt <- filter_genes(case_patho_MPC_ALPHAMISSENSE_pLI, genes_brain_filt, "brain_filt")

# Combine the filtered dataframes into a single dataframe for each case dataset
combined_case_patho_MPC_ALPHAMISSENSE <- bind_rows(
  filtered_case_patho_MPC_ALPHAMISSENSE_schema_pval,
  filtered_case_patho_MPC_ALPHAMISSENSE_bipolar,
  filtered_case_patho_MPC_ALPHAMISSENSE_sfari_non_syndromic,
  filtered_case_patho_MPC_ALPHAMISSENSE_gwas,
  filtered_case_patho_MPC_ALPHAMISSENSE_schema_or,
  filtered_case_patho_MPC_ALPHAMISSENSE_sfari_syndromic,
  filtered_case_patho_MPC_ALPHAMISSENSE_brain_ntm,
  filtered_case_patho_MPC_ALPHAMISSENSE_brain_filt
)

combined_case_PTV <- bind_rows(
  filtered_case_PTV_schema_pval,
  filtered_case_PTV_bipolar,
  filtered_case_PTV_sfari_non_syndromic,
  filtered_case_PTV_gwas,
  filtered_case_PTV_schema_or,
  filtered_case_PTV_sfari_syndromic,
  filtered_case_PTV_brain_ntm,
  filtered_case_PTV_brain_filt
)

combined_case_PTV_pLI <- bind_rows(
  filtered_case_PTV_pLI_schema_pval,
  filtered_case_PTV_pLI_bipolar,
  filtered_case_PTV_pLI_sfari_non_syndromic,
  filtered_case_PTV_pLI_gwas,
  filtered_case_PTV_pLI_schema_or,
  filtered_case_PTV_pLI_sfari_syndromic,
  filtered_case_PTV_pLI_brain_ntm,
  filtered_case_PTV_pLI_brain_filt
)

combined_case_patho_MPC_ALPHAMISSENSE_pLI <- bind_rows(
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_schema_pval,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_bipolar,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_sfari_non_syndromic,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_gwas,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_schema_or,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_sfari_syndromic,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_brain_ntm,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_brain_filt
)

# Save the combined dataframes to new TSV files
write.table(combined_case_patho_MPC_ALPHAMISSENSE, "combined_case_patho_MPC_ALPHAMISSENSE.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(combined_case_PTV, "combined_case_PTV.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(combined_case_PTV_pLI, "combined_case_PTV_pLI.tsv", sep="\t", row.names=FALSE, quote=FALSE)
write.table(combined_case_patho_MPC_ALPHAMISSENSE_pLI, "combined_case_patho_MPC_ALPHAMISSENSE_pLI.tsv", sep="\t", row.names=FALSE, quote=FALSE)


# FIND GENES WITHIN MULTIPLES
# Define the unique gene lists
genes_schema_pval <- unique(SCHEMA_pVal$Gene)
genes_gwas <- unique(GWAS_120$Gene)
genes_brain_filt <- unique(brain_gene_consensus_filtered_consensus_no_pitular$Gene)
genes_case_PTV <- unique(case_PTV$SYMBOL)

# Find the genes that are in case_PTV and brain_filt
genes_case_PTV_brain <- intersect(genes_case_PTV, genes_brain_filt)

# Find the genes that are either in schema_pval or gwas
genes_schema_or_gwas <- union(genes_schema_pval, genes_gwas)

# Find the genes that are in case_PTV_brain and either in schema_pval or gwas
common_genes <- intersect(genes_case_PTV_brain, genes_schema_or_gwas)

# Create a dataframe with the unique genes
common_genes_df <- data.frame(Gene = common_genes)

# Save the dataframe to a TSV file
write.table(common_genes_df, "common_genes_case_PTV_brain_schema_or_gwas.tsv", sep="\t", row.names=FALSE, quote=FALSE)



# !Excel document 
# case_control genes
  #Create a list of all case and control datasets
  case_control_list <- list(
    case_NON_patho = case_NON_patho,
    control_NON_patho = control_NON_patho,
    case_NON_patho_pLI = case_NON_patho_pLI,
    control_NON_patho_pLI = control_NON_patho_pLI,
    case_patho_MPC_ALPHAM_pLI = case_patho_MPC_ALPHAMISSENSE_pLI,
    control_patho_MPC_ALPHAM_pLI = control_patho_MPC_ALPHAMISSENSE_pLI,
    case_patho_MPC_ALPHAM = case_patho_MPC_ALPHAMISSENSE,
    control_patho_MPC_ALPHAM = control_patho_MPC_ALPHAMISSENSE,
    case_PTV_pLI = case_PTV_pLI,
    control_PTV_pLI = control_PTV_pLI,
    case_PTV = case_PTV,
    control_PTV = control_PTV
  )

  #Save the list to an Excel file
  write.xlsx(case_control_list, file = "case_control_datasets.xlsx")

# case control genes unqiue 
  #Create a function to filter for unique genes based on the SYMBOL column
  filter_unique_genes <- function(df) {
    df %>% distinct(SYMBOL, .keep_all = TRUE)
  }

  #Apply the function to each dataset in the list
  unique_case_control_list <- lapply(case_control_list, filter_unique_genes)

  #Save the unique datasets to another Excel file
  write.xlsx(unique_case_control_list, file = "unique_case_control_datasets.xlsx")

#  case control with disaes and brain 
#Create a list of all combined and filtered datasets
#Combine the filtered dataframes into a single dataframe for each case dataset
combined_case_patho_MPC_ALPHAMISSENSE <- bind_rows(
  filtered_case_patho_MPC_ALPHAMISSENSE_schema_pval,
  filtered_case_patho_MPC_ALPHAMISSENSE_bipolar,
  filtered_case_patho_MPC_ALPHAMISSENSE_sfari_non_syndromic,
  filtered_case_patho_MPC_ALPHAMISSENSE_gwas,
  filtered_case_patho_MPC_ALPHAMISSENSE_schema_or,
  filtered_case_patho_MPC_ALPHAMISSENSE_sfari_syndromic
)

combined_case_PTV <- bind_rows(
  filtered_case_PTV_schema_pval,
  filtered_case_PTV_bipolar,
  filtered_case_PTV_sfari_non_syndromic,
  filtered_case_PTV_gwas,
  filtered_case_PTV_schema_or,
  filtered_case_PTV_sfari_syndromic
)

combined_case_PTV_pLI <- bind_rows(
  filtered_case_PTV_pLI_schema_pval,
  filtered_case_PTV_pLI_bipolar,
  filtered_case_PTV_pLI_sfari_non_syndromic,
  filtered_case_PTV_pLI_gwas,
  filtered_case_PTV_pLI_schema_or,
  filtered_case_PTV_pLI_sfari_syndromic
)

combined_case_patho_MPC_ALPHAMISSENSE_pLI <- bind_rows(
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_schema_pval,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_bipolar,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_sfari_non_syndromic,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_gwas,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_schema_or,
  filtered_case_patho_MPC_ALPHAMISSENSE_pLI_sfari_syndromic
)

combined_filtered_list <- list(
  combined_case_patho_MPC_ALPHAM = combined_case_patho_MPC_ALPHAMISSENSE,
  filt_case_path_MPC_AM_brainntm = filtered_case_patho_MPC_ALPHAMISSENSE_brain_ntm,
  filt_case_path_MPC_AM_brainfil = filtered_case_patho_MPC_ALPHAMISSENSE_brain_filt,
  combined_case_PTV = combined_case_PTV,
  combined_case_PTV_pLI = combined_case_PTV_pLI,
  combined_case_patho_MPC_AM_pLI = combined_case_patho_MPC_ALPHAMISSENSE_pLI,
  filt_case_PTV_brain_ntm = filtered_case_PTV_brain_ntm,
  filt_case_PTV_brain_filt = filtered_case_PTV_brain_filt,
  filt_case_PTV_pLI_brain_ntm = filtered_case_PTV_pLI_brain_ntm,
  filt_case_PTV_pLI_brain_filt = filtered_case_PTV_pLI_brain_filt,
  filt_case_patho_pLI_brainntm = filtered_case_patho_MPC_ALPHAMISSENSE_pLI_brain_ntm,
  filt_case_patho_pLI_brainfil = filtered_case_patho_MPC_ALPHAMISSENSE_pLI_brain_filt
)

#Save the list to an Excel file
write.xlsx(combined_filtered_list, file = "combined_filtered_datasets.xlsx")

# case control with disease and brain unique 
#Apply the unique filter function to each dataset in the combined_filtered_list
unique_combined_filtered_list <- lapply(combined_filtered_list, filter_unique_genes)

#Save the unique datasets to another Excel file
write.xlsx(unique_combined_filtered_list, file = "unique_combined_filtered_datasets.xlsx")

