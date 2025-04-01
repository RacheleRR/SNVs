#GET GENE NAMES 



 library(dplyr)
library(readr)  # For read_tsv
library(tidyr)
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


#get labels 

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


pred_NON_patho_MPC_pLI_VEP_filtered$sample_label_1 <- sapply(pred_NON_patho_MPC_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$sample_label_1 <- sapply(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$sample_label_1 <- sapply(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
pred_PTV_pLI_VEP_filtered$sample_label_1 <- sapply(pred_PTV_pLI_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)
pred_PTV_VEP_filtered$sample_label_1 <- sapply(pred_PTV_VEP_filtered$SAMPLES,check_outlier_label,manifest_correct)


# Function to split dataframe into case and control
split_case_control <- function(df) {
case_df <- df[df$sample_label_1 == "Case", ]
control_df <- df[df$sample_label_1 == "Control", ]
return(list(case = case_df, control = control_df))
}
# Applying the function to each dataset
split_NON_patho <- split_case_control(pred_NON_patho_MPC_pLI_VEP_filtered)
split_patho_MPC_ALPHAMISSENSE_pLI <- split_case_control(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered)
split_patho_MPC_ALPHAMISSENSE <- split_case_control(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered)
split_PTV_pLI <- split_case_control(pred_PTV_pLI_VEP_filtered)
split_PTV <- split_case_control(pred_PTV_VEP_filtered)
# Accessing the case and control dataframes
case_NON_patho <- split_NON_patho$case
control_NON_patho <- split_NON_patho$control
case_patho_MPC_ALPHAMISSENSE_pLI <- split_patho_MPC_ALPHAMISSENSE_pLI$case
control_patho_MPC_ALPHAMISSENSE_pLI <- split_patho_MPC_ALPHAMISSENSE_pLI$control
case_patho_MPC_ALPHAMISSENSE <- split_patho_MPC_ALPHAMISSENSE$case
control_patho_MPC_ALPHAMISSENSE <- split_patho_MPC_ALPHAMISSENSE$control
case_PTV_pLI <- split_PTV_pLI$case
control_PTV_pLI <- split_PTV_pLI$control
case_PTV <- split_PTV$case
control_PTV <- split_PTV$control
# Function to filter unique genes and save as TSV
save_unique_genes <- function(df, filename) {
unique_df <- df[!duplicated(df$SYMBOL), ]  # Remove duplicates based on GENE column
write.table(unique_df, file = filename, sep = "\t", row.names = FALSE, quote = FALSE)
}
# Apply to case and control datasets
save_unique_genes(case_NON_patho, "case_NON_patho.tsv")
save_unique_genes(control_NON_patho, "control_NON_patho.tsv")
save_unique_genes(case_patho_MPC_ALPHAMISSENSE_pLI, "case_patho_MPC_ALPHAMISSENSE_pLI.tsv")
save_unique_genes(control_patho_MPC_ALPHAMISSENSE_pLI, "control_patho_MPC_ALPHAMISSENSE_pLI.tsv")
save_unique_genes(case_patho_MPC_ALPHAMISSENSE, "case_patho_MPC_ALPHAMISSENSE.tsv")
save_unique_genes(control_patho_MPC_ALPHAMISSENSE, "control_patho_MPC_ALPHAMISSENSE.tsv")
save_unique_genes(case_PTV_pLI, "case_PTV_pLI.tsv")
save_unique_genes(control_PTV_pLI, "control_PTV_pLI.tsv")
save_unique_genes(case_PTV, "case_PTV.tsv")
save_unique_genes(control_PTV, "control_PTV.tsv")