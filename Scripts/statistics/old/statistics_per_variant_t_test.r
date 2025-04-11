
Per variant 
 library(dplyr)
library(readr)  # For read_tsv

# Load necessary libraries
library(dplyr)
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


process_df <- function(df, name) {
    # Assume expanded_df is already created as per your code
    expanded_df <- df %>%
        separate_rows(SAMPLES, sample_label_2, sep = ",") %>%
        mutate(
            sample_id = trimws(SAMPLES),
            group = trimws(sample_label_2)
        ) %>%
        distinct()

    
    return(expanded_df)
}

# Create separate data frames for each call to process_df
result_df_NON_patho_pLI <- process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI")
result_df_patho_pLI <- process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI")
result_df_patho <- process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho")
result_df_PTV_pLI <- process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI")
result_df_PTV <- process_df(pred_PTV_VEP_filtered, "PTV")


#perform statistics 
# Step 1: Define the function for statistical analysis
perform_analysis <- function(df) {
  # Create a copy of the original dataframe to avoid overwriting it
  result_df <- df %>%
    group_by(CHROM, POS, SYMBOL) %>%
    summarise(
      a = n_distinct(sample_id[group == "Case"]),  # Cases carrying the variant
      c = n_distinct(sample_id[group == "Control"]),  # Controls carrying the variant
      total_cases = 302,  # Total number of cases (fixed)
      total_controls = 75  # Total number of controls (fixed)
    ) %>%
    ungroup()

  # Step 2: Calculate Fisher's Exact Test p-value and Odds Ratio
  result_df <- result_df %>%
    rowwise() %>%
    mutate(
      b = total_cases - a,  # Cases not carrying the variant (total_cases - a)
      d = total_controls - c,  # Controls not carrying the variant (total_controls - c)
      p_value = fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value,
      odds_ratio = (a / b) / (c / d),  # Odds ratio calculation
      # Calculate 95% CI for the odds ratio
      lower_ci = exp(log(odds_ratio) - 1.96 * sqrt(1/a + 1/b + 1/c + 1/d)),
      upper_ci = exp(log(odds_ratio) + 1.96 * sqrt(1/a + 1/b + 1/c + 1/d))
    ) %>%
    ungroup()

  # Step 3: Create a new dataframe with the new columns (keeping original data intact)
  final_result_df <- df %>%
    left_join(result_df, by = c("CHROM", "POS", "SYMBOL")) %>%
    select(
      CHROM, POS, SYMBOL, sample_id, group, 
      a, b, c, d, p_value, total_cases, total_controls, odds_ratio, lower_ci, upper_ci
    )

  return(final_result_df)
}

# Step 2: Apply the analysis function to each dataset
result_df_NON_patho_pLI_new <- perform_analysis(result_df_NON_patho_pLI)
result_df_patho_pLI_new <- perform_analysis(result_df_patho_pLI)
result_df_patho_new <- perform_analysis(result_df_patho)
result_df_PTV_pLI_new <- perform_analysis(result_df_PTV_pLI)
result_df_PTV_new <- perform_analysis(result_df_PTV)

# Step 3: Save the new dataframes to CSV files in the specified directory
write.csv(result_df_NON_patho_pLI_new, "/home/rachele/vcf_csv/result_df_NON_patho_pLI_new.csv", row.names = FALSE)
write.csv(result_df_patho_pLI_new, "/home/rachele/vcf_csv/result_df_patho_pLI_new.csv", row.names = FALSE)
write.csv(result_df_patho_new, "/home/rachele/vcf_csv/result_df_patho_new.csv", row.names = FALSE)
write.csv(result_df_PTV_pLI_new, "/home/rachele/vcf_csv/result_df_PTV_pLI_new.csv", row.names = FALSE)
write.csv(result_df_PTV_new, "/home/rachele/vcf_csv/result_df_PTV_new.csv", row.names = FALSE)

# Optional: Print confirmation
message("New result dataframes have been saved successfully in /home/rachele/vcf_csv/")
