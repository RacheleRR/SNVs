# tabel with numkber_SNVs
#try one
expanded_df <- pred_NON_patho_MPC_pLI_VEP_filtered %>%
    separate_rows(SAMPLES, sample_label_2, sep = ",") %>%
    mutate(
        sample_id = trimws(SAMPLES),
        group = trimws(sample_label_2)
    ) %>%
    distinct()


num_non_pathogenic_pLi_variants <- expanded_df %>%
    group_by(group) %>%
    summarize(num_variants = n_distinct(POS), .groups = "drop")

num_individuals_with_non_pathogenic_pLi <- expanded_df %>%
    group_by(group) %>%
    summarize(unique_individuals = n_distinct(sample_id), .groups = "drop")

result_df <- data.frame(
    row_names = c("Number of NON pathogenic variants", 
                  "Number of individuals with NON pathogenic variants"),
    
    # Get the counts for each group (Case, Control, UHR_NA)
    Case = c(
        num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Case"],
        num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "Case"]
    ),
    
    Control = c(
        num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Control"],
        num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "Control"]
    ),
    
    UHR_NA = c(
        num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "UHR_NA"],
        num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "UHR_NA"]
    )
)

# Define a function to process each dataframe and return the result
process_df <- function(df, name) {
  
  expanded_df <- df %>%
    separate_rows(SAMPLES, sample_label_2, sep = ",") %>%
    mutate(
        sample_id = trimws(SAMPLES),
        group = trimws(sample_label_2)
    ) %>%
    distinct()
  
  num_non_pathogenic_pLi_variants <- expanded_df %>%
    group_by(group) %>%
    summarize(num_variants = n_distinct(POS), .groups = "drop")
  
  num_individuals_with_non_pathogenic_pLi <- expanded_df %>%
    group_by(group) %>%
    summarize(unique_individuals = n_distinct(sample_id), .groups = "drop")




# Define the row names based on the dataframe's name
  row_names <- c(
    paste("Number of", name ,"variants"),
    paste("Number of individuals with",name ,"variants" )
  )
  
  # Create the result dataframe
  result_df <- data.frame(
    row_names = row_names,
    
    Case = c(
        num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Case"],
        num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "Case"]
    ),
    
    Control = c(
        num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "Control"],
        num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "Control"]
    ),
    
    UHR_NA = c(
        num_non_pathogenic_pLi_variants$num_variants[num_non_pathogenic_pLi_variants$group == "UHR_NA"],
        num_individuals_with_non_pathogenic_pLi$unique_individuals[num_individuals_with_non_pathogenic_pLi$group == "UHR_NA"]
    )
  )
  
  return(result_df)
}

# Process each dataframe
result_df_list <- list(
  process_df(pred_NON_patho_MPC_pLI_VEP_filtered, "NON_patho_pLI"),
  process_df(pred_patho_MPC_ALPHAMISSENSE_pLI_VEP_filtered, "patho_pLI"),
  process_df(pred_patho_MPC_ALPHAMISSENSE_VEP_filtered, "patho"),
  process_df(pred_PTV_pLI_VEP_filtered, "PTV_pLI"),
  process_df(pred_PTV_VEP_filtered, "PTV")
)

# Combine all result dataframes into one
final_result_df <- bind_rows(result_df_list)

# View the final result
print(final_result_df)