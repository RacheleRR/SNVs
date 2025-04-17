#filtering for new SNV filter criteria 


library(dplyr)
library(readr)
library(stringr)



filtered_df<- merged_VEP_output_final %>%
mutate(
gnomADe_AF = parse_number(gnomADe_AF, na = c("", "NA", ".", "None")),
gnomADg_AF = parse_number(gnomADg_AF, na = c("", "NA", ".", "None"))
) %>%
filter(
(is.na(gnomADe_AF) | gnomADe_AF <= 0.01),
(is.na(gnomADg_AF) | gnomADg_AF <= 0.01)
)


filtered_df_1 <- filtered_df %>%
    mutate(
        MPC_score = str_remove_all(MPC_score, "\\.&|&|\\.$"),
        MPC_score = parse_number(MPC_score, na = c("", "NA", ".")),
        am_pathogenicity = str_remove_all(am_pathogenicity, "\\.&|&|\\.$"),
        am_pathogenicity = parse_number(am_pathogenicity, na = c("", "NA", ".")),
        pLI_gene_value = str_remove_all(pLI_gene_value, "\\.&|&|\\.$"),
        pLI_gene_value = parse_number(pLI_gene_value, na = c("", "NA", ".")),
        CADD_phred = str_remove_all(CADD_phred, "\\.&|&|\\.$"),
        CADD_phred = parse_number(CADD_phred, na = c("", "NA", ".")),
        CADD_raw = str_remove_all(CADD_raw, "\\.&|&|\\.$"), 
        CADD_raw = parse_number(CADD_raw, na = c("", "NA", ".")),
        MVP_score = str_remove_all(MVP_score, "\\.&|&|\\.$"),
        MVP_score = parse_number(MVP_score, na = c("", "NA", ".")),
        MetaSVM_score = str_remove_all(MetaSVM_score, "\\.&|&|\\.$"),
        MetaSVM_score = parse_number(MetaSVM_score, na = c("", "NA", ".")),
        MetaSVM_pred = str_remove_all(MetaSVM_pred, "\\.&|&|\\.$"),  # Just clean, no parse_number()
        ClinPred_score = str_remove_all(ClinPred_score, "\\.&|&|\\.$"),
        ClinPred_score = parse_number(ClinPred_score, na = c("", "NA", "."))
    )



write.csv(filtered_df_1, "filtered_df_1.csv", row.names = FALSE)
write.csv(filtered_df, "filtered_df.csv", row.names = FALSE)


manifest_correct <- read.delim("~/SNVs/results/euro/manifest_correct.tsv")

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


mpc <-  filtered_df_1 %>%
filter(
MPC_score >= 2,
!is.na(MPC_score))

alpha <- filtered_df_1 %>%
filter(
am_pathogenicity >= 0.98,
!is.na(am_pathogenicity)
)

CADD_phred <- filtered_df_1 %>%
filter(
CADD_phred >= 20,
!is.na(CADD_phred)
)

MVP_score <- filtered_df_1 %>%
filter(
MVP_score >= 0.7,
!is.na(MVP_score)
)
MetaSVM_score <- filtered_df_1 %>%
filter(
MetaSVM_score >= 0,
!is.na(MetaSVM_score)
)
MetaSVM_pred <- filtered_df_1 %>%
filter(
MetaSVM_pred == "D",
!is.na(MetaSVM_pred)
)

ClinPred_score <- filtered_df_1 %>%
filter(
ClinPred_score >= 0.5,
!is.na(ClinPred_score)
)

CADD_phred_stringent <- filtered_df_1 %>%
filter(
CADD_phred >= 30,
!is.na(CADD_phred)
)



process_df <- function(df, manifest_correct) {
  df$SAMPLES <- gsub("_pool", "", df$SAMPLES)
  df$sample_label_2 <- sapply(df$SAMPLES, get_outlier_labels, manifest_correct)
  df$sample_label <- sapply(df$SAMPLES, check_outlier_label, manifest_correct)
  df <- clean_uh_na(df)
  df$sample_label_3 <- sapply(df$sample_label_2, derive_group_label)
  df <- df %>% mutate(count = sapply(strsplit(SAMPLES, ","), length))
  df <- df %>% filter(count == 1)
  df <- df %>% distinct(SYMBOL,ID, CHROM,REF,ALT, POS, .keep_all = TRUE)
  return(df)
}

mpc <- process_df(mpc, manifest_correct)
alpha <- process_df(alpha, manifest_correct)
CADD_phred <- process_df(CADD_phred, manifest_correct)
MVP_score <- process_df(MVP_score, manifest_correct)
MetaSVM_score <- process_df(MetaSVM_score, manifest_correct)
MetaSVM_pred <- process_df(MetaSVM_pred, manifest_correct)
ClinPred_score <- process_df(ClinPred_score, manifest_correct)
CADD_phred_stringent <- process_df(CADD_phred_stringent, manifest_correct)


if (!dir.exists("filter_snv_manual_new_scores")) {
dir.create("filter_snv_manual_new_scores")
}

setwd("filter_snv_manual_new_scores")
# Save each filtered dataframe
write.csv(mpc, "mpc.csv", row.names = FALSE)
write.csv(alpha, "alpha.csv", row.names = FALSE)
write.csv(CADD_phred, "CADD_phred.csv", row.names = FALSE)
write.csv(MVP_score, "MVP_score.csv", row.names = FALSE)
write.csv(MetaSVM_score, "MetaSVM_score.csv", row.names = FALSE)
write.csv(MetaSVM_pred, "MetaSVM_pred.csv", row.names = FALSE)
write.csv(ClinPred_score, "ClinPred_score.csv", row.names = FALSE)
write.csv(CADD_phred_stringent, "CADD_phred_stringent.csv", row.names = FALSE)











patho_mpc_alpha <- filtered_df_1 %>%
filter(
MPC_score >= 2,
am_pathogenicity >= 0.98,
# Optional: Explicitly exclude NAs if needed
!is.na(MPC_score),
!is.na(am_pathogenicity)
)


non_patho_mpc <- filtered_df_1 %>%
filter(
MPC_score < 2,
!is.na(MPC_score)  # Exclude NAs
)

ptv_variants <- filtered_df_1 %>%
filter(
Consequence %in% c(
"frameshift_variant",
"stop_gained",
"splice_acceptor_variant",
"splice_donor_variant",
"start_lost"
)
)


patho_mpc_alpha_pli <- filtered_df_1 %>%
filter(
pLI_gene_value > 0.9,
MPC_score >= 2,
am_pathogenicity >= 0.98,
!is.na(pLI_gene_value),
!is.na(MPC_score),
!is.na(am_pathogenicity)
)

non_patho_mpc_pli <- filtered_df_1 %>%
filter(
pLI_gene_value > 0.9,
MPC_score < 2,
!is.na(pLI_gene_value),
!is.na(MPC_score)
)

ptv_pli <- filtered_df_1 %>%
filter(
pLI_gene_value > 0.9,
Consequence %in% c(
"frameshift_variant",
"stop_gained",
"splice_acceptor_variant",
"splice_donor_variant",
"start_lost"
),
!is.na(pLI_gene_value)
)


if (!dir.exists("filter_snv_manual_sample")) {
dir.create("filter_snv_manual_sample")
}
# Save each filtered dataframe
write_tsv(patho_mpc_alpha, "filter_snv_manual_sample/patho_mpc_alpha.tsv")
write_tsv(non_patho_mpc, "filter_snv_manual_sample/non_patho_mpc.tsv")
write_tsv(ptv_variants, "filter_snv_manual_sample/ptv_variants.tsv")
write_tsv(patho_mpc_alpha_pli, "filter_snv_manual_sample/patho_mpc_alpha_pli.tsv")
write_tsv(non_patho_mpc_pli, "filter_snv_manual_sample/non_patho_mpc_pli.tsv")
write_tsv(ptv_pli, "filter_snv_manual_sample/ptv_pli.tsv")