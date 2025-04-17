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
library(tools)


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
    folder_path <- "/home/rachele/filter_snv_manual_new_scores"

# List all TSV and CSV files in the folder
all_files <- list.files(path = folder_path, pattern = "\\.(tsv|csv)$", full.names = TRUE)

# Loop through each file and assign it to an object
for (file in all_files) {
  file_name <- file_path_sans_ext(basename(file))
  
  # Choose correct reader based on file extension
  if (grepl("\\.tsv$", file, ignore.case = TRUE)) {
    df <- read_tsv(file)
  } else if (grepl("\\.csv$", file, ignore.case = TRUE)) {
    df <- read_csv(file)
  } else {
    next  # Skip any files that aren't .tsv or .csv
  }
  
  assign(file_name, df)
}


    UHR_NA_SAMPLE_IDS <- read.csv("~/UHR_NA_SAMPLE_IDS.csv", sep="")
    #example S36827

colnames(manifest_correct) <- c("Status", "sample_id", "group")


# stuff
process_df <- function(df, name, brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or) {

    
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
        )                         # Remove duplicates

    return(result_df)
}

result_mpc <- process_df(mpc, "mpc", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_alpha <- process_df(alpha, "alpha", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_CADD_phred <- process_df(CADD_phred, "cadd_phred", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_MVP_score <- process_df(MVP_score, "MVP_score", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_MetaSVM_score <- process_df(MetaSVM_score, "MetaSVM_score", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_MetaSVM_pred <- process_df(MetaSVM_pred, "MetaSVM_pred", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_ClinPred_score <- process_df(ClinPred_score, "ClinPred_score", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)
result_CADD_phred_stringent <- process_df(CADD_phred_stringent, "CADD_phred_stringent", brain_gene_consensus_filtered_consensus_no_pitular, brain_gene_consensus_ntm_consensus_no_pitular, genes_schema_pval, genes_bipolar, genes_sfari_non_syndromic, genes_schema_or)

manifest_correct <- manifest_correct %>% filter(!grepl("UHR_NA", group))

# ADD ALL SAMPLES INDEPENDENTLY OF PRESENCE/ABSENCE OF VARIANTS 
# Ensure consistency in column names and perform left join for all dataframes
complete_mpc <- manifest_correct %>%
    left_join(result_mpc, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))

complete_alpha <- manifest_correct %>%
    left_join(result_alpha, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))
complete_CADD_phred <- manifest_correct %>%
    left_join(result_CADD_phred, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))
complete_MVP_score <- manifest_correct %>%
    left_join(result_MVP_score, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))
complete_MetaSVM_score <- manifest_correct %>%
    left_join(result_MetaSVM_score, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))
complete_MetaSVM_pred <- manifest_correct %>%
    left_join(result_MetaSVM_pred, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))
complete_ClinPred_score <- manifest_correct %>%
    left_join(result_ClinPred_score, by = c("sample_id", "group")) %>%
    mutate(across(starts_with("Number_of"), ~ ifelse(is.na(.), 0, .)))
complete_CADD_phred_stringent <- manifest_correct %>%
    left_join(result_CADD_phred_stringent, by = c("sample_id", "group")) %>%
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
    wil_mpc <- analyze_metrics(complete_mpc) %>%  select(-case_zero, -control_zero)
    wil_alpha <- analyze_metrics(complete_alpha) %>%  select(-case_zero, -control_zero)
    wil_CADD_phred <- analyze_metrics(complete_CADD_phred) %>%  select(-case_zero, -control_zero)
    wil_MVP_score <- analyze_metrics(complete_MVP_score) %>%  select(-case_zero, -control_zero)
    wil_ClinPred_score <- analyze_metrics(complete_ClinPred_score) %>%  select(-case_zero, -control_zero)
    wil_MetaSVM_score <- analyze_metrics(complete_MetaSVM_score) %>%  select(-case_zero, -control_zero)
    wil_MetaSVM_pred <- analyze_metrics(complete_MetaSVM_pred) %>%  select(-case_zero, -control_zero)
    wil_CADD_phred_stringent <- analyze_metrics(complete_CADD_phred_stringent) %>%  select(-case_zero, -control_zero)




# VISUALIZATION
setwd("/home/rachele/SNVs/results/stats/new_scores")
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

            plot_mpc <- plot(
            complete_mpc ,
            wil_mpc,
            "MPC"
            )

            plot_alpha <- plot(
            complete_alpha, 
            wil_alpha,
            "Alpha"
            )
            plot_CADD_phred <- plot(
            complete_CADD_phred,
            wil_CADD_phred,
            "CADD_phred"
            )
            plot_MVP_score <- plot(
            complete_MVP_score,
            wil_MVP_score,
            "MVP_score"
            )
            plot_MetaSVM_score <- plot(
            complete_MetaSVM_score,
            wil_MetaSVM_score,
            "MetaSVM_score"
            )
            plot_MetaSVM_pred <- plot(
            complete_MetaSVM_pred,
            wil_MetaSVM_pred,
            "MetaSVM_pred"
            )
            plot_ClinPred_score <- plot(
            complete_ClinPred_score,
            wil_ClinPred_score,
            "ClinPred_score"
            )
            plot_CADD_phred_stringent <- plot(
            complete_CADD_phred_stringent,
            wil_CADD_phred_stringent,
            "CADD_phred_stringent"
            )

            # To save all plots:
            plots <- list( plot_mpc, plot_alpha, plot_CADD_phred, plot_MVP_score, plot_MetaSVM_score, plot_MetaSVM_pred, plot_ClinPred_score, plot_CADD_phred_stringent)


            walk2(plots, c("mpc","alpha", "CADD_phred", "MVP_score", "MetaSVM_score", "MetaSVM_pred", "ClinPred_score", "CADD_phred_stringent"), 
                ~ ggsave(paste0("variant_comparison_", .y, ".png"), .x, width = 10, height = 8))


        #save tabels 
        combined_results_non_norm <- bind_rows(
        "mpc" = wil_mpc,
        "alpha" = wil_alpha,
        "CADD_phred" = wil_CADD_phred,
        "MVP_score" = wil_MVP_score,
        "MetaSVM_score" = wil_MetaSVM_score,
        "MetaSVM_pred" = wil_MetaSVM_pred,
        "ClinPred_score" = wil_ClinPred_score,
        "CADD_phred_stringent" = wil_CADD_phred_stringent,
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
    plot_significant(complete_mpc, wil_mpc, "MPC"),
    plot_significant(complete_alpha, wil_alpha, "Alpha"),
    plot_significant(complete_CADD_phred, wil_CADD_phred, "CADD_phred"),
    plot_significant(complete_MVP_score, wil_MVP_score, "MVP_score"),
    plot_significant(complete_MetaSVM_score, wil_MetaSVM_score, "MetaSVM_score"),
    plot_significant(complete_MetaSVM_pred, wil_MetaSVM_pred, "MetaSVM_pred"),
    plot_significant(complete_ClinPred_score, wil_ClinPred_score, "ClinPred_score"),
    plot_significant(complete_CADD_phred_stringent, wil_CADD_phred_stringent, "CADD_phred_stringent")
  
    ) %>% 
    compact()  # Remove NULL entries (non-significant results)

    # Save significant plots
    walk2(
    sig_plots, 
    c("mpc","alpha", "CADD_phred", "MVP_score", "MetaSVM_score", "MetaSVM_pred", "ClinPred_score", "CADD_phred_stringent")[1:length(sig_plots)],
    ~ ggsave(
        paste0("SIGNIFICANT_variant_comparison_", .y, ".png"), 
        .x, 
        width = 10, 
        height = 6
    )
    )


# Save the combined results to a CSV file
write.csv(combined_results_non_norm, "combined_results_non_norm.csv", row.names = FALSE)

