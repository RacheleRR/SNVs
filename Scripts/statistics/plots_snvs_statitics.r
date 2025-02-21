plot_fisher_volcano <- function(test_results_df) {
    ggplot(test_results_df, aes(x = fisher_odds_ratio, y = -log10(fisher_p_value), color = Test)) +
        geom_point(size = 3) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
        labs(title = "Fisher Test Volcano Plot", x = "Odds Ratio", y = "-log10(p-value)") +
        theme_minimal()
}
plot_fisher_volcano(fisher_results_ind_df)
plot_fisher_volcano(result_var_conf_df)
plot_fisher_volcano(result_var_conf_norm_df)
plot_fisher_volcano(result_var_conf_enrich_df)

plot_fisher_forest <- function(test_results_df) {
    ggplot(test_results_df, aes(x = fisher_odds_ratio, y = Row, color = Test)) +
        geom_point(size = 3) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
        labs(title = "Forest Plot of Fisher Odds Ratios", x = "Odds Ratio", y = "Row") +
        theme_minimal()
}
plot_fisher_forest(fisher_results_ind_df)
plot_fisher_forest(result_var_conf_df)
plot_fisher_forest(result_var_conf_norm_df)
plot_fisher_forest(result_var_conf_enrich_df)

library(ggplot2)

# Assuming test_results_df is the data frame containing the Chi-Square test results
# If not already done, create the test_results_df as shown in the previous code



# Create the dot plot
plot_pval_dot <- function(test_results_df, Row, p_value, Test) {
ggplot(test_results_df, aes(x = Row, y = p_value, color = Test)) +
  geom_point(size = 3) +
  labs(title = "Dot Plot of Test P-Values",
       x = "Variant Type",
       y = "P-Value",
       color = "Test") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red")  # Add a line for significance threshold
}

plot_pval_dot(result_var_conf_df, result_var_conf_df$Row, result_var_conf_df$chi_squared_p_value, result_var_conf_df$Test)
plot_pval_dot(wilcoxon_results_df_no_lev, wilcoxon_results_df_no_lev$Column, wilcoxon_results_df_no_lev$p_value_wilcox, wilcoxon_results_df_no_lev$DataFrame)
plot_pval_dot(result_var_conf_enrich_df,result_var_conf_enrich_df$Row,result_var_conf_enrich_df$chi_squared_p_value,result_var_conf_enrich_df$Test)