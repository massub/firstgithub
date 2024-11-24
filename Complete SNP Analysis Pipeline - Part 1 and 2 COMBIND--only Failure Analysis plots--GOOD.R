# Required packages
library(tidyverse)
library(readr)
library(ggplot2)
library(viridis)

# File path
snp_file <- "D:/USDA/Small Grains Research Unit/Projects/Project--SNP position and cluster charecterization/Project_Root/00_Raw_Data/Oat_Datasets/SDSU2023Caffe-O3K_01-08-FDT.txt"

# Output directories
output_dir <- "D:/USDA/Small Grains Research Unit/Projects/Project--SNP position and cluster charecterization/Project_Root/03_Results"
dir.create(file.path(output_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)

# Data processing function
process_snp_data <- function(file_path) {
  message("Reading SNP data...")
  snp_data <- read_delim(file_path, delim = "\t", show_col_types = FALSE)
  
  # Get all column names
  all_cols <- names(snp_data)
  
  # Extract unique sample names from GT columns
  sample_names <- unique(sub("\\.GT.*$", "", grep("\\.GT", all_cols, value = TRUE)))
  n_samples <- length(sample_names)
  
  message("Number of SNPs: ", nrow(snp_data))
  message("Number of samples: ", n_samples)
  
  # Initialize result dataframe
  result <- tibble(
    Index = rep(snp_data$Index, times = n_samples),
    Name = rep(snp_data$Name, times = n_samples),
    Address = rep(snp_data$Address, times = n_samples),
    Chr = rep(snp_data$Chr, times = n_samples),
    Position = rep(snp_data$Position, times = n_samples),
    Sample = rep(sample_names, each = nrow(snp_data)),
    GT = NA_character_,
    Theta = NA_real_,
    R = NA_real_
  )
  
  # Process each sample
  for(i in seq_along(sample_names)) {
    sample <- sample_names[i]
    start_idx <- ((i-1) * nrow(snp_data) + 1)
    end_idx <- i * nrow(snp_data)
    
    # Get column names for this sample
    gt_col <- grep(paste0("^", sample, "\\.GT"), all_cols, value = TRUE)[1]
    theta_col <- grep(paste0("^", sample, "\\.Theta"), all_cols, value = TRUE)[1]
    r_col <- grep(paste0("^", sample, "\\.R"), all_cols, value = TRUE)[1]
    
    # Fill data
    result$GT[start_idx:end_idx] <- snp_data[[gt_col]]
    result$Theta[start_idx:end_idx] <- snp_data[[theta_col]]
    result$R[start_idx:end_idx] <- snp_data[[r_col]]
  }
  
  return(result)
}

# SNP stability analysis function
analyze_snp_stability <- function(snp_long) {
  message("Analyzing SNP stability...")
  
  snp_metrics <- snp_long %>%
    group_by(Name, Chr, Position) %>%
    summarise(
      # Basic counts
      n_samples = n(),
      n_AA = sum(GT == "AA", na.rm = TRUE),
      n_AB = sum(GT == "AB", na.rm = TRUE),
      n_BB = sum(GT == "BB", na.rm = TRUE),
      n_NC = sum(is.na(GT)),
      
      # Signal metrics
      mean_theta = mean(Theta, na.rm = TRUE),
      sd_theta = sd(Theta, na.rm = TRUE),
      mean_R = mean(R, na.rm = TRUE),
      sd_R = sd(R, na.rm = TRUE),
      
      # Cluster metrics
      theta_AA = mean(Theta[GT == "AA"], na.rm = TRUE),
      theta_AB = mean(Theta[GT == "AB"], na.rm = TRUE),
      theta_BB = mean(Theta[GT == "BB"], na.rm = TRUE),
      
      R_AA = mean(R[GT == "AA"], na.rm = TRUE),
      R_AB = mean(R[GT == "AB"], na.rm = TRUE),
      R_BB = mean(R[GT == "BB"], na.rm = TRUE),
      
      # Quality metrics
      missing_rate = n_NC / n_samples,
      
      # Calculate cluster separation
      cluster_sep = case_when(
        all(c(n_AA, n_BB) > 0) ~ abs(theta_AA - theta_BB),
        all(c(n_AA, n_AB) > 0) ~ abs(theta_AA - theta_AB),
        all(c(n_AB, n_BB) > 0) ~ abs(theta_AB - theta_BB),
        TRUE ~ NA_real_
      ),
      
      # Population metrics
      maf = pmin(n_AA, n_BB) / (n_AA + n_AB + n_BB),
      het_rate = n_AB / (n_AA + n_AB + n_BB),
      
      .groups = "drop"
    ) %>%
    mutate(
      # Quality flags
      flag_missing = missing_rate > 0.2,
      flag_theta_var = sd_theta > 0.2,
      flag_R_var = sd_R > 1.0,
      flag_cluster = (cluster_sep < 0.1 & !is.na(cluster_sep)),
      flag_het = (het_rate > 0.9 & !is.na(het_rate)),
      
      # Overall problematic status
      is_problematic = flag_missing | flag_theta_var | flag_R_var | 
        flag_cluster | flag_het
    )
  
  return(snp_metrics)
}

# Analysis and visualization functions
analyze_failure_patterns <- function(snp_stability) {
  # Single issue analysis
  single_issue_patterns <- snp_stability %>%
    mutate(
      failure_category = case_when(
        !is_problematic ~ "Pass",
        flag_missing & !flag_theta_var & !flag_R_var & !flag_cluster & !flag_het ~ "Missing Rate Only",
        !flag_missing & flag_theta_var & !flag_R_var & !flag_cluster & !flag_het ~ "Theta Variation Only",
        !flag_missing & !flag_theta_var & flag_R_var & !flag_cluster & !flag_het ~ "R Variation Only",
        !flag_missing & !flag_theta_var & !flag_R_var & flag_cluster & !flag_het ~ "Cluster Separation Only",
        !flag_missing & !flag_theta_var & !flag_R_var & !flag_cluster & flag_het ~ "Heterozygosity Only",
        TRUE ~ "Multiple Issues"
      )
    )
  
  # Multiple issue breakdown - fixed vectorization
  multiple_issues <- snp_stability %>%
    filter(is_problematic) %>%
    mutate(
      issues = paste0(
        ifelse(flag_missing, "Missing;", ""),
        ifelse(flag_theta_var, "Theta;", ""),
        ifelse(flag_R_var, "R;", ""),
        ifelse(flag_cluster, "Cluster;", ""),
        ifelse(flag_het, "Het;", "")
      )
    ) %>%
    mutate(issues = str_remove(issues, ";$"))
  
  return(list(
    single_issues = single_issue_patterns,
    multiple_issues = multiple_issues
  ))
}

plot_failure_distributions <- function(failure_analysis, output_dir) {
  # Overall failure distribution
  p1 <- failure_analysis$single_issues %>%
    ggplot(aes(x = failure_category, fill = failure_category)) +
    geom_bar() +
    coord_flip() +
    labs(title = "Distribution of SNP Quality Issues",
         x = "Category",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Multiple issues combinations
  p2 <- failure_analysis$multiple_issues %>%
    count(issues) %>%
    arrange(desc(n)) %>%
    slice_head(n = 10) %>%
    ggplot(aes(x = reorder(issues, n), y = n)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Top 10 Multiple Issue Combinations",
         x = "Issue Combination",
         y = "Count") +
    theme_minimal()
  
  # Relationship plots
  p3 <- failure_analysis$single_issues %>%
    ggplot(aes(x = sd_theta, y = sd_R, color = failure_category)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
    labs(title = "Signal Variation Patterns",
         x = "Theta SD",
         y = "R SD") +
    theme_minimal()
  
  p4 <- failure_analysis$single_issues %>%
    ggplot(aes(x = missing_rate, y = het_rate, color = failure_category)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    labs(title = "Missing Rate vs Heterozygosity",
         x = "Missing Rate",
         y = "Heterozygosity Rate") +
    theme_minimal()
  
  p5 <- failure_analysis$single_issues %>%
    ggplot(aes(x = Chr, fill = failure_category)) +
    geom_bar(position = "fill") +
    labs(title = "Proportion of Issues by Chromosome",
         x = "Chromosome",
         y = "Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save plots
  ggsave(file.path(output_dir, "Figures", "failure_distribution.png"), p1, width = 10, height = 8)
  ggsave(file.path(output_dir, "Figures", "multiple_issues.png"), p2, width = 10, height = 8)
  ggsave(file.path(output_dir, "Figures", "signal_variation_patterns.png"), p3, width = 10, height = 8)
  ggsave(file.path(output_dir, "Figures", "missing_vs_het.png"), p4, width = 10, height = 8)
  ggsave(file.path(output_dir, "Figures", "chromosome_distribution.png"), p5, width = 12, height = 8)
  
  # Create and save summary tables
  failure_summary <- failure_analysis$single_issues %>%
    count(failure_category) %>%
    mutate(percentage = n/sum(n) * 100)
  
  multiple_issue_summary <- failure_analysis$multiple_issues %>%
    count(issues) %>%
    arrange(desc(n)) %>%
    mutate(percentage = n/sum(n) * 100)
  
  write_csv(failure_summary, file.path(output_dir, "Tables", "failure_category_summary.csv"))
  write_csv(multiple_issue_summary, file.path(output_dir, "Tables", "multiple_issues_summary.csv"))
  
  return(list(
    overall_distribution = p1,
    multiple_issues = p2,
    signal_patterns = p3,
    missing_het_relation = p4,
    chromosome_distribution = p5,
    failure_summary = failure_summary,
    multiple_issue_summary = multiple_issue_summary
  ))
}

# Main analysis function
main_analysis <- function(snp_file) {
  # Process data
  snp_data <- process_snp_data(snp_file)
  
  # Analyze SNP stability
  snp_stability <- analyze_snp_stability(snp_data)
  
  # Analyze failure patterns
  failure_patterns <- analyze_failure_patterns(snp_stability)
  failure_plots <- plot_failure_distributions(failure_patterns, output_dir)
  
  # Save results
  write_csv(snp_stability, file.path(output_dir, "Tables", "snp_stability_metrics.csv"))
  write_csv(filter(snp_stability, is_problematic), 
            file.path(output_dir, "Tables", "problematic_snps.csv"))
  
  return(list(
    snp_stability = snp_stability,
    failure_patterns = failure_patterns,
    failure_plots = failure_plots
  ))
}

# Run analysis
results <- try({
  main_analysis(snp_file)
})

if(inherits(results, "try-error")) {
  message("Error in analysis. Please check the error message above.")
} else {
  message("Analysis completed successfully!")
  message("\nFailure Summary:")
  print(results$failure_patterns$single_issues %>% 
          count(failure_category) %>% 
          mutate(percentage = n/sum(n) * 100))
}
