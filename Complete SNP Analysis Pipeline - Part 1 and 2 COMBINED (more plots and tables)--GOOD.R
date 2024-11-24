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

create_basic_summary_plots <- function(snp_stats) {
  # 1. Missing Rate Distribution
  p1 <- ggplot(snp_stats, aes(x = missing_rate)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 0.2, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Missing Rate",
         subtitle = "Red line: 20% threshold",
         x = "Missing Rate",
         y = "Count") +
    theme_minimal()
  
  # 2. Signal Variation
  p2 <- ggplot(snp_stats, aes(x = sd_theta)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 0.2, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Theta Standard Deviation",
         subtitle = "Red line: 0.2 threshold",
         x = "Theta SD",
         y = "Count") +
    theme_minimal()
  
  # 3. Cluster Separation
  p3 <- ggplot(snp_stats, aes(x = cluster_sep)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 0.1, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Cluster Separation",
         subtitle = "Red line: 0.1 threshold",
         x = "Cluster Separation",
         y = "Count") +
    theme_minimal()
  
  # 4. Heterozygosity Rate
  p4 <- ggplot(snp_stats, aes(x = het_rate)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 0.9, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Heterozygosity Rate",
         subtitle = "Red line: 90% threshold",
         x = "Heterozygosity Rate",
         y = "Count") +
    theme_minimal()
  
  # 5. Signal Intensity Variation
  p5 <- ggplot(snp_stats, aes(x = sd_R)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 1.0, color = "red", linetype = "dashed") +
    labs(title = "Distribution of R Standard Deviation",
         subtitle = "Red line: 1.0 threshold",
         x = "R SD",
         y = "Count") +
    theme_minimal()
  
  # Save plots
  ggsave(file.path(output_dir, "Figures", "missing_rate_dist_basic.png"), p1, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "theta_sd_dist_basic.png"), p2, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "cluster_sep_dist_basic.png"), p3, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "het_rate_dist_basic.png"), p4, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "r_sd_dist_basic.png"), p5, width = 8, height = 6)
  
  return(list(
    missing_rate_dist = p1,
    theta_sd_dist = p2,
    cluster_sep_dist = p3,
    het_rate_dist = p4,
    r_sd_dist = p5
  ))
}

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
  
  # Multiple issue breakdown
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
  
  # Signal Variation Patterns
  p3 <- failure_analysis$single_issues %>%
    ggplot(aes(x = sd_theta, y = sd_R, color = failure_category)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 1.0, linetype = "dashed", color = "red") +
    labs(title = "Signal Variation Patterns",
         x = "Theta SD",
         y = "R SD") +
    theme_minimal()
  
  # Missing vs Heterozygosity
  p4 <- failure_analysis$single_issues %>%
    ggplot(aes(x = missing_rate, y = het_rate, color = failure_category)) +
    geom_point(alpha = 0.6) +
    geom_vline(xintercept = 0.2, linetype = "dashed", color = "red") +
    geom_hline(yintercept = 0.9, linetype = "dashed", color = "red") +
    labs(title = "Missing Rate vs Heterozygosity",
         x = "Missing Rate",
         y = "Heterozygosity Rate") +
    theme_minimal()
  
  # Chromosome distribution
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
  
  # Create summary tables
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

plot_snp_clusters <- function(snp_long, snp_stats, n_examples = 5) {
  # Sample problematic SNPs with different issues
  problem_snps <- snp_stats %>%
    filter(is_problematic) %>%
    mutate(
      issue = case_when(
        flag_missing ~ "High Missing Rate",
        flag_theta_var ~ "High Theta Variation",
        flag_R_var ~ "High R Variation",
        flag_cluster ~ "Poor Cluster Separation",
        flag_het ~ "High Heterozygosity",
        TRUE ~ "Multiple Issues"
      )
    ) %>%
    group_by(issue) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  # Create cluster plots
  for(i in 1:nrow(problem_snps)) {
    snp <- problem_snps$Name[i]
    issue <- problem_snps$issue[i]
    
    snp_data <- snp_long %>%
      filter(Name == snp)
    
    p <- ggplot(snp_data, aes(x = Theta, y = R, color = GT)) +
      geom_point(alpha = 0.6) +
      scale_color_viridis_d() +
      labs(title = paste("Cluster Plot for", snp),
           subtitle = paste("Issue:", issue, "\nN =", nrow(snp_data)),
           x = "Theta",
           y = "R") +
      theme_minimal()
    
    ggsave(file.path(output_dir, "Figures", 
                     paste0("cluster_plot_", snp, "_", 
                            gsub(" ", "_", issue), ".png")), 
           p, width = 8, height = 6)
  }
}


# Main analysis function
main_analysis <- function(snp_file) {
  # Process data
  snp_data <- process_snp_data(snp_file)
  
  # Analyze SNP stability
  snp_stability <- analyze_snp_stability(snp_data)
  
  # Generate basic summary plots
  basic_plots <- create_basic_summary_plots(snp_stability)
  
  # Analyze failure patterns and create additional plots
  failure_patterns <- analyze_failure_patterns(snp_stability)
  failure_plots <- plot_failure_distributions(failure_patterns, output_dir)
  
  # Create cluster plots for problematic SNPs
  plot_snp_clusters(snp_data, snp_stability)
  
  # Save results
  write_csv(snp_stability, file.path(output_dir, "Tables", "snp_stability_metrics.csv"))
  write_csv(filter(snp_stability, is_problematic), 
            file.path(output_dir, "Tables", "problematic_snps.csv"))
  
  # Return all results
  return(list(
    snp_stability = snp_stability,
    basic_plots = basic_plots,
    failure_patterns = failure_patterns,
    failure_plots = failure_plots
  ))
}

# Run the analysis
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
