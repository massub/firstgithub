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
      # Updated criteria for problematic SNPs
      is_problematic = missing_rate > 0.2 |           # >20% missing data
        sd_theta > 0.2 |                # High theta variation
        sd_R > 1.0 |                    # High R variation
        (cluster_sep < 0.1 & !is.na(cluster_sep)) | # Poor separation
        (het_rate > 0.9 & !is.na(het_rate)),        # Unrealistic heterozygosity
      
      # Add flags for specific issues
      flag_missing = missing_rate > 0.2,
      flag_theta_var = sd_theta > 0.2,
      flag_R_var = sd_R > 1.0,
      flag_cluster = (cluster_sep < 0.1 & !is.na(cluster_sep)),
      flag_het = (het_rate > 0.9 & !is.na(het_rate))
    )
  
  return(snp_metrics)
}

create_summary_plots <- function(snp_stats) {
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
  ggsave(file.path(output_dir, "Figures", "missing_rate_dist.png"), p1, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "theta_sd_dist.png"), p2, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "cluster_sep_dist.png"), p3, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "het_rate_dist.png"), p4, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "r_sd_dist.png"), p5, width = 8, height = 6)
  
  return(list(
    missing_rate_dist = p1,
    theta_sd_dist = p2,
    cluster_sep_dist = p3,
    het_rate_dist = p4,
    r_sd_dist = p5
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

main_analysis <- function(snp_file) {
  # Process data
  snp_data <- process_snp_data(snp_file)
  
  # Analyze SNP stability
  snp_stability <- analyze_snp_stability(snp_data)
  
  # Identify problematic SNPs
  problematic_snps <- snp_stability %>%
    filter(is_problematic)
  
  # Create detailed problem summary
  problem_summary <- snp_stability %>%
    summarise(
      total_snps = n(),
      total_problematic = sum(is_problematic),
      missing_rate_issues = sum(flag_missing),
      theta_var_issues = sum(flag_theta_var),
      r_var_issues = sum(flag_R_var),
      cluster_issues = sum(flag_cluster),
      het_issues = sum(flag_het)
    )
  
  # Create summary plots
  summary_plots <- create_summary_plots(snp_stability)
  
  # Create cluster plots for problematic SNPs
  plot_snp_clusters(snp_data, snp_stability)
  
  # Save results
  write_csv(snp_stability, file.path(output_dir, "Tables", "snp_stability_metrics.csv"))
  write_csv(problematic_snps, file.path(output_dir, "Tables", "problematic_snps.csv"))
  write_csv(problem_summary, file.path(output_dir, "Tables", "problem_summary.csv"))
  
  return(list(
    snp_stability = snp_stability,
    problematic_snps = problematic_snps,
    problem_summary = problem_summary,
    summary_plots = summary_plots
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
  message("\nProblem Summary:")
  print(results$problem_summary)
}