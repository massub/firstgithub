# Required packages installation and loading
required_packages <- c(
  "tidyverse", "ggplot2", "gridExtra", 
  "corrplot", "viridis", "knitr",
  "rmarkdown"
)

# Install and load packages
for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Core Data Processing Functions
read_fdt_data <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Input file does not exist: ", file_path)
  }
  
  message("Reading raw data...")
  raw_data <- read.delim(file_path, sep = "\t", header = TRUE, check.names = FALSE)
  
  sample_cols <- grep("\\.GT$", names(raw_data), value = TRUE)
  if (length(sample_cols) == 0) {
    stop("No genotype columns found in input file")
  }
  
  sample_names <- gsub("\\.GT$", "", sample_cols)
  message(sprintf("Found %d samples", length(sample_names)))
  
  long_data_list <- list()
  
  message("Processing samples...")
  for(i in seq_along(sample_names)) {
    sample <- sample_names[i]
    message(sprintf("Processing sample %d of %d: %s", i, length(sample_names), sample))
    
    theta_col <- paste0(sample, ".Theta")
    r_col <- paste0(sample, ".R")
    gt_col <- paste0(sample, ".GT")
    
    required_cols <- c(theta_col, r_col, gt_col)
    missing_cols <- required_cols[!required_cols %in% names(raw_data)]
    
    if (length(missing_cols) > 0) {
      stop("Missing required columns for sample ", sample, ": ", 
           paste(missing_cols, collapse = ", "))
    }
    
    sample_df <- data.frame(
      Name = raw_data$Name,
      Chr = raw_data$Chr,
      Position = raw_data$Position,
      Sample = sample,
      Theta = raw_data[[theta_col]],
      R = raw_data[[r_col]],
      GT = raw_data[[gt_col]],
      stringsAsFactors = FALSE
    )
    
    long_data_list[[i]] <- sample_df
  }
  
  message("Combining sample data...")
  long_data <- do.call(rbind, long_data_list)
  
  return(list(
    long_data = long_data,
    metadata = data.frame(
      Name = raw_data$Name,
      Chr = raw_data$Chr,
      Position = raw_data$Position,
      stringsAsFactors = FALSE
    ),
    sample_count = length(sample_names)
  ))
}

calculate_snp_metrics <- function(data) {
  message("Calculating metrics for each SNP...")
  
  metrics <- data$long_data %>%
    group_by(Name) %>%
    summarise(
      # Sample counts
      total_samples = n(),
      no_calls = sum(GT == "NULL" | is.na(GT)),
      aa_calls = sum(GT == "AA", na.rm = TRUE),
      ab_calls = sum(GT == "AB", na.rm = TRUE),
      bb_calls = sum(GT == "BB", na.rm = TRUE),
      
      # Missing rate calculation (>20% threshold)
      missing_rate = no_calls / total_samples,
      
      # Theta statistics (>0.2 threshold)
      mean_theta = mean(Theta, na.rm = TRUE),
      sd_theta = sd(Theta, na.rm = TRUE),
      
      # R statistics (>1.0 threshold)
      mean_r = mean(R, na.rm = TRUE),
      sd_r = sd(R, na.rm = TRUE),
      
      # Cluster separation (<0.1 threshold)
      cluster_sep = case_when(
        sum(GT == "AA", na.rm = TRUE) > 0 & sum(GT == "BB", na.rm = TRUE) > 0 ~ {
          aa_theta <- mean(Theta[GT == "AA"], na.rm = TRUE)
          bb_theta <- mean(Theta[GT == "BB"], na.rm = TRUE)
          abs(bb_theta - aa_theta)
        },
        TRUE ~ NA_real_
      ),
      
      # Heterozygosity rate (>90% threshold)
      het_rate = ab_calls / total_samples,
      
      .groups = "drop"
    ) %>%
    left_join(data$metadata, by = "Name")
  
  return(metrics)
}

analyze_snp_quality <- function(metrics) {
  message("Analyzing SNP quality based on defined criteria...")
  
  metrics %>%
    mutate(
      # Quality flags based on specific criteria
      flag_missing = missing_rate > 0.20,  # >20% missing data
      flag_theta_var = sd_theta > 0.20,    # >0.2 theta SD
      flag_r_var = sd_r > 1.0,             # >1.0 R SD
      flag_cluster = !is.na(cluster_sep) & cluster_sep < 0.10,  # <0.1 cluster separation
      flag_het = het_rate > 0.90,          # >90% heterozygosity
      
      # Overall problematic status
      is_problematic = flag_missing | flag_theta_var | flag_r_var | flag_cluster | flag_het,
      
      # Detailed quality classification
      quality_class = case_when(
        !is_problematic ~ "Pass",
        flag_missing & !flag_theta_var & !flag_r_var & !flag_cluster & !flag_het ~ 
          sprintf("High Missing Rate (%.1f%%)", missing_rate * 100),
        !flag_missing & flag_theta_var & !flag_r_var & !flag_cluster & !flag_het ~
          sprintf("High Theta Variation (SD=%.3f)", sd_theta),
        !flag_missing & !flag_theta_var & flag_r_var & !flag_cluster & !flag_het ~
          sprintf("High R Variation (SD=%.3f)", sd_r),
        !flag_missing & !flag_theta_var & !flag_r_var & flag_cluster & !flag_het ~
          sprintf("Poor Cluster Separation (%.3f)", cluster_sep),
        !flag_missing & !flag_theta_var & !flag_r_var & !flag_cluster & flag_het ~
          sprintf("High Heterozygosity (%.1f%%)", het_rate * 100),
        TRUE ~ "Multiple Issues"
      )
    )
}

# Visualization Functions
create_quality_plots <- function(metrics, output_dir) {
  message("Creating quality metric plots...")
  plots <- list()
  
  # 1. Missing Rate Distribution
  p1 <- ggplot(metrics, aes(x = missing_rate * 100)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 20, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Missing Rate",
         subtitle = "Red line: 20% threshold",
         x = "Missing Rate (%)",
         y = "Count") +
    theme_minimal()
  plots$missing_rate <- p1
  
  # 2. Theta Variation
  p2 <- ggplot(metrics, aes(x = sd_theta)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 0.2, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Theta Standard Deviation",
         subtitle = "Red line: 0.2 threshold",
         x = "Theta SD",
         y = "Count") +
    theme_minimal()
  plots$theta_var <- p2
  
  # 3. R Variation
  p3 <- ggplot(metrics, aes(x = sd_r)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 1.0, color = "red", linetype = "dashed") +
    labs(title = "Distribution of R Standard Deviation",
         subtitle = "Red line: 1.0 threshold",
         x = "R SD",
         y = "Count") +
    theme_minimal()
  plots$r_var <- p3
  
  # 4. Cluster Separation
  p4 <- ggplot(metrics, aes(x = cluster_sep)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 0.1, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Cluster Separation",
         subtitle = "Red line: 0.1 threshold",
         x = "Cluster Separation",
         y = "Count") +
    theme_minimal()
  plots$cluster_sep <- p4
  
  # 5. Heterozygosity Rate
  p5 <- ggplot(metrics, aes(x = het_rate * 100)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "black") +
    geom_vline(xintercept = 90, color = "red", linetype = "dashed") +
    labs(title = "Distribution of Heterozygosity Rate",
         subtitle = "Red line: 90% threshold",
         x = "Heterozygosity Rate (%)",
         y = "Count") +
    theme_minimal()
  plots$het_rate <- p5
  
  # Save all plots
  ggsave(file.path(output_dir, "Figures", "missing_rate_dist.png"), p1, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "theta_var_dist.png"), p2, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "r_var_dist.png"), p3, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "cluster_sep_dist.png"), p4, width = 8, height = 6)
  ggsave(file.path(output_dir, "Figures", "het_rate_dist.png"), p5, width = 8, height = 6)
  
  return(plots)
}

# Problematic SNP Visualization
plot_problematic_examples <- function(data, metrics, output_dir) {
  message("Creating example plots for problematic SNPs...")
  
  # Get problematic SNPs for each category
  issue_categories <- c(
    "High Missing Rate", "High Theta Variation", 
    "High R Variation", "Poor Cluster Separation", 
    "High Heterozygosity", "Multiple Issues"
  )
  
  plots <- list()
  
  for(issue in issue_categories) {
    # Select SNPs for each issue type
    snps_with_issue <- metrics %>%
      filter(grepl(issue, quality_class)) %>%
      arrange(desc(missing_rate)) %>%
      slice_head(n = 2)  # Take top 2 examples per issue
    
    for(i in 1:nrow(snps_with_issue)) {
      snp <- snps_with_issue$Name[i]
      
      # Get SNP data
      snp_data <- data$long_data %>%
        filter(Name == snp) %>%
        filter(!is.na(Theta) & !is.na(R))  # Remove missing values
      
      # Get metrics for subtitle
      snp_metrics <- snps_with_issue[i, ]
      
      # Create detailed subtitle
      subtitle <- sprintf(
        "Issue: %s\nMissing Rate: %.1f%% | Theta SD: %.3f | R SD: %.3f\nCluster Sep: %.3f | Het Rate: %.1f%%",
        issue,
        snp_metrics$missing_rate * 100,
        snp_metrics$sd_theta,
        snp_metrics$sd_r,
        if(is.na(snp_metrics$cluster_sep)) 0 else snp_metrics$cluster_sep,
        snp_metrics$het_rate * 100
      )
      
      # Create plot
      p <- ggplot(snp_data, aes(x = Theta, y = R, color = GT)) +
        geom_point(alpha = 0.6) +
        scale_color_viridis_d() +
        labs(title = paste("Cluster Plot for", snp),
             subtitle = subtitle,
             x = "Theta",
             y = "R") +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10),
          legend.position = "bottom"
        )
      
      # Save plot
      plot_name <- paste0("problem_snp_", make.names(snp), "_", 
                          gsub(" ", "_", issue), ".png")
      tryCatch({
        ggsave(file.path(output_dir, "Figures", plot_name), 
               p, width = 10, height = 8)
        plots[[paste0(issue, "_", i)]] <- p
      }, error = function(e) {
        message("Warning: Could not save plot for SNP ", snp, ": ", e$message)
      })
    }
  }
  
  return(plots)
}

# Report Generation Function
generate_html_report <- function(data, metrics, output_dir) {
  message("Generating comprehensive HTML report...")
  
  # Calculate summary statistics
  summary_stats <- metrics %>%
    summarise(
      total_snps = n(),
      problematic_snps = sum(is_problematic),
      missing_rate_issues = sum(flag_missing),
      theta_var_issues = sum(flag_theta_var),
      r_var_issues = sum(flag_r_var),
      cluster_sep_issues = sum(flag_cluster),
      het_rate_issues = sum(flag_het)
    )
  
  # Create report content
  rmd_content <- c(
    "---",
    "title: SNP Quality Analysis Report",
    "author: SNP Analysis Pipeline",
    paste("date:", Sys.Date()),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: cosmo",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
    "library(knitr)",
    "library(dplyr)",
    "```",
    "",
    "# Executive Summary",
    "",
    sprintf("* Total SNPs analyzed: %d", summary_stats$total_snps),
    sprintf("* Total samples: %d", data$sample_count),
    sprintf("* Problematic SNPs identified: %d (%.1f%%)", 
            summary_stats$problematic_snps,
            summary_stats$problematic_snps/summary_stats$total_snps * 100),
    "",
    "## Quality Issues Breakdown",
    "",
    sprintf("* Missing Rate (>20%%): %d SNPs", summary_stats$missing_rate_issues),
    sprintf("* Theta Variation (SD >0.2): %d SNPs", summary_stats$theta_var_issues),
    sprintf("* R Variation (SD >1.0): %d SNPs", summary_stats$r_var_issues),
    sprintf("* Poor Cluster Separation (<0.1): %d SNPs", summary_stats$cluster_sep_issues),
    sprintf("* High Heterozygosity (>90%%): %d SNPs", summary_stats$het_rate_issues),
    "",
    "# Quality Distributions",
    "",
    "## Missing Rate Distribution",
    "![](Figures/missing_rate_dist.png)",
    "",
    "## Theta Variation Distribution",
    "![](Figures/theta_var_dist.png)",
    "",
    "## R Variation Distribution",
    "![](Figures/r_var_dist.png)",
    "",
    "## Cluster Separation Distribution",
    "![](Figures/cluster_sep_dist.png)",
    "",
    "## Heterozygosity Rate Distribution",
    "![](Figures/het_rate_dist.png)"
  )
  
  # Save and render report
  rmd_file <- file.path(output_dir, "snp_quality_report.Rmd")
  writeLines(rmd_content, rmd_file)
  
  tryCatch({
    rmarkdown::render(rmd_file, 
                      output_file = file.path(output_dir, "snp_quality_report.html"),
                      quiet = TRUE)
  }, error = function(e) {
    message("Warning: Error in report generation: ", e$message)
  })
}

# Results Export Function
export_results <- function(metrics, output_dir) {
  message("Exporting detailed results...")
  
  # Create directories
  dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
  
  # Export all SNP metrics
  write.csv(metrics, 
            file.path(output_dir, "Tables", "all_snp_metrics.csv"), 
            row.names = FALSE)
  
  # Export problematic SNPs
  problematic_snps <- metrics %>%
    filter(is_problematic) %>%
    arrange(quality_class, desc(missing_rate))
  
  write.csv(problematic_snps,
            file.path(output_dir, "Tables", "problematic_snps.csv"),
            row.names = FALSE)
  
  # Export chromosome summary if available
  if("Chr" %in% names(metrics)) {
    chr_summary <- metrics %>%
      group_by(Chr) %>%
      summarise(
        total_snps = n(),
        problematic_snps = sum(is_problematic),
        high_missing = sum(flag_missing),
        high_theta_var = sum(flag_theta_var),
        high_r_var = sum(flag_r_var),
        poor_clusters = sum(flag_cluster),
        high_het = sum(flag_het),
        .groups = 'drop'
      )
    
    write.csv(chr_summary,
              file.path(output_dir, "Tables", "chromosome_summary.csv"),
              row.names = FALSE)
  }
}

# Main Execution Function
main <- function(snp_file, output_dir, sample_limit = NULL) {
  start_time <- Sys.time()
  
  tryCatch({
    # Create output directories
    dir.create(file.path(output_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
    
    # Step 1: Read and process data
    message("\nStep 1: Reading and processing data...")
    data <- read_fdt_data(snp_file)
    
    if (!is.null(sample_limit)) {
      message(sprintf("Limiting analysis to %d samples for testing", sample_limit))
      unique_samples <- unique(data$long_data$Sample)
      selected_samples <- unique_samples[1:min(sample_limit, length(unique_samples))]
      data$long_data <- data$long_data %>% 
        filter(Sample %in% selected_samples)
      data$sample_count <- length(selected_samples)
    }
    
    # Step 2: Calculate metrics
    message("\nStep 2: Calculating metrics...")
    metrics <- calculate_snp_metrics(data)
    
    # Step 3: Analyze quality
    message("\nStep 3: Analyzing SNP quality...")
    metrics <- analyze_snp_quality(metrics)
    
    # Step 4: Create quality plots
    message("\nStep 4: Creating quality plots...")
    quality_plots <- create_quality_plots(metrics, output_dir)
    
    # Step 5: Create example plots
    message("\nStep 5: Creating example plots for problematic SNPs...")
    problem_plots <- plot_problematic_examples(data, metrics, output_dir)
    
    # Step 6: Generate report
    message("\nStep 6: Generating comprehensive report...")
    generate_html_report(data, metrics, output_dir)
    
    # Step 7: Export results
    message("\nStep 7: Exporting detailed results...")
    export_results(metrics, output_dir)
    
    # Calculate execution time
    end_time <- Sys.time()
    execution_time <- difftime(end_time, start_time, units = "mins")
    
    message(sprintf("\nAnalysis completed in %.1f minutes", as.numeric(execution_time)))
    
    return(list(
      metrics = metrics,
      execution_time = execution_time
    ))
    
  }, error = function(e) {
    message("\nERROR: Analysis failed with the following error:")
    message(e$message)
    return(structure(e$message, class = "try-error"))
  })
}

# Set file paths and run analysis
snp_file <- "D:/USDA/Small Grains Research Unit/Projects/Project--SNP position and cluster charecterization/Project_Root/00_Raw_Data/Oat_Datasets/SDSU2023Caffe-O3K_01-08-FDT.txt"
output_dir <- "D:/USDA/Small Grains Research Unit/Projects/Project--SNP position and cluster charecterization/Project_Root/03_Results"

# Run the analysis
results <- main(snp_file, output_dir)

# Check results
if(inherits(results, "try-error")) {
  message("Analysis failed. Please check the error messages above.")
} else {
  message("Analysis completed successfully!")
}


############# PART 2 ##############
##################################
##################################

# Function for analyzing problematic theta values
analyze_theta_problems <- function(snp_file, metrics, output_dir) {
  # Initialize
  if (!require("moments")) install.packages("moments")
  if (!require("dplyr")) install.packages("dplyr")
  library(moments)
  library(dplyr)
  
  message("Starting analysis...")
  theta_dir <- file.path(output_dir, "Theta_Analysis")
  dir.create(theta_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(theta_dir, "Problem_Examples"), recursive = TRUE, showWarnings = FALSE)
  
  # Load data
  raw_data <- read.delim(snp_file, sep="\t", header=TRUE)
  theta_problems <- metrics %>%
    dplyr::filter(flag_theta_var) %>%
    arrange(desc(sd_theta))
  
  # Get sample data
  sample_cols <- grep("\\.GT$", names(raw_data), value = TRUE)
  sample_names <- gsub("\\.GT$", "", sample_cols)
  
  # Process data
  message("Processing samples...")
  problem_data <- do.call(rbind, lapply(sample_names, function(sample) {
    data.frame(
      Name = raw_data$Name,
      GT = raw_data[[paste0(sample, ".GT")]],
      Theta = raw_data[[paste0(sample, ".Theta")]],
      R = raw_data[[paste0(sample, ".R")]],
      Sample = sample,
      stringsAsFactors = FALSE
    )
  })) %>%
    dplyr::filter(Name %in% theta_problems$Name)
  
  # Calculate metrics
  message("Analyzing theta distributions...")
  detailed_theta <- problem_data %>%
    group_by(Name) %>%
    summarise(
      # Basic stats
      total_samples = n(),
      valid_samples = sum(!is.na(Theta)),
      mean_theta = mean(Theta, na.rm = TRUE),
      median_theta = median(Theta, na.rm = TRUE),
      sd_theta = sd(Theta, na.rm = TRUE),
      
      # Genotype counts
      aa_count = sum(GT == "AA", na.rm = TRUE),
      ab_count = sum(GT == "AB", na.rm = TRUE),
      bb_count = sum(GT == "BB", na.rm = TRUE),
      
      # Cluster stats
      aa_theta_mean = if(sum(GT == "AA", na.rm = TRUE) >= 3) mean(Theta[GT == "AA"], na.rm = TRUE) else NA_real_,
      ab_theta_mean = if(sum(GT == "AB", na.rm = TRUE) >= 3) mean(Theta[GT == "AB"], na.rm = TRUE) else NA_real_,
      bb_theta_mean = if(sum(GT == "BB", na.rm = TRUE) >= 3) mean(Theta[GT == "BB"], na.rm = TRUE) else NA_real_,
      
      aa_theta_sd = if(sum(GT == "AA", na.rm = TRUE) >= 3) sd(Theta[GT == "AA"], na.rm = TRUE) else NA_real_,
      ab_theta_sd = if(sum(GT == "AB", na.rm = TRUE) >= 3) sd(Theta[GT == "AB"], na.rm = TRUE) else NA_real_,
      bb_theta_sd = if(sum(GT == "BB", na.rm = TRUE) >= 3) sd(Theta[GT == "BB"], na.rm = TRUE) else NA_real_,
      
      # Distribution shape
      skewness = if(sum(!is.na(Theta)) >= 3) moments::skewness(Theta, na.rm = TRUE) else NA_real_,
      kurtosis = if(sum(!is.na(Theta)) >= 4) moments::kurtosis(Theta, na.rm = TRUE) else NA_real_
    ) %>%
    # Join with original metrics data
    mutate(
      Chr = metrics$Chr[match(Name, metrics$Name)],
      Position = metrics$Position[match(Name, metrics$Name)]
    ) %>%
    mutate(
      problem_type = case_when(
        !is.na(ab_theta_sd) & ab_theta_sd > 0.1 & ab_count >= 5 ~ "High Heterozygous Variance",
        (!is.na(aa_theta_sd) & aa_theta_sd > 0.1 & aa_count >= 5) |
          (!is.na(bb_theta_sd) & bb_theta_sd > 0.1 & bb_count >= 5) ~ "High Homozygous Variance",
        !is.na(skewness) & abs(skewness) > 1 ~ "Skewed Distribution",
        TRUE ~ "Other Theta Issues"
      )
    )
  
  # Visualizations
  message("Creating plots...")
  
  # Chromosome distribution
  p1 <- ggplot(detailed_theta, aes(x = as.factor(Chr), y = sd_theta)) +
    geom_boxplot(fill = "steelblue", alpha = 0.7) +
    labs(title = "Distribution of Theta SD by Chromosome",
         x = "Chromosome",
         y = "Theta Standard Deviation") +
    theme_minimal()
  ggsave(file.path(theta_dir, "theta_sd_by_chromosome.png"), p1, width = 10, height = 6)
  
  # Problem types
  p2 <- ggplot(detailed_theta, aes(x = problem_type)) +
    geom_bar(fill = "steelblue") +
    labs(title = "Distribution of Theta Problem Types",
         x = "Problem Type",
         y = "Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(theta_dir, "theta_problem_types.png"), p2, width = 10, height = 6)
  
  # Example plots
  message("Generating example plots...")
  for(problem in unique(detailed_theta$problem_type)) {
    examples <- detailed_theta %>%
      dplyr::filter(problem_type == problem) %>%
      arrange(desc(sd_theta)) %>%
      slice_head(n = 3)
    
    for(i in 1:nrow(examples)) {
      snp_data <- problem_data %>%
        dplyr::filter(Name == examples$Name[i])
      
      p <- ggplot(snp_data) +
        geom_density(aes(x = Theta, fill = GT), alpha = 0.5) +
        geom_point(aes(x = Theta, y = -0.1, color = GT), 
                   position = position_jitter(height = 0.05), alpha = 0.5) +
        scale_fill_viridis_d() +
        scale_color_viridis_d() +
        labs(title = paste("Theta Distribution for", examples$Name[i]),
             subtitle = sprintf("Problem Type: %s\nSD: %.3f, Skewness: %.3f",
                                problem,
                                examples$sd_theta[i],
                                examples$skewness[i]),
             y = "Density") +
        theme_minimal() +
        ylim(-0.2, NA)
      
      ggsave(file.path(theta_dir, "Problem_Examples",
                       paste0("theta_problem_", make.names(examples$Name[i]), ".png")),
             p, width = 8, height = 6)
    }
  }
  
  # Export results
  message("Saving results...")
  write.csv(detailed_theta,
            file.path(theta_dir, "detailed_theta_analysis.csv"),
            row.names = FALSE)
  
  # Summary
  summary_stats <- detailed_theta %>%
    summarise(
      total_problems = n(),
      mean_sd = mean(sd_theta, na.rm = TRUE),
      max_sd = max(sd_theta, na.rm = TRUE),
      high_het_var = sum(problem_type == "High Heterozygous Variance", na.rm = TRUE),
      high_homo_var = sum(problem_type == "High Homozygous Variance", na.rm = TRUE),
      skewed = sum(problem_type == "Skewed Distribution", na.rm = TRUE),
      other = sum(problem_type == "Other Theta Issues", na.rm = TRUE)
    )
  
  write.csv(summary_stats,
            file.path(theta_dir, "theta_summary_stats.csv"),
            row.names = FALSE)
  
  message("\nSummary of Theta Problems:")
  print(summary_stats)
  
  message("\nAnalysis complete! Results are in the 'Theta_Analysis' folder.")
  
  return(list(
    detailed_metrics = detailed_theta,
    summary_stats = summary_stats,
    raw_problem_data = problem_data
  ))
}


metrics <- results$metrics
theta_analysis <- analyze_theta_problems(snp_file, metrics, output_dir)

# Function to create comparative visualizations
create_comparison_plots <- function(theta_analysis, raw_data, output_dir) {
  # Load required packages
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  library(viridis)
  
  message("Creating comparative visualizations...")
  
  # Create directory for comparison plots
  comp_dir <- file.path(output_dir, "Theta_Analysis", "Comparisons")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Get examples of good SNPs (low variance in all clusters)
  good_snps <- theta_analysis$detailed_metrics %>%
    dplyr::filter(!is.na(aa_theta_sd) & !is.na(ab_theta_sd) & !is.na(bb_theta_sd)) %>%
    dplyr::filter(aa_theta_sd < 0.05, ab_theta_sd < 0.05, bb_theta_sd < 0.05) %>%
    dplyr::filter(aa_count >= 10, ab_count >= 10, bb_count >= 10) %>%
    arrange(sd_theta) %>%
    head(3)
  
  # Get worst examples of each problem type
  het_problems <- theta_analysis$detailed_metrics %>%
    dplyr::filter(problem_type == "High Heterozygous Variance") %>%
    arrange(desc(ab_theta_sd)) %>%
    head(3)
  
  homo_problems <- theta_analysis$detailed_metrics %>%
    dplyr::filter(problem_type == "High Homozygous Variance") %>%
    arrange(desc(pmax(aa_theta_sd, bb_theta_sd, na.rm=TRUE))) %>%
    head(3)
  
  skewed_problems <- theta_analysis$detailed_metrics %>%
    dplyr::filter(problem_type == "Skewed Distribution") %>%
    arrange(desc(abs(skewness))) %>%
    head(3)
  
  # Function to create detailed SNP plot
  create_detailed_plot <- function(snp_data, title, metrics, snp_name) {
    # Create main plot with both density and points
    p <- ggplot(snp_data) +
      # Density plots for each genotype
      geom_density(aes(x = Theta, fill = GT), alpha = 0.4) +
      # Add jittered points below
      geom_point(aes(x = Theta, y = -0.1, color = GT), 
                 position = position_jitter(height = 0.05), alpha = 0.5) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() +
      labs(title = title,
           subtitle = sprintf("SNP: %s\nSD: AA=%.3f, AB=%.3f, BB=%.3f\nCounts: AA=%d, AB=%d, BB=%d\nSkewness: %.3f",
                              snp_name,
                              metrics$aa_theta_sd, metrics$ab_theta_sd, metrics$bb_theta_sd,
                              metrics$aa_count, metrics$ab_count, metrics$bb_count,
                              metrics$skewness),
           x = "Theta",
           y = "Density") +
      theme_minimal() +
      ylim(-0.2, NA)
    
    return(p)
  }
  
  # 1. Create Good vs Bad Heterozygous Comparison
  message("Creating heterozygous comparison plots...")
  good_het <- good_snps[1,]
  bad_het <- het_problems[1,]
  
  p1 <- create_detailed_plot(
    theta_analysis$raw_problem_data %>% dplyr::filter(Name == good_het$Name),
    "Good SNP: Well-clustered Heterozygous Calls",
    good_het,
    good_het$Name
  )
  
  p2 <- create_detailed_plot(
    theta_analysis$raw_problem_data %>% dplyr::filter(Name == bad_het$Name),
    "Problem SNP: High Heterozygous Variance",
    bad_het,
    bad_het$Name
  )
  
  # Save combined plot directly
  png(file.path(comp_dir, "heterozygous_comparison.png"), 
      width = 15, height = 7, units = "in", res = 300)
  grid.arrange(p1, p2, ncol=2)
  dev.off()
  
  # 2. Create Good vs Bad Homozygous Comparison
  message("Creating homozygous comparison plots...")
  good_homo <- good_snps[2,]
  bad_homo <- homo_problems[1,]
  
  p3 <- create_detailed_plot(
    theta_analysis$raw_problem_data %>% dplyr::filter(Name == good_homo$Name),
    "Good SNP: Well-clustered Homozygous Calls",
    good_homo,
    good_homo$Name
  )
  
  p4 <- create_detailed_plot(
    theta_analysis$raw_problem_data %>% dplyr::filter(Name == bad_homo$Name),
    "Problem SNP: High Homozygous Variance",
    bad_homo,
    bad_homo$Name
  )
  
  # Save combined plot directly
  png(file.path(comp_dir, "homozygous_comparison.png"), 
      width = 15, height = 7, units = "in", res = 300)
  grid.arrange(p3, p4, ncol=2)
  dev.off()
  
  # 3. Create Normal vs Skewed Distribution Comparison
  message("Creating skewness comparison plots...")
  good_dist <- good_snps[3,]
  bad_dist <- skewed_problems[1,]
  
  p5 <- create_detailed_plot(
    theta_analysis$raw_problem_data %>% dplyr::filter(Name == good_dist$Name),
    "Good SNP: Normal Distribution",
    good_dist,
    good_dist$Name
  )
  
  p6 <- create_detailed_plot(
    theta_analysis$raw_problem_data %>% dplyr::filter(Name == bad_dist$Name),
    "Problem SNP: Skewed Distribution",
    bad_dist,
    bad_dist$Name
  )
  
  # Save combined plot directly
  png(file.path(comp_dir, "distribution_comparison.png"), 
      width = 15, height = 7, units = "in", res = 300)
  grid.arrange(p5, p6, ncol=2)
  dev.off()
  
  # Create summary table of examples with more details
  example_summary <- bind_rows(
    mutate(good_snps, category = "Good SNPs"),
    mutate(het_problems, category = "High Het Variance"),
    mutate(homo_problems, category = "High Homo Variance"),
    mutate(skewed_problems, category = "Skewed Distribution")
  ) %>%
    dplyr::select(category, Name, Chr, Position, sd_theta, aa_theta_sd, ab_theta_sd, bb_theta_sd,
                  aa_count, ab_count, bb_count, skewness)
  
  # Export summary
  write.csv(example_summary,
            file.path(comp_dir, "example_snps_summary.csv"),
            row.names = FALSE)
  
  message("Comparison visualizations complete! Check the 'Comparisons' folder.")
  
  return(example_summary)
}

comparison_results <- create_comparison_plots(theta_analysis, raw_data, output_dir)


############# PART 3 ##############
##################################
##################################

# Function for advanced theta analysis
analyze_theta_patterns <- function(theta_analysis, metrics, output_dir) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  
  message("Performing advanced theta analysis...")
  
  adv_dir <- file.path(output_dir, "Theta_Analysis", "Advanced_Analysis")
  dir.create(adv_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Chromosome-level analysis
  chr_analysis <- theta_analysis$detailed_metrics %>%
    group_by(Chr) %>%
    summarise(
      total_snps = n(),
      high_het_var = sum(problem_type == "High Heterozygous Variance"),
      high_homo_var = sum(problem_type == "High Homozygous Variance"),
      skewed = sum(problem_type == "Skewed Distribution"),
      mean_sd_theta = mean(sd_theta, na.rm = TRUE),
      median_sd_theta = median(sd_theta, na.rm = TRUE),
      problem_rate = total_snps / sum(metrics$flag_theta_var) * 100
    ) %>%
    arrange(desc(problem_rate))
  
  # 2. Position effect analysis
  position_analysis <- theta_analysis$detailed_metrics %>%
    mutate(
      position_mb = Position / 1000000,  # Convert to Mb
      chr_position = paste(Chr, position_mb)
    ) %>%
    arrange(Chr, Position)
  
  # Create chromosome map of problems
  p1 <- ggplot(position_analysis, 
               aes(x = position_mb, y = as.factor(Chr), color = sd_theta)) +
    geom_point(alpha = 0.6) +
    scale_color_viridis() +
    labs(title = "Chromosome-wide Distribution of Theta SD",
         x = "Position (Mb)",
         y = "Chromosome",
         color = "Theta SD") +
    theme_minimal()
  
  # 3. Cluster pattern analysis
  cluster_patterns <- theta_analysis$detailed_metrics %>%
    mutate(
      # Calculate cluster ratios
      ab_ratio = ab_count / (aa_count + ab_count + bb_count),
      aa_bb_ratio = aa_count / (aa_count + bb_count),
      # Identify potential underlying issues
      potential_issue = case_when(
        ab_ratio > 0.5 ~ "Excess Heterozygotes",
        ab_ratio < 0.1 ~ "Missing Heterozygotes",
        aa_bb_ratio > 0.8 ~ "AA Bias",
        aa_bb_ratio < 0.2 ~ "BB Bias",
        TRUE ~ "Normal Distribution"
      )
    )
  
  # 4. Create detailed summary report
  problem_summary <- data.frame(
    # Overall statistics
    total_analyzed = nrow(theta_analysis$detailed_metrics),
    high_het_problems = sum(cluster_patterns$problem_type == "High Heterozygous Variance"),
    high_homo_problems = sum(cluster_patterns$problem_type == "High Homozygous Variance"),
    skewed_problems = sum(cluster_patterns$problem_type == "Skewed Distribution"),
    
    # Cluster patterns
    excess_het = sum(cluster_patterns$potential_issue == "Excess Heterozygotes"),
    missing_het = sum(cluster_patterns$potential_issue == "Missing Heterozygotes"),
    aa_bias = sum(cluster_patterns$potential_issue == "AA Bias"),
    bb_bias = sum(cluster_patterns$potential_issue == "BB Bias"),
    
    # Average metrics
    mean_theta_sd = mean(theta_analysis$detailed_metrics$sd_theta, na.rm = TRUE),
    median_theta_sd = median(theta_analysis$detailed_metrics$sd_theta, na.rm = TRUE),
    worst_theta_sd = max(theta_analysis$detailed_metrics$sd_theta, na.rm = TRUE)
  )
  
  # 5. Identify potential systematic issues
  systematic_issues <- cluster_patterns %>%
    group_by(Chr) %>%
    summarise(
      total_problems = n(),
      excess_het_rate = mean(potential_issue == "Excess Heterozygotes"),
      missing_het_rate = mean(potential_issue == "Missing Heterozygotes"),
      aa_bias_rate = mean(potential_issue == "AA Bias"),
      bb_bias_rate = mean(potential_issue == "BB Bias"),
      mean_ab_ratio = mean(ab_ratio, na.rm = TRUE),
      mean_aa_bb_ratio = mean(aa_bb_ratio, na.rm = TRUE)
    ) %>%
    mutate(
      potential_systematic_issue = case_when(
        excess_het_rate > 0.4 ~ "Possible paralog issues",
        missing_het_rate > 0.4 ~ "Possible null allele issues",
        aa_bias_rate > 0.4 ~ "Possible probe bias towards reference",
        bb_bias_rate > 0.4 ~ "Possible probe bias towards alternate",
        TRUE ~ "No systematic issues detected"
      )
    )
  
  # Save all results
  write.csv(chr_analysis, 
            file.path(adv_dir, "chromosome_level_analysis.csv"), 
            row.names = FALSE)
  
  write.csv(cluster_patterns, 
            file.path(adv_dir, "cluster_pattern_analysis.csv"), 
            row.names = FALSE)
  
  write.csv(systematic_issues, 
            file.path(adv_dir, "systematic_issues_analysis.csv"), 
            row.names = FALSE)
  
  write.csv(problem_summary, 
            file.path(adv_dir, "problem_summary.csv"), 
            row.names = FALSE)
  
  # Save chromosome map plot
  png(file.path(adv_dir, "chromosome_theta_map.png"), 
      width = 12, height = 8, units = "in", res = 300)
  print(p1)
  dev.off()
  
  message("Advanced analysis complete! Check the 'Advanced_Analysis' folder.")
  
  return(list(
    chr_analysis = chr_analysis,
    cluster_patterns = cluster_patterns,
    systematic_issues = systematic_issues,
    problem_summary = problem_summary
  ))
}

advanced_results <- analyze_theta_patterns(theta_analysis, metrics, output_dir)

analyze_null_alleles <- function(theta_analysis, metrics, output_dir) {
  library(dplyr)
  library(ggplot2)
  
  message("Analyzing null allele patterns and problem rates...")
  
  null_dir <- file.path(output_dir, "Theta_Analysis", "Null_Allele_Analysis")
  dir.create(null_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate expected Hardy-Weinberg proportions
  null_patterns <- theta_analysis$detailed_metrics %>%
    mutate(
      total_calls = aa_count + ab_count + bb_count,
      # Calculate observed versus expected heterozygosity
      obs_het_freq = ab_count / total_calls,
      exp_het_freq = 2 * (aa_count/total_calls) * (bb_count/total_calls),
      het_deficit = exp_het_freq - obs_het_freq,
      
      # Classify severity of heterozygote deficit
      null_allele_status = case_when(
        het_deficit > 0.3 ~ "Severe heterozygote deficit",
        het_deficit > 0.2 ~ "Moderate heterozygote deficit",
        het_deficit > 0.1 ~ "Mild heterozygote deficit",
        TRUE ~ "Normal"
      ),
      
      # Calculate homozygote ratio
      homo_ratio = (aa_count + bb_count) / total_calls,
      
      # Additional flags
      high_missing = (total_samples - total_calls)/total_samples > 0.1,
      excess_homo = homo_ratio > 0.9
    ) %>%
    arrange(desc(het_deficit))
  
  # Chromosome-wise problem rate analysis
  chr_problems <- null_patterns %>%
    group_by(Chr) %>%
    summarise(
      total_snps = n(),
      null_allele_candidates = sum(null_allele_status != "Normal"),
      severe_nulls = sum(null_allele_status == "Severe heterozygote deficit"),
      problem_rate = (null_allele_candidates / n()) * 100,
      mean_het_deficit = mean(het_deficit, na.rm = TRUE),
      mean_homo_ratio = mean(homo_ratio, na.rm = TRUE),
      high_missing_rate = sum(high_missing),
      .groups = 'drop'
    ) %>%
    arrange(desc(problem_rate))
  
  # Create visualization of heterozygote deficit by chromosome
  p1 <- ggplot(null_patterns, aes(x = as.factor(Chr), y = het_deficit)) +
    geom_boxplot(fill = "steelblue", alpha = 0.6) +
    geom_hline(yintercept = c(0.1, 0.2, 0.3), 
               linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "Heterozygote Deficit by Chromosome",
         subtitle = "Dashed lines indicate mild, moderate, and severe deficit thresholds",
         x = "Chromosome",
         y = "Heterozygote Deficit") +
    theme_minimal()
  
  # Plot relationship between heterozygote deficit and theta SD
  p2 <- ggplot(null_patterns, 
               aes(x = het_deficit, y = sd_theta, color = null_allele_status)) +
    geom_point(alpha = 0.6) +
    scale_color_viridis_d() +
    labs(title = "Relationship between Heterozygote Deficit and Theta SD",
         x = "Heterozygote Deficit",
         y = "Theta SD") +
    theme_minimal()
  
  # Save plots
  png(file.path(null_dir, "het_deficit_by_chr.png"), 
      width = 12, height = 8, units = "in", res = 300)
  print(p1)
  dev.off()
  
  png(file.path(null_dir, "het_deficit_vs_theta_sd.png"), 
      width = 12, height = 8, units = "in", res = 300)
  print(p2)
  dev.off()
  
  # Create summary report
  summary_stats <- data.frame(
    total_analyzed = nrow(null_patterns),
    severe_nulls = sum(null_patterns$null_allele_status == "Severe heterozygote deficit"),
    moderate_nulls = sum(null_patterns$null_allele_status == "Moderate heterozygote deficit"),
    mild_nulls = sum(null_patterns$null_allele_status == "Mild heterozygote deficit"),
    high_missing_rate = sum(null_patterns$high_missing),
    excess_homo = sum(null_patterns$excess_homo)
  )
  
  # Export results
  write.csv(null_patterns, 
            file.path(null_dir, "detailed_null_allele_analysis.csv"), 
            row.names = FALSE)
  
  write.csv(chr_problems, 
            file.path(null_dir, "chromosome_null_patterns.csv"), 
            row.names = FALSE)
  
  write.csv(summary_stats, 
            file.path(null_dir, "null_allele_summary.csv"), 
            row.names = FALSE)
  
  message("Null allele analysis complete! Check the 'Null_Allele_Analysis' folder.")
  
  return(list(
    null_patterns = null_patterns,
    chr_problems = chr_problems,
    summary_stats = summary_stats
  ))
}

null_allele_results <- analyze_null_alleles(theta_analysis, metrics, output_dir)

############# PART 4 ##############
##################################
##################################

# Function to analyze homozygous cluster quality
analyze_homozygous_clusters <- function(theta_analysis, metrics, output_dir) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  
  message("Analyzing homozygous cluster quality...")
  
  homo_dir <- file.path(output_dir, "Theta_Analysis", "Homozygous_Clusters")
  dir.create(homo_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Calculate cluster quality metrics
  cluster_quality <- theta_analysis$detailed_metrics %>%
    mutate(
      # Homozygous call rates
      total_calls = aa_count + ab_count + bb_count,
      homo_call_rate = (aa_count + bb_count) / total_calls,
      
      # Cluster separation quality - handle NAs
      aa_bb_separation = if_else(is.na(aa_theta_mean) | is.na(bb_theta_mean),
                                 NA_real_,
                                 abs(aa_theta_mean - bb_theta_mean)),
      
      # Cluster tightness
      aa_cluster_quality = if_else(aa_count >= 5, aa_theta_sd, NA_real_),
      bb_cluster_quality = if_else(bb_count >= 5, bb_theta_sd, NA_real_),
      
      # Quality classifications
      cluster_quality = case_when(
        is.na(aa_bb_separation) ~ "Insufficient Data",
        aa_bb_separation < 0.5 ~ "Poor Separation",
        (aa_cluster_quality > 0.1 | bb_cluster_quality > 0.1) ~ "Diffuse Clusters",
        (aa_theta_mean > 0.1 | bb_theta_mean < 0.9) ~ "Shifted Clusters",
        TRUE ~ "Good Quality"
      ),
      
      # Signal strength (using R values if available)
      signal_strength = case_when(
        cluster_quality == "Good Quality" & homo_call_rate > 0.9 ~ "Strong",
        cluster_quality == "Good Quality" ~ "Moderate",
        TRUE ~ "Weak"
      ),
      
      # Add recommendation
      recommendation = case_when(
        cluster_quality == "Good Quality" & signal_strength == "Strong" ~ 
          "High confidence - Keep",
        cluster_quality == "Good Quality" & signal_strength == "Moderate" ~
          "Moderate confidence - Review if critical",
        cluster_quality == "Poor Separation" ~ 
          "Consider redesign - Poor cluster separation",
        cluster_quality == "Diffuse Clusters" ~
          "Consider redesign - High cluster variance",
        cluster_quality == "Shifted Clusters" ~
          "Review probe sequence - Possible binding issues",
        TRUE ~ "Review manually"
      )
    )
  
  # Create summary visualizations
  # 1. Cluster separation plot
  p1 <- ggplot(cluster_quality %>% filter(!is.na(aa_bb_separation)), 
               aes(x = aa_bb_separation)) +
    geom_histogram(bins = 50, fill = "steelblue") +
    geom_vline(xintercept = 0.5, color = "red", linetype = "dashed") +
    labs(title = "Distribution of AA/BB Cluster Separation",
         subtitle = "Red line indicates minimum recommended separation (0.5)",
         x = "AA to BB Cluster Separation",
         y = "Count") +
    theme_minimal()
  
  # 2. Cluster quality by chromosome
  p2 <- ggplot(cluster_quality, 
               aes(x = as.factor(Chr), fill = cluster_quality)) +
    geom_bar(position = "fill") +
    scale_fill_viridis_d() +
    labs(title = "SNP Cluster Quality by Chromosome",
         x = "Chromosome",
         y = "Proportion",
         fill = "Cluster Quality") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # 3. Cluster mean positions
  p3 <- ggplot(cluster_quality %>% 
                 filter(!is.na(aa_theta_mean) & !is.na(bb_theta_mean))) +
    geom_point(aes(x = aa_theta_mean, y = bb_theta_mean, 
                   color = cluster_quality), alpha = 0.6) +
    scale_color_viridis_d() +
    labs(title = "AA vs BB Cluster Mean Positions",
         x = "AA Cluster Mean",
         y = "BB Cluster Mean") +
    theme_minimal()
  
  # Save plots
  png(file.path(homo_dir, "cluster_separation_dist.png"), 
      width = 10, height = 6, units = "in", res = 300)
  print(p1)
  dev.off()
  
  png(file.path(homo_dir, "cluster_quality_by_chr.png"), 
      width = 12, height = 7, units = "in", res = 300)
  print(p2)
  dev.off()
  
  png(file.path(homo_dir, "cluster_positions.png"), 
      width = 10, height = 8, units = "in", res = 300)
  print(p3)
  dev.off()
  
  # Create detailed quality summary
  quality_summary <- cluster_quality %>%
    group_by(Chr) %>%
    summarise(
      total_snps = n(),
      poor_separation = sum(cluster_quality == "Poor Separation", na.rm = TRUE),
      diffuse_clusters = sum(cluster_quality == "Diffuse Clusters", na.rm = TRUE),
      shifted_clusters = sum(cluster_quality == "Shifted Clusters", na.rm = TRUE),
      good_quality = sum(cluster_quality == "Good Quality", na.rm = TRUE),
      insufficient_data = sum(cluster_quality == "Insufficient Data", na.rm = TRUE),
      
      # Averages
      mean_separation = mean(aa_bb_separation, na.rm = TRUE),
      mean_aa_sd = mean(aa_cluster_quality, na.rm = TRUE),
      mean_bb_sd = mean(bb_cluster_quality, na.rm = TRUE),
      
      # Signal quality
      strong_signals = sum(signal_strength == "Strong", na.rm = TRUE),
      moderate_signals = sum(signal_strength == "Moderate", na.rm = TRUE),
      weak_signals = sum(signal_strength == "Weak", na.rm = TRUE),
      
      # Success rate
      success_rate = (good_quality / total_snps) * 100
    ) %>%
    arrange(desc(success_rate))
  
  # Create SNP recommendations
  snp_recommendations <- cluster_quality %>%
    dplyr::select(Name, Chr, Position, cluster_quality, signal_strength, 
                  aa_bb_separation, aa_cluster_quality, bb_cluster_quality,
                  homo_call_rate, recommendation)
  
  # Export results
  write.csv(quality_summary,
            file.path(homo_dir, "chromosome_quality_summary.csv"),
            row.names = FALSE)
  
  write.csv(snp_recommendations,
            file.path(homo_dir, "snp_recommendations.csv"),
            row.names = FALSE)
  
  # Print summary
  message("\nSNP Quality Summary:")
  print(table(cluster_quality$cluster_quality))
  
  message("\nSignal Strength Summary:")
  print(table(cluster_quality$signal_strength))
  
  message("\nAnalysis complete! Check the 'Homozygous_Clusters' folder.")
  
  return(list(
    cluster_quality = cluster_quality,
    quality_summary = quality_summary,
    recommendations = snp_recommendations
  ))
}

homo_results <- analyze_homozygous_clusters(theta_analysis, metrics, output_dir)

############# PART 5 ##############
##################################
##################################

# Function to create detailed SNP comparisons with separate figures
analyze_snp_examples <- function(theta_analysis, homo_results, output_dir) {
  library(dplyr)
  library(ggplot2)
  library(gridExtra)
  
  message("Creating detailed SNP comparison analysis...")
  
  # Create directory for detailed examples within Homozygous_Clusters
  example_dir <- file.path(output_dir, "Theta_Analysis", "Homozygous_Clusters", "Detailed_Examples")
  dir.create(example_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Find interesting example cases
  examples <- homo_results$cluster_quality %>%
    mutate(
      case_type = case_when(
        # Case 1: Failed first filter but good clusters
        sd_theta > 0.2 & 
          aa_bb_separation > 0.7 & 
          (!is.na(aa_cluster_quality) & aa_cluster_quality < 0.05) &
          (!is.na(bb_cluster_quality) & bb_cluster_quality < 0.05) ~ 
          "High_SD_but_Good_Clusters",
        
        # Case 2: Failed both filters
        sd_theta > 0.2 & 
          (aa_bb_separation < 0.5 | 
             ((!is.na(aa_cluster_quality) & aa_cluster_quality > 0.1) |
                (!is.na(bb_cluster_quality) & bb_cluster_quality > 0.1))) ~
          "Failed_Both_Analyses",
        
        # Case 3: Borderline cases
        sd_theta > 0.18 & sd_theta < 0.22 &
          aa_bb_separation > 0.6 & aa_bb_separation < 0.8 ~
          "Borderline_Case",
        
        TRUE ~ "Other"
      )
    )
  
  # Filter for examples
  examples <- examples %>%
    filter(case_type != "Other") %>%
    group_by(case_type) %>%
    slice_head(n = 2) %>%  # Get top 2 examples of each type
    ungroup()
  
  # Create summary of what we found
  message("\nExample cases found:")
  print(table(examples$case_type))
  
  # Function to create detailed plot for each SNP
  create_detailed_plot <- function(snp_data, metrics, title) {
    p <- ggplot(snp_data) +
      # Density plot
      geom_density(aes(x = Theta, fill = GT), alpha = 0.4) +
      # Add points
      geom_point(aes(x = Theta, y = -0.1, color = GT), 
                 position = position_jitter(height = 0.05), alpha = 0.5) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() +
      labs(title = title,
           subtitle = sprintf(
             "Overall SD: %.3f\nAA-BB Separation: %.3f\nCluster SDs: AA=%.3f, BB=%.3f\nCounts: AA=%d, AB=%d, BB=%d",
             metrics$sd_theta, 
             if(is.na(metrics$aa_bb_separation)) 0 else metrics$aa_bb_separation,
             if(is.na(metrics$aa_cluster_quality)) 0 else metrics$aa_cluster_quality,
             if(is.na(metrics$bb_cluster_quality)) 0 else metrics$bb_cluster_quality,
             metrics$aa_count, metrics$ab_count, metrics$bb_count
           )) +
      theme_minimal() +
      ylim(-0.2, NA)
    return(p)
  }
  
  # Create and save individual plots
  for(i in 1:nrow(examples)) {
    snp_data <- theta_analysis$raw_problem_data %>%
      filter(Name == examples$Name[i])
    
    p <- create_detailed_plot(
      snp_data,
      examples[i,],
      sprintf("SNP: %s\n%s", 
              examples$Name[i],
              gsub("_", " ", examples$case_type[i]))
    )
    
    # Save individual plot
    png(file.path(example_dir, 
                  sprintf("%s_SNP_%s.png", 
                          examples$case_type[i],
                          make.names(examples$Name[i]))),
        width = 10, height = 8, units = "in", res = 300)
    print(p)
    dev.off()
  }
  
  # Also create a combined comparison plot
  all_plots <- list()
  for(i in 1:nrow(examples)) {
    snp_data <- theta_analysis$raw_problem_data %>%
      filter(Name == examples$Name[i])
    
    all_plots[[i]] <- create_detailed_plot(
      snp_data,
      examples[i,],
      sprintf("SNP: %s\n%s", 
              examples$Name[i],
              gsub("_", " ", examples$case_type[i]))
    )
  }
  
  # Save combined plot
  png(file.path(example_dir, "all_examples_comparison.png"),
      width = 15, height = 12, units = "in", res = 300)
  do.call(grid.arrange, c(all_plots, ncol=2))
  dev.off()
  
  # Create detailed summary table
  example_summary <- examples %>%
    dplyr::select(dplyr::any_of(c("Name", "case_type", "sd_theta", "aa_bb_separation",
                                  "aa_cluster_quality", "bb_cluster_quality",
                                  "aa_count", "ab_count", "bb_count",
                                  "cluster_quality", "recommendation")))
  
  write.csv(example_summary,
            file.path(example_dir, "detailed_examples.csv"),
            row.names = FALSE)
  
  message("\nExample analysis complete! Check the 'Homozygous_Clusters/Detailed_Examples' folder.")
  message("Individual plots and a combined comparison plot have been generated.")
  
  return(list(
    examples = example_summary,
    n_examples = nrow(examples)
  ))
}

example_results <- analyze_snp_examples(theta_analysis, homo_results, output_dir)


############# PART 6 ##############
##################################
##################################

analyze_theta_patterns_advanced <- function(theta_analysis, homo_results, metrics, output_dir) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(viridis)
  library(reshape2)
  
  message("Performing advanced theta pattern analysis...")
  
  # Create directory for advanced analysis
  adv_dir <- file.path(output_dir, "Theta_Analysis", "Advanced_Insights")
  dir.create(adv_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Function to create a cleaner correlation matrix
  clean_cor_matrix <- function(data) {
    # Get numeric columns
    numeric_cols <- data %>%
      dplyr::select(where(is.numeric)) %>%
      dplyr::select(-Position) # exclude position
    
    # Remove columns with zero variance
    var_check <- apply(numeric_cols, 2, var, na.rm = TRUE)
    valid_cols <- names(var_check[!is.na(var_check) & var_check > 0])
    
    if(length(valid_cols) >= 2) {
      cor_data <- cor(numeric_cols[valid_cols], use = "pairwise.complete.obs")
      return(cor_data)
    } else {
      message("Not enough valid variables for correlation analysis")
      return(NULL)
    }
  }
  
  # 1. Position-based Analysis
  pos_analysis <- theta_analysis$detailed_metrics %>%
    mutate(
      position_mb = Position / 1000000,  # Convert to Mb
      region = cut(position_mb, breaks = seq(0, max(position_mb, na.rm=TRUE), by=10),
                   labels = FALSE)
    ) %>%
    group_by(Chr, region) %>%
    summarise(
      start_mb = min(position_mb),
      end_mb = max(position_mb),
      n_snps = n(),
      mean_theta_sd = mean(sd_theta, na.rm=TRUE),
      problem_rate = sum(sd_theta > 0.2, na.rm=TRUE) / n(),
      .groups = 'drop'
    )
  
  # Plot chromosomal heat map
  p1 <- ggplot(pos_analysis, 
               aes(x = start_mb, y = as.factor(Chr), fill = problem_rate)) +
    geom_tile(aes(width = end_mb - start_mb)) +
    scale_fill_viridis() +
    labs(title = "Theta Problem Distribution Across Chromosomes",
         x = "Position (Mb)",
         y = "Chromosome",
         fill = "Problem Rate") +
    theme_minimal()
  
  # 2. Modified Correlation Analysis
  message("Performing correlation analysis on valid metrics...")
  cor_data <- clean_cor_matrix(theta_analysis$detailed_metrics)
  
  if(!is.null(cor_data)) {
    p2 <- ggplot(data = reshape2::melt(cor_data),
                 aes(x = Var1, y = Var2, fill = value)) +
      geom_tile() +
      scale_fill_viridis() +
      geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
      labs(title = "Correlation Between Valid Theta Metrics") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(adv_dir, "metric_correlations.png"), p2, 
           width = 10, height = 10)
  }
  
  # 3. Cluster Pattern Analysis
  cluster_patterns <- theta_analysis$detailed_metrics %>%
    mutate(
      aa_bb_ratio = aa_count / bb_count,
      cluster_pattern = case_when(
        is.na(aa_bb_ratio) ~ "Missing Data",
        aa_bb_ratio > 2 ~ "AA Dominant",
        aa_bb_ratio < 0.5 ~ "BB Dominant",
        TRUE ~ "Balanced"
      )
    ) %>%
    group_by(Chr, cluster_pattern) %>%
    summarise(
      count = n(),
      mean_theta_sd = mean(sd_theta, na.rm=TRUE),
      .groups = 'drop'
    )
  
  # 4. Technical Quality Metrics
  tech_metrics <- theta_analysis$detailed_metrics %>%
    summarise(
      total_snps = n(),
      high_sd = sum(sd_theta > 0.2, na.rm=TRUE),
      moderate_sd = sum(sd_theta >= 0.15 & sd_theta <= 0.2, na.rm=TRUE),
      low_sd = sum(sd_theta < 0.15, na.rm=TRUE),
      missing_data = sum(is.na(sd_theta)),
      
      # Calculate percentages
      pct_high_sd = high_sd/total_snps * 100,
      pct_moderate_sd = moderate_sd/total_snps * 100,
      pct_low_sd = low_sd/total_snps * 100,
      pct_missing = missing_data/total_snps * 100
    )
  
  # 5. Pattern Recognition
  pattern_analysis <- theta_analysis$detailed_metrics %>%
    mutate(
      theta_pattern = case_when(
        sd_theta > 0.2 ~ "High Variance",
        sd_theta > 0.15 & sd_theta <= 0.2 ~ "Moderate Variance",
        sd_theta <= 0.15 ~ "Low Variance",
        TRUE ~ "Missing Data"
      )
    ) %>%
    group_by(Chr, theta_pattern) %>%
    summarise(
      count = n(),
      mean_sd = mean(sd_theta, na.rm=TRUE),
      .groups = 'drop'
    )
  
  # 6. Chromosome-specific Analysis
  chr_analysis <- theta_analysis$detailed_metrics %>%
    group_by(Chr) %>%
    summarise(
      total_snps = n(),
      mean_sd = mean(sd_theta, na.rm=TRUE),
      median_sd = median(sd_theta, na.rm=TRUE),
      problem_snps = sum(sd_theta > 0.2, na.rm=TRUE),
      problem_rate = problem_snps/total_snps * 100,
      .groups = 'drop'
    ) %>%
    arrange(desc(problem_rate))
  
  # Save all plots
  ggsave(file.path(adv_dir, "chromosomal_distribution.png"), p1, 
         width = 12, height = 8)
  
  # Save all analysis results
  write.csv(pos_analysis, 
            file.path(adv_dir, "position_analysis.csv"), 
            row.names = FALSE)
  
  write.csv(cluster_patterns, 
            file.path(adv_dir, "cluster_patterns.csv"), 
            row.names = FALSE)
  
  write.csv(tech_metrics, 
            file.path(adv_dir, "technical_metrics.csv"), 
            row.names = FALSE)
  
  write.csv(pattern_analysis, 
            file.path(adv_dir, "pattern_analysis.csv"), 
            row.names = FALSE)
  
  write.csv(chr_analysis,
            file.path(adv_dir, "chromosome_analysis.csv"),
            row.names = FALSE)
  
  # Create comprehensive summary report
  summary_report <- data.frame(
    Metric = c(
      "Total SNPs Analyzed",
      "High Variance SNPs (%)",
      "Moderate Variance SNPs (%)",
      "Low Variance SNPs (%)",
      "Most Common Pattern",
      "Worst Affected Chromosome",
      "Average Theta SD",
      "Maximum Theta SD",
      "AA Dominant Patterns",
      "BB Dominant Patterns",
      "Balanced Patterns"
    ),
    Value = c(
      tech_metrics$total_snps,
      sprintf("%.1f%%", tech_metrics$pct_high_sd),
      sprintf("%.1f%%", tech_metrics$pct_moderate_sd),
      sprintf("%.1f%%", tech_metrics$pct_low_sd),
      names(which.max(table(pattern_analysis$theta_pattern))),
      chr_analysis$Chr[1],
      sprintf("%.3f", mean(theta_analysis$detailed_metrics$sd_theta, na.rm=TRUE)),
      sprintf("%.3f", max(theta_analysis$detailed_metrics$sd_theta, na.rm=TRUE)),
      sum(cluster_patterns$cluster_pattern == "AA Dominant"),
      sum(cluster_patterns$cluster_pattern == "BB Dominant"),
      sum(cluster_patterns$cluster_pattern == "Balanced")
    )
  )
  
  write.csv(summary_report,
            file.path(adv_dir, "summary_report.csv"),
            row.names = FALSE)
  
  message("\nAdvanced analysis complete! Check the 'Advanced_Insights' folder.")
  message("\nKey findings from the analysis:")
  print(summary_report)
  
  return(list(
    position_analysis = pos_analysis,
    correlations = if(!is.null(cor_data)) cor_data else "No valid correlations",
    cluster_patterns = cluster_patterns,
    technical_metrics = tech_metrics,
    pattern_analysis = pattern_analysis,
    chromosome_analysis = chr_analysis,
    summary_report = summary_report
  ))
}

advanced_insights <- analyze_theta_patterns_advanced(theta_analysis, homo_results, metrics, output_dir)

############# PART 7 ##############
##################################
##################################

analyze_theta_final_insights <- function(theta_analysis, homo_results, metrics, output_dir) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  message("Performing final detailed theta analysis...")
  
  # Create directory
  final_dir <- file.path(output_dir, "Theta_Analysis", "Final_Insights")
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 1. Advanced Theta Distribution Analysis
  theta_distribution <- theta_analysis$raw_problem_data %>%
    group_by(Name) %>%
    summarise(
      # Basic theta stats
      mean_theta = mean(Theta, na.rm = TRUE),
      median_theta = median(Theta, na.rm = TRUE),
      q1_theta = quantile(Theta, 0.25, na.rm = TRUE),
      q3_theta = quantile(Theta, 0.75, na.rm = TRUE),
      iqr_theta = q3_theta - q1_theta,
      
      # Cluster-specific stats
      aa_spread = IQR(Theta[GT == "AA"], na.rm = TRUE),
      bb_spread = IQR(Theta[GT == "BB"], na.rm = TRUE),
      
      # Distribution shape
      skewness = if(sum(!is.na(Theta)) >= 3) moments::skewness(Theta, na.rm = TRUE) else NA_real_,
      kurtosis = if(sum(!is.na(Theta)) >= 4) moments::kurtosis(Theta, na.rm = TRUE) else NA_real_
    ) %>%
    mutate(
      distribution_type = case_when(
        abs(skewness) > 1 ~ "Skewed",
        kurtosis > 4 ~ "Heavy-tailed",
        kurtosis < 2 ~ "Light-tailed",
        TRUE ~ "Normal-like"
      )
    )
  
  # 2. Linkage Analysis
  linkage_patterns <- theta_analysis$detailed_metrics %>%
    arrange(Chr, Position) %>%
    group_by(Chr) %>%
    mutate(
      next_snp_dist = c(diff(Position), NA),
      prev_prob = lag(sd_theta > 0.2),
      next_prob = lead(sd_theta > 0.2),
      prob_cluster = (sd_theta > 0.2 & (prev_prob | next_prob))
    ) %>%
    summarise(
      total_snps = n(),
      clustered_problems = sum(prob_cluster, na.rm = TRUE),
      avg_snp_spacing = mean(next_snp_dist, na.rm = TRUE),
      problem_clusters = sum(prob_cluster & !lag(prob_cluster), na.rm = TRUE),
      .groups = 'drop'
    )
  
  # 3. Sample-based Analysis
  sample_patterns <- theta_analysis$raw_problem_data %>%
    group_by(Sample) %>%
    summarise(
      n_snps = n(),
      mean_theta = mean(Theta, na.rm = TRUE),
      theta_sd = sd(Theta, na.rm = TRUE),
      aa_prop = sum(GT == "AA", na.rm = TRUE) / n(),
      bb_prop = sum(GT == "BB", na.rm = TRUE) / n(),
      ab_prop = sum(GT == "AB", na.rm = TRUE) / n(),
      missing_prop = sum(is.na(GT)) / n()
    ) %>%
    mutate(
      sample_quality = case_when(
        missing_prop > 0.2 ~ "Poor",
        theta_sd > 0.3 ~ "Poor",
        theta_sd > 0.2 ~ "Moderate",
        TRUE ~ "Good"
      )
    )
  
  # Create visualizations
  
  # 1. Theta Distribution Plot
  p1 <- ggplot(theta_distribution, aes(x = mean_theta, fill = distribution_type)) +
    geom_histogram(bins = 50, alpha = 0.7) +
    facet_wrap(~distribution_type) +
    labs(title = "Theta Distribution Patterns",
         x = "Mean Theta",
         y = "Count") +
    theme_minimal()
  
  # 2. Problem Clustering Plot
  p2 <- ggplot(linkage_patterns, 
               aes(x = as.factor(Chr), y = clustered_problems/total_snps)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(title = "Proportion of Clustered Problem SNPs by Chromosome",
         x = "Chromosome",
         y = "Proportion of Clustered Problems") +
    theme_minimal()
  
  # 3. Sample Quality Plot
  p3 <- ggplot(sample_patterns, aes(x = theta_sd, y = missing_prop, 
                                    color = sample_quality)) +
    geom_point(alpha = 0.6) +
    labs(title = "Sample Quality Distribution",
         x = "Theta SD",
         y = "Missing Data Proportion") +
    theme_minimal()
  
  # Save plots
  ggsave(file.path(final_dir, "theta_distributions.png"), p1, width = 12, height = 8)
  ggsave(file.path(final_dir, "problem_clustering.png"), p2, width = 10, height = 6)
  ggsave(file.path(final_dir, "sample_quality.png"), p3, width = 10, height = 6)
  
  # Save analysis results
  write.csv(theta_distribution,
            file.path(final_dir, "detailed_theta_distribution.csv"),
            row.names = FALSE)
  
  write.csv(linkage_patterns,
            file.path(final_dir, "linkage_patterns.csv"),
            row.names = FALSE)
  
  write.csv(sample_patterns,
            file.path(final_dir, "sample_patterns.csv"),
            row.names = FALSE)
  
  # Create summary report
  summary_stats <- data.frame(
    Metric = c(
      "Distribution Types Found",
      "Most Common Distribution",
      "Chromosomes with Problem Clusters",
      "Average Problem Cluster Size",
      "Samples with Poor Quality",
      "Best Performing Sample",
      "Worst Performing Sample"
    ),
    Value = c(
      paste(unique(theta_distribution$distribution_type), collapse = ", "),
      names(which.max(table(theta_distribution$distribution_type))),
      sum(linkage_patterns$problem_clusters > 0),
      sprintf("%.1f SNPs", mean(linkage_patterns$clustered_problems/
                                  linkage_patterns$problem_clusters, na.rm = TRUE)),
      sum(sample_patterns$sample_quality == "Poor"),
      sample_patterns$Sample[which.min(sample_patterns$theta_sd)],
      sample_patterns$Sample[which.max(sample_patterns$theta_sd)]
    )
  )
  
  write.csv(summary_stats,
            file.path(final_dir, "final_summary.csv"),
            row.names = FALSE)
  
  message("\nFinal analysis complete! Check the 'Final_Insights' folder.")
  message("\nKey findings from the final analysis:")
  print(summary_stats)
  
  return(list(
    theta_distribution = theta_distribution,
    linkage_patterns = linkage_patterns,
    sample_patterns = sample_patterns,
    summary_stats = summary_stats
  ))
}

final_insights <- analyze_theta_final_insights(theta_analysis, homo_results, metrics, output_dir)
