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
