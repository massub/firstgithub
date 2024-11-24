# Part 1: Data Processing and Basic Analysis
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

#' Read and process GenomeStudio FDT format
#' @param file_path Path to the FDT file
read_fdt_data <- function(file_path) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("Input file does not exist: ", file_path)
  }
  
  message("Reading raw data...")
  raw_data <- read.delim(file_path, sep = "\t", header = TRUE, check.names = FALSE)
  
  # Get sample names
  sample_cols <- grep("\\.GT$", names(raw_data), value = TRUE)
  if (length(sample_cols) == 0) {
    stop("No genotype columns found in input file")
  }
  
  sample_names <- gsub("\\.GT$", "", sample_cols)
  message(sprintf("Found %d samples", length(sample_names)))
  
  # Initialize list for storing processed data
  long_data_list <- list()
  
  # Process each sample
  message("Processing samples...")
  for(i in seq_along(sample_names)) {
    sample <- sample_names[i]
    message(sprintf("Processing sample %d of %d: %s", i, length(sample_names), sample))
    
    # Check for required columns
    theta_col <- paste0(sample, ".Theta")
    r_col <- paste0(sample, ".R")
    gt_col <- paste0(sample, ".GT")
    
    required_cols <- c(theta_col, r_col, gt_col)
    missing_cols <- required_cols[!required_cols %in% names(raw_data)]
    
    if (length(missing_cols) > 0) {
      stop("Missing required columns for sample ", sample, ": ", 
           paste(missing_cols, collapse = ", "))
    }
    
    # Create data frame for current sample
    sample_df <- data.frame(
      SNP_Name = raw_data$Name,
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
  
  message("Creating SNP metadata...")
  snp_metadata <- data.frame(
    SNP_Name = raw_data$Name,
    Chr = raw_data$Chr,
    Position = raw_data$Position,
    stringsAsFactors = FALSE
  )
  
  message(sprintf("Processed %d SNPs across %d samples", 
                  nrow(snp_metadata), 
                  length(sample_names)))
  
  return(list(
    long_data = long_data,
    metadata = snp_metadata,
    sample_count = length(sample_names)
  ))
}

#' Calculate SNP metrics
#' @param data Processed data list
calculate_snp_metrics <- function(data) {
  message("Calculating metrics for each SNP...")
  
  # Calculate metrics
  metrics <- data$long_data %>%
    group_by(SNP_Name) %>%
    summarise(
      # Basic statistics
      mean_theta = mean(Theta, na.rm = TRUE),
      sd_theta = sd(Theta, na.rm = TRUE),
      mean_r = mean(R, na.rm = TRUE),
      sd_r = sd(R, na.rm = TRUE),
      
      # Call rates
      total_calls = n(),
      null_calls = sum(GT == "NULL" | is.na(GT)),
      call_rate = (total_calls - null_calls) / total_calls * 100,
      
      # Genotype frequencies
      aa_freq = sum(GT == "AA", na.rm = TRUE) / total_calls * 100,
      ab_freq = sum(GT == "AB", na.rm = TRUE) / total_calls * 100,
      bb_freq = sum(GT == "BB", na.rm = TRUE) / total_calls * 100,
      
      # Distribution metrics
      theta_range = diff(range(Theta, na.rm = TRUE)),
      
      .groups = "drop"  # This prevents the grouping warning
    ) %>%
    # Join with metadata
    left_join(data$metadata, by = "SNP_Name")
  
  return(metrics)
}

#' Classify SNP patterns
#' @param metrics Data frame with calculated metrics
classify_patterns <- function(metrics) {
  message("Classifying SNP patterns...")
  
  metrics %>%
    mutate(
      quality_class = case_when(
        call_rate >= 95 & sd_theta < 0.1 ~ "Excellent",
        call_rate >= 90 & sd_theta < 0.15 ~ "Good",
        call_rate >= 85 & sd_theta < 0.2 ~ "Moderate",
        call_rate >= 80 ~ "Poor",
        TRUE ~ "Failed"
      ),
      
      cluster_pattern = case_when(
        sd_theta < 0.1 & call_rate >= 95 ~ "Well-separated",
        sd_theta > 0.2 ~ "Diffuse clusters",
        call_rate < 80 ~ "Poor clustering",
        TRUE ~ "Moderate separation"
      )
    )
}

# Part 2: Visualization and Report Generation Functions

#' Generate comprehensive visualizations
#' @param data Processed data list
#' @param metrics Calculated metrics
#' @param output_dir Output directory
generate_visualizations <- function(data, metrics, output_dir) {
  message("Generating visualizations...")
  fig_dir <- file.path(output_dir, "Figures")
  
  # List to store plots
  plots <- list()
  
  tryCatch({
    # 1. Quality Distribution Plot
    message("Creating quality distribution plot...")
    p1 <- ggplot(metrics, aes(x = call_rate, y = sd_theta)) +
      geom_point(aes(color = quality_class), alpha = 0.6) +
      scale_color_viridis_d() +
      theme_minimal() +
      labs(title = "SNP Quality Distribution",
           x = "Call Rate (%)",
           y = "Theta Standard Deviation") +
      theme(legend.position = "bottom")
    
    ggsave(file.path(fig_dir, "quality_distribution.png"), p1, width = 10, height = 8)
    plots$quality_dist <- p1
    
    # 2. Chromosome Distribution
    if("Chr" %in% names(metrics)) {
      message("Creating chromosome quality plot...")
      p2 <- ggplot(metrics, aes(x = Chr)) +
        geom_bar(aes(fill = quality_class)) +
        scale_fill_viridis_d() +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = "bottom") +
        labs(title = "SNP Distribution by Chromosome",
             x = "Chromosome",
             y = "Count")
      
      ggsave(file.path(fig_dir, "chromosome_distribution.png"), p2, width = 12, height = 8)
      plots$chr_dist <- p2
    }
    
    # 3. Call Rate Distribution
    message("Creating call rate distribution plot...")
    p3 <- ggplot(metrics, aes(x = call_rate)) +
      geom_histogram(aes(fill = quality_class), bins = 50) +
      scale_fill_viridis_d() +
      theme_minimal() +
      labs(title = "Distribution of Call Rates",
           x = "Call Rate (%)",
           y = "Count") +
      theme(legend.position = "bottom")
    
    ggsave(file.path(fig_dir, "call_rate_distribution.png"), p3, width = 10, height = 8)
    plots$call_rate_dist <- p3
    
    # 4. Genotype Frequencies
    message("Creating genotype frequencies plot...")
    genotype_data <- metrics %>%
      select(SNP_Name, aa_freq, ab_freq, bb_freq) %>%
      gather(genotype, frequency, -SNP_Name)
    
    p4 <- ggplot(genotype_data, aes(x = frequency, fill = genotype)) +
      geom_histogram(bins = 50, position = "identity", alpha = 0.6) +
      scale_fill_viridis_d() +
      theme_minimal() +
      labs(title = "Distribution of Genotype Frequencies",
           x = "Frequency (%)",
           y = "Count") +
      theme(legend.position = "bottom")
    
    ggsave(file.path(fig_dir, "genotype_frequencies.png"), p4, width = 10, height = 8)
    plots$genotype_freq <- p4
    
  }, error = function(e) {
    message("Warning: Error in plot generation: ", e$message)
    message("Continuing with analysis...")
  })
  
  return(plots)
}

#' Generate HTML report
#' @param data Processed data list
#' @param metrics Calculated metrics
#' @param plots Generated plots
#' @param output_dir Output directory
generate_html_report <- function(data, metrics, plots, output_dir) {
  message("Generating HTML report...")
  
  # Create report content
  rmd_content <- c(
    "---",
    "title: SNP Quality Analysis Report",
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
    "## Analysis Overview",
    sprintf("* Total SNPs analyzed: %d", nrow(metrics)),
    sprintf("* Total samples: %d", data$sample_count),
    "",
    "## Quality Distribution",
    "```{r}",
    "quality_summary <- data.frame(",
    "  Quality = names(table(metrics$quality_class)),",
    "  Count = as.numeric(table(metrics$quality_class)),",
    "  Percentage = round(100 * as.numeric(table(metrics$quality_class))/length(metrics$quality_class), 2)",
    ")",
    "kable(quality_summary, col.names = c('Quality Class', 'Count', 'Percentage (%)'))",
    "```",
    "",
    "## Key Metrics",
    "```{r}",
    "summary_stats <- metrics %>%",
    "  summarise(",
    "    'Mean Call Rate (%)' = mean(call_rate, na.rm = TRUE),",
    "    'Median Call Rate (%)' = median(call_rate, na.rm = TRUE),",
    "    'Mean Theta SD' = mean(sd_theta, na.rm = TRUE),",
    "    'Failed SNPs (%)' = 100 * sum(quality_class == 'Failed')/n()",
    "  ) %>%",
    "  round(2)",
    "kable(t(summary_stats), col.names = 'Value')",
    "```",
    "",
    "# Detailed Analysis",
    "",
    "## Quality Distribution",
    "![](Figures/quality_distribution.png)",
    "",
    "## Chromosome Distribution",
    "![](Figures/chromosome_distribution.png)",
    "",
    "## Call Rate Distribution",
    "![](Figures/call_rate_distribution.png)",
    "",
    "## Genotype Frequencies",
    "![](Figures/genotype_frequencies.png)",
    "",
    "# SNP Quality Analysis",
    "",
    "## Chromosome-wise Quality Distribution",
    "```{r}",
    "if('Chr' %in% names(metrics)) {",
    "  chr_summary <- metrics %>%",
    "    group_by(Chr) %>%",
    "    summarise(",
    "      Total = n(),",
    "      'Excellent (%)' = round(100 * sum(quality_class == 'Excellent')/n(), 2),",
    "      'Failed (%)' = round(100 * sum(quality_class == 'Failed')/n(), 2),",
    "      'Mean Call Rate' = round(mean(call_rate, na.rm = TRUE), 2)",
    "    )",
    "  kable(chr_summary)",
    "}",
    "```",
    "",
    "# Recommendations",
    "",
    sprintf("1. %d SNPs marked as 'Failed' require review", sum(metrics$quality_class == 'Failed')),
    sprintf("2. %d SNPs marked as 'Poor' should be evaluated", sum(metrics$quality_class == 'Poor')),
    sprintf("3. Overall data quality: %.1f%% SNPs passed quality thresholds", 
            100 * sum(metrics$quality_class %in% c('Excellent', 'Good', 'Moderate'))/nrow(metrics)),
    "",
    "# Analysis Details",
    sprintf("* Analysis Date: %s", Sys.Date()),
    sprintf("* Total Processing Time: [Calculated in final report]")
  )
  
  # Write and render report
  rmd_file <- file.path(output_dir, "snp_analysis_report.Rmd")
  writeLines(rmd_content, rmd_file)
  
  tryCatch({
    rmarkdown::render(rmd_file, 
                      output_file = file.path(output_dir, "snp_analysis_report.html"),
                      quiet = TRUE)
  }, error = function(e) {
    message("Warning: Error in HTML report generation: ", e$message)
    message("Check the Rmd file at: ", rmd_file)
  })
}

# Part 3: Results Export and Main Execution

#' Export detailed results
#' @param data Processed data list
#' @param metrics Calculated metrics
#' @param output_dir Output directory
export_results <- function(data, metrics, output_dir) {
  message("Exporting results...")
  tables_dir <- file.path(output_dir, "Tables")
  
  tryCatch({
    # 1. Detailed metrics export
    message("Exporting detailed metrics...")
    write.csv(metrics, 
              file.path(tables_dir, "snp_detailed_metrics.csv"), 
              row.names = FALSE)
    
    # 2. Quality summary
    message("Creating quality summary...")
    quality_summary <- metrics %>%
      group_by(quality_class) %>%
      summarise(
        count = n(),
        mean_call_rate = mean(call_rate, na.rm = TRUE),
        mean_theta_sd = mean(sd_theta, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(percentage = count / sum(count) * 100)
    
    write.csv(quality_summary, 
              file.path(tables_dir, "quality_summary.csv"), 
              row.names = FALSE)
    
    # 3. Failed SNPs report
    message("Generating failed SNPs report...")
    failed_snps <- metrics %>%
      filter(quality_class %in% c("Poor", "Failed")) %>%
      arrange(call_rate) %>%
      select(SNP_Name, Chr, Position, quality_class, call_rate, 
             sd_theta, aa_freq, ab_freq, bb_freq)
    
    write.csv(failed_snps, 
              file.path(tables_dir, "failed_snps.csv"), 
              row.names = FALSE)
    
    # 4. Chromosome summary if chromosome information exists
    if("Chr" %in% names(metrics)) {
      message("Creating chromosome summary...")
      chr_summary <- metrics %>%
        group_by(Chr) %>%
        summarise(
          total_snps = n(),
          failed_snps = sum(quality_class == "Failed"),
          mean_call_rate = mean(call_rate, na.rm = TRUE),
          .groups = 'drop'
        )
      
      write.csv(chr_summary, 
                file.path(tables_dir, "chromosome_summary.csv"), 
                row.names = FALSE)
    }
    
  }, error = function(e) {
    message("Error in results export: ", e$message)
  })
}

#' Display analysis summary
#' @param metrics Calculated metrics
#' @param data Processed data list
#' @param output_dir Output directory
display_analysis_summary <- function(metrics, data, output_dir) {
  message("\nGenerating analysis summary...")
  
  tryCatch({
    # Calculate summary statistics
    summary_stats <- list(
      total_snps = nrow(metrics),
      total_samples = data$sample_count,
      quality_dist = table(metrics$quality_class),
      mean_call_rate = mean(metrics$call_rate, na.rm = TRUE),
      failed_snps = sum(metrics$quality_class %in% c("Poor", "Failed"))
    )
    
    # Create summary text
    summary_text <- c(
      "\n=== SNP ANALYSIS SUMMARY ===\n",
      sprintf("Total SNPs analyzed: %d", summary_stats$total_snps),
      sprintf("Total samples: %d", summary_stats$total_samples),
      "\nQUALITY DISTRIBUTION:"
    )
    
    for(qual in names(summary_stats$quality_dist)) {
      count <- summary_stats$quality_dist[qual]
      percentage <- (count/summary_stats$total_snps) * 100
      summary_text <- c(summary_text,
                        sprintf("%s: %d (%.1f%%)", qual, count, percentage))
    }
    
    summary_text <- c(
      summary_text,
      sprintf("\nMean call rate: %.1f%%", summary_stats$mean_call_rate),
      sprintf("SNPs requiring attention: %d", summary_stats$failed_snps),
      sprintf("\nResults location: %s", output_dir)
    )
    
    # Print to console
    cat(paste(summary_text, collapse = "\n"))
    
    # Save to file
    writeLines(summary_text, file.path(output_dir, "Tables", "analysis_summary.txt"))
    
  }, error = function(e) {
    message("Error in summary generation: ", e$message)
  })
}

#' Main execution function
#' @param snp_file Path to FDT file
#' @param output_dir Output directory
#' @param sample_limit Optional limit on number of samples to process (for testing)
main <- function(snp_file, output_dir, sample_limit = NULL) {
  # Record start time
  start_time <- Sys.time()
  
  tryCatch({
    # Create output directories
    message("\nCreating output directories...")
    dir.create(file.path(output_dir, "Figures"), recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(output_dir, "Tables"), recursive = TRUE, showWarnings = FALSE)
    
    # Step 1: Read and process data
    message("\nStep 1: Reading and processing data...")
    data <- read_fdt_data(snp_file)
    
    # Apply sample limit if specified
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
    
    # Step 3: Classify patterns
    message("\nStep 3: Classifying patterns...")
    metrics <- classify_patterns(metrics)
    
    # Step 4: Generate visualizations
    message("\nStep 4: Generating visualizations...")
    plots <- generate_visualizations(data, metrics, output_dir)
    
    # Step 5: Export results
    message("\nStep 5: Exporting results...")
    export_results(data, metrics, output_dir)
    
    # Step 6: Generate HTML report
    message("\nStep 6: Generating HTML report...")
    generate_html_report(data, metrics, plots, output_dir)
    
    # Display summary
    display_analysis_summary(metrics, data, output_dir)
    
    # Calculate execution time
    end_time <- Sys.time()
    execution_time <- difftime(end_time, start_time, units = "mins")
    message(sprintf("\nAnalysis completed in %.1f minutes", as.numeric(execution_time)))
    
    # Return results
    return(list(
      metrics = metrics,
      plots = plots,
      summary = list(
        total_snps = nrow(metrics),
        total_samples = data$sample_count,
        quality_distribution = table(metrics$quality_class),
        execution_time = execution_time
      )
    ))
    
  }, error = function(e) {
    message("\nERROR: Analysis failed with the following error:")
    message(e$message)
    message("\nStack trace:")
    print(sys.calls())
    message("\nPlease check:")
    message("1. Input file exists and is readable")
    message("2. Input file format is correct")
    message("3. Output directory is writable")
    message("4. Sufficient memory is available")
    return(structure(e$message, class = "try-error"))
  })
}

# Example usage:
# Set your file paths
snp_file <- "D:/USDA/Small Grains Research Unit/Projects/Project--SNP position and cluster charecterization/Project_Root/00_Raw_Data/Oat_Datasets/SDSU2023Caffe-O3K_01-08-FDT.txt"
output_dir <- "D:/USDA/Small Grains Research Unit/Projects/Project--SNP position and cluster charecterization/Project_Root/03_Results"

# Run the analysis
results <- main(snp_file, output_dir)

# Check results
if(inherits(results, "try-error")) {
  message("Error in analysis. Please check the error message above.")
} else {
  message("Analysis completed successfully!")
  message("\nQuality Distribution:")
  print(results$summary$quality_distribution)
}