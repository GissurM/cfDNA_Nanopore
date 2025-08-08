# This code is simple to use and does not rely on many dependencies. 
# It runs with a clear verbose output which explains in real time how the run is going.
# The output directory location is poorly defined and it is likely to appear in your user home directory on windows probably something like "C:\Users\"your-name"\"
library(Rsamtools)  
library(ggplot2)    
library(dplyr)      

# Define the directory containing BAM files
base_dir <- "/your/bam/dir/path"  # Adjust this to your BAM files directory

# Filter settings
min_qwidth_filtered <- 100   # Minimum fragment size for filtered analysis
max_qwidth_filtered <- 1000  # Maximum fragment size for filtered analysis

# cfDNA size classification standards (in base pairs)
mono_min <- 50; mono_max <- 200      # Mononucleosome fragments
di_min <- 201; di_max <- 400         # Dinucleosome fragments  
tri_min <- 401; tri_max <- 600       # Trinucleosome fragments
hmw_min <- 601; hmw_max <- 2000      # High molecular weight fragments

# Helper function to list available BAM files (excluding index files)
list_bam_files <- function(dir = base_dir) {
  all_files <- list.files(dir, pattern = "\\.bam$", full.names = TRUE)
  bam_files <- all_files[!endsWith(all_files, ".bai")]
  return(bam_files)
}

# Function to classify cfDNA fragments
classify_fragments <- function(qwidths) {
  mono <- sum(qwidths >= mono_min & qwidths <= mono_max, na.rm = TRUE)
  di <- sum(qwidths >= di_min & qwidths <= di_max, na.rm = TRUE)
  tri <- sum(qwidths >= tri_min & qwidths <= tri_max, na.rm = TRUE)
  hmw <- sum(qwidths >= hmw_min & qwidths <= hmw_max, na.rm = TRUE)
  
  total <- length(qwidths[!is.na(qwidths)])
  
  return(list(
    mono = mono,
    di = di,
    tri = tri,
    hmw = hmw,
    mono_pct = (mono / total) * 100,
    di_pct = (di / total) * 100,
    tri_pct = (tri / total) * 100,
    hmw_pct = (hmw / total) * 100
  ))
}

# Function to extract qwidth values and create analysis for a single BAM file
analyze_bam_file <- function(file_path) {
  file_name <- tools::file_path_sans_ext(basename(file_path))
  
  # Check if file exists
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  # Extra check to ensure it's actually a BAM file and not an index file
  if (endsWith(file_path, ".bai") || !endsWith(file_path, ".bam")) {
    warning(paste("Not a BAM file or possibly an index file:", file_path))
    return(NULL)
  }
  
  cat("Processing:", file_name, "\n")
  
  # Read BAM file - we only need the qwidth values
  param <- ScanBamParam(what=c("qwidth"))
  tryCatch({
    bam_data <- scanBam(file_path, param=param)[[1]]
    qwidths_all <- bam_data$qwidth[!is.na(bam_data$qwidth)]
    
    # Filter qwidth values to the specified range for focused analysis
    qwidths_filtered <- qwidths_all[qwidths_all >= min_qwidth_filtered & qwidths_all <= max_qwidth_filtered]
    
    # Log the filtering results
    total_reads <- length(qwidths_all)
    filtered_reads <- length(qwidths_filtered)
    filtered_pct <- (filtered_reads / total_reads) * 100
    
    cat(sprintf("  %d total reads, %d (%.1f%%) in %d-%d bp range\n", 
                total_reads, filtered_reads, filtered_pct, min_qwidth_filtered, max_qwidth_filtered))
    
    # Calculate statistics for both all data and filtered data
    stats_all <- list(
      mean = mean(qwidths_all, na.rm = TRUE),
      median = median(qwidths_all, na.rm = TRUE),
      sd = sd(qwidths_all, na.rm = TRUE),
      count = length(qwidths_all)
    )
    
    stats_filtered <- list(
      mean = mean(qwidths_filtered, na.rm = TRUE),
      median = median(qwidths_filtered, na.rm = TRUE),
      sd = sd(qwidths_filtered, na.rm = TRUE),
      count = length(qwidths_filtered)
    )
    
    # Classify fragments (using all data for biological relevance)
    classification <- classify_fragments(qwidths_all)
    
    # Create histogram data for plotting (5 bp bins from 100-1000 bp)
    bins <- seq(min_qwidth_filtered, max_qwidth_filtered, by = 5)
    hist_data <- hist(qwidths_filtered, breaks = bins, plot = FALSE)
    
    # Create bar plot
    plot_data <- data.frame(
      size = hist_data$mids,
      count = hist_data$counts
    )
    
    p <- ggplot(plot_data, aes(x = size, y = count)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7, width = 4) +
      labs(
        title = paste("cfDNA Fragment Size Distribution:", file_name),
        subtitle = paste("Filtered to", min_qwidth_filtered, "-", max_qwidth_filtered, "bp (5 bp bins)"),
        x = "Fragment Size (bp)",
        y = "Count"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      scale_x_continuous(breaks = seq(100, 1000, by = 100))
    
    # Save the plot
    output_dir <- "cfDNA_individual_analysis"
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    
    plot_file <- file.path(output_dir, paste0(file_name, "_size_distribution.png"))
    ggsave(plot_file, plot = p, width = 10, height = 6, dpi = 300)
    cat("  Saved plot:", basename(plot_file), "\n")
    
    # Return all the data for compilation
    return(list(
      file_name = file_name,
      stats_all = stats_all,
      stats_filtered = stats_filtered,
      classification = classification,
      plot_data = plot_data
    ))
    
  }, error = function(e) {
    warning(paste("Error reading file:", file_path, "-", e$message))
    return(NULL)
  })
}

# Main analysis
cat("Starting cfDNA individual file analysis...\n")
cat("Looking for BAM files in directory:", base_dir, "\n")

available_bam_files <- list_bam_files()
if (length(available_bam_files) == 0) {
  stop("No BAM files found in the directory. Please check the path.")
}

cat("Found", length(available_bam_files), "BAM files:\n")
for (file in available_bam_files) {
  cat("  ", basename(file), "\n")
}

# Process all files
cat("\nProcessing files...\n")
results_list <- list()

for (file_path in available_bam_files) {
  result <- analyze_bam_file(file_path)
  if (!is.null(result)) {
    results_list[[result$file_name]] <- result
  }
}

# Compile results into a comprehensive CSV
if (length(results_list) > 0) {
  cat("\nCompiling results...\n")
  
  # Create data frame for CSV output
  compiled_data <- data.frame(
    File_Name = character(),
    
    # All data statistics
    Mean_All = numeric(),
    Median_All = numeric(),
    SD_All = numeric(),
    Count_All = numeric(),
    
    # Filtered data statistics (100-1000 bp)
    Mean_Filtered = numeric(),
    Median_Filtered = numeric(),
    SD_Filtered = numeric(),
    Count_Filtered = numeric(),
    
    # cfDNA classification counts
    Mono_Count = numeric(),
    Di_Count = numeric(),
    Tri_Count = numeric(),
    HMW_Count = numeric(),
    
    # cfDNA classification percentages
    Mono_Percent = numeric(),
    Di_Percent = numeric(),
    Tri_Percent = numeric(),
    HMW_Percent = numeric(),
    
    stringsAsFactors = FALSE
  )
  
  for (file_name in names(results_list)) {
    result <- results_list[[file_name]]
    
    new_row <- data.frame(
      File_Name = file_name,
      
      # All data statistics
      Mean_All = round(result$stats_all$mean, 2),
      Median_All = result$stats_all$median,
      SD_All = round(result$stats_all$sd, 2),
      Count_All = result$stats_all$count,
      
      # Filtered data statistics
      Mean_Filtered = round(result$stats_filtered$mean, 2),
      Median_Filtered = result$stats_filtered$median,
      SD_Filtered = round(result$stats_filtered$sd, 2),
      Count_Filtered = result$stats_filtered$count,
      
      # cfDNA classification counts
      Mono_Count = result$classification$mono,
      Di_Count = result$classification$di,
      Tri_Count = result$classification$tri,
      HMW_Count = result$classification$hmw,
      
      # cfDNA classification percentages
      Mono_Percent = round(result$classification$mono_pct, 2),
      Di_Percent = round(result$classification$di_pct, 2),
      Tri_Percent = round(result$classification$tri_pct, 2),
      HMW_Percent = round(result$classification$hmw_pct, 2),
      
      stringsAsFactors = FALSE
    )
    
    compiled_data <- rbind(compiled_data, new_row)
  }
  
  # Save the compiled CSV
  output_dir <- "cfDNA_individual_analysis"
  csv_file <- file.path(output_dir, "cfDNA_compiled_analysis.csv")
  write.csv(compiled_data, csv_file, row.names = FALSE)
  
  cat("Saved compiled analysis to:", csv_file, "\n")
  
  # Print summary
  cat("\n=== ANALYSIS SUMMARY ===\n")
  cat("Processed", nrow(compiled_data), "BAM files\n")
  cat("cfDNA size classifications used:\n")
  cat("  Mononucleosome:", mono_min, "-", mono_max, "bp\n")
  cat("  Dinucleosome:", di_min, "-", di_max, "bp\n")
  cat("  Trinucleosome:", tri_min, "-", tri_max, "bp\n")
  cat("  High Molecular Weight:", hmw_min, "-", hmw_max, "bp\n")
  cat("\nFiltered analysis range:", min_qwidth_filtered, "-", max_qwidth_filtered, "bp\n")
  cat("Bar plot bin size: 5 bp\n")
  cat("\nOutput files saved in:", output_dir, "/\n")
  cat("- Individual PNG plots for each BAM file\n")
  cat("- Compiled CSV with all statistics\n")
  
} else {
  cat("No data was successfully extracted from any BAM files.\n")
}

cat("\nAnalysis complete!\n")
