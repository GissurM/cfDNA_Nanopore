# Load necessary libraries
library(Rsamtools)
library(ggplot2)

# Define the directory containing BAM files
bam_dir <- "C:/Users/gissu/Documents/hg38_1-24"  # Replace with the path to your BAM files
bam_files <- list.files(bam_dir, pattern = "\\.bam$", full.names = TRUE)  # List all BAM files (ignore .bam.bai)

# Loop through each BAM file
for (bam_file in bam_files) {
    # Extract a unique identifier from the file name (use basename without extension)
    barcode <- tools::file_path_sans_ext(basename(bam_file))

    # Read BAM file using Rsamtools
    cat("Processing file:", basename(bam_file), "\n")
    
    # Open BAM file
    bf <- BamFile(bam_file)
    
    # Read all alignments and extract chromosome information
    param <- ScanBamParam(what = "rname")
    bam_data <- scanBam(bf, param = param)
    
    # Extract chromosome names (rname = reference name)
    chromosomes <- bam_data[[1]]$rname
    
    # Remove NA values (unmapped reads)
    chromosomes <- chromosomes[!is.na(chromosomes)]
    
    # Count occurrences of each chromosome
    chromosome_table <- table(chromosomes)
    
    # Convert to data frame
    chromosome_counts <- data.frame(
        Chromosome = names(chromosome_table),
        Count = as.integer(chromosome_table),
        stringsAsFactors = FALSE
    )

    # Add a column for the percentage of total counts
    total_count <- sum(chromosome_counts$Count)
   chromosome_counts$Percentage <- round((chromosome_counts$Count / total_count) * 100, digits = 2)

    # Add column names to the data frame
    colnames(chromosome_counts) <- c("Chromosome", "Count", "Percentage")

    # Aggregate counts to ensure only one bar per chromosome
    chromosome_counts <- aggregate(cbind(Count, Percentage) ~ Chromosome, data = chromosome_counts, sum)

    # Update chromosome order to use uppercase names
    chromosome_order <- c(paste0("CHR", 1:22), "CHRX", "CHRY", "CHRM", "CHRUN", "CHREB")

    # Normalize chromosome names to uppercase
    chromosome_counts$Chromosome <- toupper(chromosome_counts$Chromosome)

    # Group unrecognized contigs into "Other"
    chromosome_counts$Chromosome <- ifelse(chromosome_counts$Chromosome %in% chromosome_order,
                                           chromosome_counts$Chromosome, "Other")

    # Order chromosomes
    chromosome_counts$Chromosome <- factor(chromosome_counts$Chromosome, levels = c(chromosome_order, "Other"))
    chromosome_counts <- chromosome_counts[order(chromosome_counts$Chromosome), ]

    # Plot the histogram for the current BAM file
    plot <- ggplot(chromosome_counts, aes(x = Chromosome, y = Count)) +
        geom_bar(stat = "identity", fill = "steelblue", width = 0.8) +
        theme_minimal() +
        labs(title = paste("cfDNA Reads per Chromosome -", barcode), 
             x = "Chromosome", y = "Read Count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    # Save the plot to the directory with the filename-based identifier
    output_file <- file.path(bam_dir, paste0("chromosome_counts_", barcode, ".png"))
    ggsave(output_file, plot)

    # Save a table showing the total amount of each chromosome present in the BAM file
    output_table_file <- file.path(bam_dir, paste0("chromosome_counts_table_", barcode, ".csv"))
    write.table(chromosome_counts, output_table_file, row.names = FALSE, sep = ";", dec = ",", col.names = TRUE, quote = FALSE)
}
