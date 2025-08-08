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

    # Convert Windows path to WSL/Linux path
    linux_path <- system(paste("wsl wslpath", shQuote(bam_file)), intern = TRUE)

    # Construct the samtools command to directly count chromosomes
    command <- paste("wsl samtools view", shQuote(linux_path), "| awk '{print $3}' | sort | uniq -c")

    # Execute the command and capture the output
    output <- system(command, intern = TRUE)

    # Parse the output to create a data frame with separate columns for Chromosome and Count
    chromosome_counts <- do.call(rbind, lapply(output, function(line) {
        parts <- strsplit(line, " ")[[1]]
        parts <- parts[parts != ""]  # Remove empty strings
        data.frame(Chromosome = parts[2], Count = as.integer(parts[1]), stringsAsFactors = FALSE)
    }))

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
        labs(title = paste("cfDNA Reads per Chromosome - Barcode", barcode), 
             x = "Chromosome", y = "Read Count") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

    # Save the plot to the directory with the name chromosome_counts_barcodeXX.png
    output_file <- file.path(bam_dir, paste0("chromosome_counts_barcode", barcode, ".png"))
    ggsave(output_file, plot)

    # Save a table showing the total amount of each header present in the BAM file
    output_table_file <- file.path(bam_dir, paste0("chromosome_counts_table_barcode", barcode, ".csv"))
    write.table(chromosome_counts, output_table_file, row.names = FALSE, sep = ";", dec = ",", col.names = TRUE, quote = FALSE)
}
