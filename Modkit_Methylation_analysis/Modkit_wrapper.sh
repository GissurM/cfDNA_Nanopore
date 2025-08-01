#!/bin/bash

#SBATCH --job-name=base
#SBATCH --partition=mimir  # request node from a specific partition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32 
#SBATCH --hint=multithread     # 48 cores per node (96 in total)
#SBATCH --output=messages.mod1.out.txt   
#SBATCH --error=messages.mod1.err.txt 
#SBATCH --time=05:00:00 

pwd 
date

# Check reference genome and modkit binary
OUTPUT_DIR="/path/to/your/output" # The path should remain in quotes this also applies for REF_GENOME and MODKIT_BIN.
REF_GENOME="/path/to/your/reference.fa"
INPUT_DIR="/path/to/your/bamfiledirectory" # The wrapper is made to run through all .bam files in a directory and requires a directory input a simpler command is possible if just running one file. 
MODKIT_BIN="/path/to/your/modkit.exe" #not needed if running modkit as a module


# Check if reference genome index exists, create it if it doesn't
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Reference genome index not found: ${REF_GENOME}.fai"
    echo "Attempting to create index with samtools faidx..."
    
    # Check if samtools is available
    if ! command -v samtools &> /dev/null; then
        # Try to load samtools module if direct command not available
        if command -v module &> /dev/null; then
            echo "Loading samtools module..."
            module load samtools || true
        fi
    fi
    
    # Check again if samtools is available
    if command -v samtools &> /dev/null; then
        echo "Creating index for reference genome..."
        samtools faidx "$REF_GENOME"
        
        if [ ! -f "${REF_GENOME}.fai" ]; then
            echo "ERROR: Failed to create index file. Manual indexing required."
            echo "Please run: samtools faidx $REF_GENOME"
            exit 1
        else
            echo "Successfully created reference genome index."
        fi
    else
        echo "ERROR: samtools not found. Cannot create reference genome index."
        echo "Please manually index the reference genome:"
        echo "module load samtools  # if available"
        echo "samtools faidx $REF_GENOME"
        exit 1
    fi
fi

# Check input BAM files
BAM_COUNT=$(ls "$INPUT_DIR"/*.bam 2>/dev/null | wc -l)
if [ "$BAM_COUNT" -eq 0 ]; then
    echo "ERROR: No BAM files found in: $INPUT_DIR"
    echo "Please ensure the input directory contains BAM files."
    exit 1
fi

echo "Found $BAM_COUNT BAM files to process."
echo "All required files and directories verified."

# Iterate over all .bam files in the directory, excluding .bam.bai files
echo "Starting processing of BAM files..."
for bamfile in "$INPUT_DIR"/*.bam; do
    if [[ $bamfile != *.bam.bai ]]; then
        # Extract the base name of the BAM file (without extension)
        base_name=$(basename "$bamfile" .bam)
        
        echo "Processing: $base_name"
        
        # Run modkit pileup for each BAM file
        time "$MODKIT_BIN" pileup "$bamfile" "$OUTPUT_DIR/${base_name}.bed" \
            --ref "$REF_GENOME" \
            --preset traditional  # The traditional preset is one of many presets that Modkit has. It is useful when you are only interested in 5mC methylation at CpG sites
        
        # Check if the command executed successfully
        if [ $? -eq 0 ]; then
            echo "Successfully processed: $base_name"
        else
            echo "ERROR: Failed to process: $base_name"
        fi
    fi
done

date 
