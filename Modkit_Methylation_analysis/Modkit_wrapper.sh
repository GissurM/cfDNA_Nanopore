#!/bin/bash

#SBATCH --job-name=base
#SBATCH --partition=mimir  # request node from a specific partition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8 
#SBATCH --hint=multithread     # 48 cores per node (96 in total)
#SBATCH --output=messages.mod1.out.txt   
#SBATCH --error=messages.mod1.err.txt 
#SBATCH --time=15:00:00 

pwd 
date

# Verify that required directories and files exist
echo "Checking required files and directories..."

# Check if output directory exists and create it if not
OUTPUT_DIR="/hpcdata/Mimir/gmi6/coronary_cfDNA/mod_BED_coronary"
if [ ! -d "$OUTPUT_DIR" ]; then
    echo "Creating output directory: $OUTPUT_DIR"
    mkdir -p "$OUTPUT_DIR"
fi

# Check reference genome and modkit binary
REF_GENOME="/hpchome/gmi6/ref/hg38ref/hg38.fa"
MODKIT_BIN="/hpchome/gmi6/ref/dist_modkit_v0.6.1_481e3c9/modkit"

# Check if reference genome exists
if [ ! -f "$REF_GENOME" ]; then
    echo "ERROR: Reference genome file not found: $REF_GENOME"
    echo "Please ensure the reference genome file is available."
    exit 1
fi

# Check if reference genome index exists (commented out - fa.fai already exists)
# if [ ! -f "${REF_GENOME}.fai" ]; then
#     echo "Reference genome index not found: ${REF_GENOME}.fai"
#     echo "Attempting to create index with samtools faidx..."
#     
#     # Check if samtools is available
#     if ! command -v samtools &> /dev/null; then
#         # Try to load samtools module if direct command not available
#         if command -v module &> /dev/null; then
#             echo "Loading samtools module..."
#             module load samtools || true
#         fi
#     fi
#     
#     # Check again if samtools is available
#     if command -v samtools &> /dev/null; then
#         echo "Creating index for reference genome..."
#         samtools faidx "$REF_GENOME"
#         
#         if [ ! -f "${REF_GENOME}.fai" ]; then
#             echo "ERROR: Failed to create index file. Manual indexing required."
#             echo "Please run: samtools faidx $REF_GENOME"
#             exit 1
#         else
#             echo "Successfully created reference genome index."
#         fi
#     else
#         echo "ERROR: samtools not found. Cannot create reference genome index."
#         echo "Please manually index the reference genome:"
#         echo "module load samtools  # if available"
#         echo "samtools faidx $REF_GENOME"
#         exit 1
#     fi
# fi

echo "Reference genome index exists: ${REF_GENOME}.fai"

# Check if modkit binary exists and is executable
if [ ! -x "$MODKIT_BIN" ]; then
    echo "ERROR: modkit binary not found or not executable: $MODKIT_BIN"
    echo "Please ensure the modkit binary is available and has execute permissions."
    exit 1
fi

# Check input BAM files using wildcard pattern
BAM_PATTERN="/hpcdata/Mimir/gmi6/coronary_cfDNA/*.bam"
echo "Searching for BAM files matching pattern: $BAM_PATTERN"

# Count BAM files
BAM_FILES=($BAM_PATTERN)
BAM_COUNT=${#BAM_FILES[@]}

if [ "$BAM_COUNT" -eq 0 ] || [ ! -f "${BAM_FILES[0]}" ]; then
    echo "ERROR: No BAM files found matching pattern: $BAM_PATTERN"
    echo "Please ensure the BAM files exist at the specified locations."
    exit 1
fi

echo "Found $BAM_COUNT BAM files to process."
echo "All required files and directories verified."

# Iterate over all .bam files matching the pattern, excluding .bam.bai files
echo "Starting processing of BAM files..."
for bamfile in $BAM_PATTERN; do
    if [[ $bamfile != *.bam.bai ]]; then
        # Extract the base name of the BAM file (without extension)
        base_name=$(basename "$bamfile" .bam)
        
        echo "Processing: $base_name"
        
        # Run modkit pileup for each BAM file
        # Note: --preset traditional removed in v0.6.0+
        # Using recommended replacement for comparable output
        time "$MODKIT_BIN" pileup "$bamfile" "$OUTPUT_DIR/${base_name}.bed" \
            --ref "$REF_GENOME" \
            --modified-bases 5mC 5hmC \
            --combine-mods \
            --cpg \
            --combine-strands \
            --threads 8
        
        # Check if the command executed successfully
        if [ $? -eq 0 ]; then
            echo "Successfully processed: $base_name"
        else
            echo "ERROR: Failed to process: $base_name"
        fi
    fi
done

date 
