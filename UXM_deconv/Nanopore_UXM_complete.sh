#!/bin/bash
# Enhanced Nanopore UXM Deconvolution Pipeline
# Handles Nanopore-specific methylation calling issues including read-end trimming

set -eu

echo "=== Nanopore-Enhanced UXM Deconvolution Pipeline ==="
echo "Date: $(date)"
echo ""

# Configuration
UXM_DIR="/home/gissu/UXM_deconv"
OUTPUT_DIR="nanopore_uxm_results_$(date +%Y%m%d_%H%M%S)"

# Nanopore-specific parameters
DEFAULT_MIN_CPGS=1      # Lowered for Nanopore - more lenient for longer reads
DEFAULT_MIN_MAPQ=10
DEFAULT_NP_THRESH=0.3  # Nanopore probability threshold
DEFAULT_TRIM_ENDS=27    # Trim 27nt from read ends as recommended

# Function to show usage
show_usage() {
    echo "Usage: $0 [options] input_directory OR individual_files..."
    echo ""
    echo "Input options:"
    echo "  directory/                Directory containing .bam or .pat.gz files"
    echo "  file1.bam file2.bam      Individual BAM files"
    echo "  file1.pat.gz file2.pat.gz Individual PAT.gz files"
    echo ""
    echo "Nanopore-specific options:"
    echo "  -o, --output DIR          Output directory (default: auto-generated)"
    echo "  -a, --atlas FILE          Atlas file (default: auto-detect hg38)"
    echo "  --min-cpgs NUM            Minimum CpGs per read (default: $DEFAULT_MIN_CPGS)"
    echo "  --min-mapq NUM            Minimum mapping quality (default: $DEFAULT_MIN_MAPQ)"
    echo "  --np-thresh NUM           Nanopore probability threshold 0-1 (default: $DEFAULT_NP_THRESH)"
    echo "  --trim-ends NUM           Bases to trim from read ends (default: $DEFAULT_TRIM_ENDS)"
    echo "  --no-trim                 Skip read-end trimming (not recommended)"
    echo "  --skip-bam2pat            Skip BAM to PAT conversion (input must be PAT.gz files)"
    echo "  --pattern PATTERN         File pattern for directory processing (default: *.bam or *.pat.gz)"
    echo "  -h, --help                Show this help"
    echo ""
    echo "Examples:"
    echo "  $0 /path/to/bam_directory/                       # Process all BAMs in directory"
    echo "  $0 --pattern '*BRCA*.bam' /data/bams/            # Process only BRCA BAMs"
    echo "  $0 --trim-ends 30 --np-thresh 0.8 /data/        # High-quality directory processing"
    echo "  $0 --skip-bam2pat /data/pat_files/               # Process PAT.gz directory"
    echo "  $0 sample1.bam sample2.bam                       # Individual files (legacy mode)"
    echo ""
    echo "Note: This pipeline uses wgbstools Nanopore mode (-np) and implements"
    echo "      recommended read-end trimming to improve methylation accuracy."
    echo ""
}

# Parse arguments
ATLAS_FILE=""
MIN_CPGS=$DEFAULT_MIN_CPGS
MIN_MAPQ=$DEFAULT_MIN_MAPQ
NP_THRESH=$DEFAULT_NP_THRESH
TRIM_ENDS=$DEFAULT_TRIM_ENDS
NO_TRIM=0
SKIP_BAM2PAT=0
FILE_PATTERN=""  # For directory processing

while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -a|--atlas)
            ATLAS_FILE="$2"
            shift 2
            ;;
        --min-cpgs)
            MIN_CPGS="$2"
            shift 2
            ;;
        --min-mapq)
            MIN_MAPQ="$2"
            shift 2
            ;;
        --np-thresh)
            NP_THRESH="$2"
            shift 2
            ;;
        --trim-ends)
            TRIM_ENDS="$2"
            shift 2
            ;;
        --no-trim)
            NO_TRIM=1
            shift
            ;;
        --skip-bam2pat)
            SKIP_BAM2PAT=1
            shift
            ;;
        --pattern)
            FILE_PATTERN="$2"
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        -*)
            echo "Unknown option: $1"
            show_usage
            exit 1
            ;;
        *)
            break
            ;;
    esac
done

# Check if we have input files
if [[ $# -eq 0 ]]; then
    echo "Error: No input files or directory specified"
    show_usage
    exit 1
fi

# Function to collect files from directory or individual files
collect_input_files() {
    local input_files=()
    
    for arg in "$@"; do
        if [[ -d "$arg" ]]; then
            echo "Processing directory: $arg"
            
            # Determine file pattern if not specified
            if [[ -z "$FILE_PATTERN" ]]; then
                # Auto-detect based on skip-bam2pat setting
                if [[ "$SKIP_BAM2PAT" -eq 1 ]]; then
                    FILE_PATTERN="*.pat.gz"
                else
                    FILE_PATTERN="*.bam"
                fi
            fi
            
            echo "Using file pattern: $FILE_PATTERN"
            
            # Find files matching pattern (only regular files, exclude directories)
            while IFS= read -r -d '' file; do
                # Double-check it's a regular file and has correct extension
                if [[ -f "$file" && ("$file" == *.bam || "$file" == *.pat.gz) ]]; then
                    # Convert to absolute path immediately
                    abs_file=$(realpath "$file")
                    input_files+=("$abs_file")
                fi
            done < <(find "$arg" -maxdepth 1 -name "$FILE_PATTERN" -type f -print0 2>/dev/null)
            
            if [[ ${#input_files[@]} -eq 0 ]]; then
                echo "Warning: No files matching pattern '$FILE_PATTERN' found in directory: $arg"
            else
                echo "Found ${#input_files[@]} files matching pattern '$FILE_PATTERN'"
            fi
            
        elif [[ -f "$arg" ]]; then
            echo "Adding individual file: $arg"
            # Convert to absolute path immediately
            abs_file=$(realpath "$arg")
            input_files+=("$abs_file")
        else
            echo "Warning: Input not found or not accessible: $arg"
        fi
    done
    
    if [[ ${#input_files[@]} -eq 0 ]]; then
        echo "Error: No valid input files found"
        echo "Check your directory path and file pattern"
        exit 1
    fi
    
    # Return the collected files
    printf '%s\n' "${input_files[@]}"
}

# Collect all input files
echo "Collecting input files..."
mapfile -t ALL_INPUT_FILES < <(collect_input_files "$@")

# Filter out any directory names that might have been included
FILTERED_FILES=()
for file in "${ALL_INPUT_FILES[@]}"; do
    if [[ -f "$file" && ("$file" == *.bam || "$file" == *.pat.gz) ]]; then
        FILTERED_FILES+=("$file")
    fi
done
ALL_INPUT_FILES=("${FILTERED_FILES[@]}")

echo "Total files to process: ${#ALL_INPUT_FILES[@]}"
for file in "${ALL_INPUT_FILES[@]}"; do
    echo "  - $(basename "$file")"
done
echo ""

# Determine input file type from collected files
INPUT_TYPE=""
for input_file in "${ALL_INPUT_FILES[@]}"; do
    if [[ "$input_file" == *.bam ]]; then
        if [[ -z "$INPUT_TYPE" ]]; then
            INPUT_TYPE="bam"
        elif [[ "$INPUT_TYPE" != "bam" ]]; then
            echo "Error: Mixed input file types not supported"
            exit 1
        fi
    elif [[ "$input_file" == *.pat.gz ]]; then
        if [[ -z "$INPUT_TYPE" ]]; then
            INPUT_TYPE="pat"
        elif [[ "$INPUT_TYPE" != "pat" ]]; then
            echo "Error: Mixed input file types not supported"
            echo "Found both .bam and .pat.gz files"
            exit 1
        fi
    else
        echo "Error: Unsupported file type: $input_file (must be .bam or .pat.gz)"
        exit 1
    fi
done

echo "Input file type detected: $INPUT_TYPE"
echo "Processing ${#ALL_INPUT_FILES[@]} $INPUT_TYPE files"
echo ""

# Validate skip-bam2pat option
if [[ "$SKIP_BAM2PAT" -eq 1 && "$INPUT_TYPE" == "bam" ]]; then
    echo "Error: Cannot skip bam2pat conversion with BAM input files"
    exit 1
elif [[ "$SKIP_BAM2PAT" -eq 0 && "$INPUT_TYPE" == "pat" ]]; then
    echo "Info: Input files are already PAT.gz, automatically skipping bam2pat conversion"
    SKIP_BAM2PAT=1
fi

# Check that UXM_deconv exists
if [[ ! -d "$UXM_DIR" ]]; then
    echo "Error: UXM_deconv directory not found at $UXM_DIR"
    echo "Please ensure UXM_deconv is installed or update UXM_DIR variable"
    exit 1
fi

# Setup environment
echo "Setting up environment..."
conda activate pysam_env || {
    echo "Warning: Could not activate pysam_env, continuing with current environment"
}

# Verify wgbstools is available
if ! command -v wgbstools &> /dev/null; then
    echo "Error: wgbstools not found in PATH"
    echo "Please ensure wgbstools is properly installed and in PATH"
    exit 1
fi

echo "✓ wgbstools found: $(which wgbstools)"

# Display Nanopore-specific parameters
if [[ "$SKIP_BAM2PAT" -eq 0 ]]; then
    echo ""
    echo "Nanopore BAM to PAT conversion parameters:"
    echo "  --min_cpg: $MIN_CPGS"
    echo "  --mapq: $MIN_MAPQ"
    echo "  --nanopore: enabled"
    echo "  --np_thresh: $NP_THRESH"
    if [[ "$NO_TRIM" -eq 0 ]]; then
        echo "  --clip: $TRIM_ENDS (trimming read ends)"
    else
        echo "  --clip: 0 (no trimming - not recommended)"
    fi
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
cd "$OUTPUT_DIR"
echo "Working in: $(pwd)"
echo ""

# Auto-detect atlas if not specified
if [[ -z "$ATLAS_FILE" ]]; then
    ATLAS_FILE="$UXM_DIR/supplemental/Atlas.U25.l4.hg38.tsv"
    if [[ ! -f "$ATLAS_FILE" ]]; then
        echo "Warning: Default hg38 atlas not found, trying hg19..."
        ATLAS_FILE="$UXM_DIR/supplemental/Atlas.U25.l4.hg19.tsv"
    fi
fi

echo "Using atlas: $ATLAS_FILE"
echo ""

# Process each input file
INPUT_FILES=()
PAT_FILES=()
FAILED_FILES=()

echo "=== File Processing Phase ==="
file_counter=0

for input_file in "${ALL_INPUT_FILES[@]}"; do
    file_counter=$((file_counter + 1))
    echo "[$file_counter/${#ALL_INPUT_FILES[@]}] Processing: $(basename "$input_file")"
    
    if [[ ! -f "$input_file" ]]; then
        echo "  ❌ File not found: $input_file"
        FAILED_FILES+=("$input_file")
        continue
    fi
    
    # Convert to absolute path
    abs_path=$(realpath "$input_file")
    INPUT_FILES+=("$abs_path")
    
    # Handle BAM to PAT conversion if needed
    if [[ "$SKIP_BAM2PAT" -eq 0 ]]; then
        # Extract filename without extension for output
        basename_file=$(basename "$input_file" .bam)
        pat_output="$OUTPUT_DIR/${basename_file}.pat.gz"
        
        echo "  Converting to PAT format with Nanopore-specific settings..."
        
        # Check if BAM contains Nanopore methylation tags
        has_meth_tags=$(samtools view "$abs_path" | head -1000 | grep -c "MM:" || echo "0")
        
        if [[ "$has_meth_tags" -eq 0 ]]; then
            echo "  ⚠️  WARNING: No MM/ML methylation tags found in BAM file"
            echo "     This suggests methylation tags were lost during alignment"
            echo "     Switching to bisulfite-compatible mode (not ideal for Nanopore)"
            USE_NANOPORE_MODE=0
        else
            echo "  ✅ Methylation tags detected, using full Nanopore mode"
            USE_NANOPORE_MODE=1
        fi
        
        # Ensure output directory exists for this specific conversion
        mkdir -p "$OUTPUT_DIR"
        
        # Build bam2pat command based on available data
        bam2pat_cmd="wgbstools bam2pat"
        bam2pat_cmd+=" --min_cpg $MIN_CPGS"
        bam2pat_cmd+=" -q $MIN_MAPQ"
        
        if [[ "$USE_NANOPORE_MODE" -eq 1 ]]; then
            bam2pat_cmd+=" --nanopore"  # Enable Nanopore mode
            bam2pat_cmd+=" --np_thresh $NP_THRESH"
        else
            echo "     Using bisulfite-compatible mode due to missing methylation tags"
        fi
        
        bam2pat_cmd+=" -f"  # Force overwrite existing files
        
        # Apply read-end trimming unless disabled
        if [[ "$NO_TRIM" -eq 0 ]]; then
            bam2pat_cmd+=" --clip $TRIM_ENDS"
        else
            bam2pat_cmd+=" --clip 0"
        fi
        
        bam2pat_cmd+=" -o $OUTPUT_DIR"
        bam2pat_cmd+=" $abs_path"
        
        # Execute bam2pat conversion (show progress for batch processing)  
        if eval $bam2pat_cmd; then
            if [[ -f "$pat_output" ]]; then
                echo "  ✅ Conversion completed: $(basename "$pat_output")"
                PAT_FILES+=("$pat_output")
                
                # Quick stats for batch processing
                total_reads=$(zcat "$pat_output" | wc -l)
                echo "     Reads: $total_reads"
            else
                echo "  ❌ PAT file not created: $(basename "$pat_output")"
                FAILED_FILES+=("$input_file")
            fi
        else
            echo "  ❌ Conversion failed for: $(basename "$input_file")"
            FAILED_FILES+=("$input_file")
        fi
    else
        # Input is already PAT.gz
        echo "  ✅ Using existing PAT file"
        PAT_FILES+=("$abs_path")
    fi
    echo ""
done

# Summary of file processing
echo "=== File Processing Summary ==="
echo "Total input files: ${#ALL_INPUT_FILES[@]}"
echo "Successfully processed: ${#PAT_FILES[@]}"
echo "Failed files: ${#FAILED_FILES[@]}"

if [[ ${#FAILED_FILES[@]} -gt 0 ]]; then
    echo "Failed files:"
    for failed_file in "${FAILED_FILES[@]}"; do
        echo "  - $(basename "$failed_file")"
    done
fi
echo ""

if [[ ${#PAT_FILES[@]} -eq 0 ]]; then
    echo "Error: No valid PAT files available for UXM processing"
    echo "All file conversions failed"
    exit 1
fi

echo "PAT files ready for UXM processing:"
for pat_file in "${PAT_FILES[@]}"; do
    echo "  $(basename "$pat_file")"
done
echo ""

# Run UXM deconvolution
echo "=== UXM Deconvolution Phase ==="
echo "Running UXM deconvolution on ${#PAT_FILES[@]} files..."
echo ""

cd "$UXM_DIR"
python3 src/uxm.py deconv "${PAT_FILES[@]}" --atlas "$ATLAS_FILE" --output "$OUTPUT_DIR/results.csv" --verbose

# Check if results were created
if [[ -f "$OUTPUT_DIR/results.csv" ]]; then
    echo ""
    echo "✓ UXM deconvolution completed successfully!"
    echo "Results saved to: $OUTPUT_DIR/results.csv"
    echo ""
    echo "Results preview:"
    head -10 "$OUTPUT_DIR/results.csv"
    echo ""
    
    # Nanopore-specific result analysis
    echo "Nanopore-specific result analysis:"
    non_zero_tissues=$(tail -n +2 "$OUTPUT_DIR/results.csv" | awk -F, '$2 > 0.001' | wc -l)
    max_proportion=$(tail -n +2 "$OUTPUT_DIR/results.csv" | cut -d, -f2 | sort -nr | head -1)
    total_prop=$(tail -n +2 "$OUTPUT_DIR/results.csv" | cut -d, -f2 | awk '{sum+=$1} END {printf "%.3f", sum}')
    
    echo "  Tissues with >0.1% contribution: $non_zero_tissues"
    echo "  Maximum tissue proportion: $max_proportion"
    echo "  Total proportion sum: $total_prop"
    
    if (( $(echo "$max_proportion > 0.8" | bc -l) )); then
        echo "  ⚠️  Warning: Single tissue dominance may indicate sample issues"
    fi
    
    # Optional: Create visualization (fix the parameter name)
    echo ""
    echo "Creating visualization..."
    python3 src/uxm.py plot "$OUTPUT_DIR/results.csv" --outpath "$OUTPUT_DIR/results_plot.pdf" 2>/dev/null || {
        echo "Note: Visualization failed (parameter issue with UXM plot command)"
    }
    
    if [[ -f "$OUTPUT_DIR/results_plot.pdf" ]]; then
        echo "✓ Plot saved to: $OUTPUT_DIR/results_plot.pdf"
    fi
    
else
    echo "✗ UXM deconvolution failed - no results file created"
    exit 1
fi

echo ""
echo "=== Nanopore Directory Pipeline Complete ==="
echo "Output directory: $OUTPUT_DIR"
echo ""
echo "Processing Summary:"
echo "  Total input files: ${#ALL_INPUT_FILES[@]}"
echo "  Successfully processed: ${#PAT_FILES[@]}"
echo "  Failed files: ${#FAILED_FILES[@]}"
echo ""
echo "Key improvements for Nanopore data:"
echo "  ✓ Used --nanopore flag for direct methylation calling"
echo "  ✓ Applied --np_thresh $NP_THRESH for probability filtering"
if [[ "$NO_TRIM" -eq 0 ]]; then
    echo "  ✓ Trimmed $TRIM_ENDS bp from read ends (recommended by Nanopore)"
else
    echo "  ⚠️  Skipped read-end trimming (not recommended)"
fi
echo ""
echo "Files created:"
ls -la "$OUTPUT_DIR/"
