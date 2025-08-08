#!/bin/bash

#SBATCH --job-name=NanoTelo
#SBATCH --partition=mimir  # request node from a specific partition
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32 
#SBATCH --mem=32G  # Request 32GB memory 
#SBATCH --hint=multithread     # 48 cores per node (96 in total)
#SBATCH --output=messages.Telo.out.txt   
#SBATCH --error=messages.Telo.err.txt 
#SBATCH --time=08:00:00 

# HPC Wrapper Script for Telomere Analysis
# Usage: ./Run_Telomeric_analysis.sh [options]

# Set default values
INPUT_DIR=""
OUTPUT_DIR=""
THREADS=8
CHUNK_SIZE=1000
MIN_TELOMERE_LENGTH=40
MAX_ERROR_RATE=0.1
PREFIX="STRTelomere"
PYTHON_VERSION="3.9"
VENV_NAME="telomere_env"

# Determine script directory - handle SLURM job environment
if [[ -n "$SLURM_SUBMIT_DIR" ]]; then
    # In SLURM job, use the submission directory
    SCRIPT_DIR="$SLURM_SUBMIT_DIR"
elif [[ -n "${BASH_SOURCE[0]}" ]]; then
    # Try standard method
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
else
    # Fallback to current directory
    SCRIPT_DIR="$(pwd)"
fi

# Use absolute paths for the Python script and requirements
SCRIPT_PATH="${SCRIPT_DIR}/Telomeric_analysis.py"
REQUIREMENTS_PATH="${SCRIPT_DIR}/requirements.txt"

# Function to display usage
usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required Options:
    -i, --input-dir DIR          Directory containing BAM files
    -o, --output-dir DIR         Output directory for results

Optional Parameters:
    -t, --threads INT            Number of processing threads (default: 8)
    -c, --chunk-size INT         Number of reads per chunk (default: 1000)
    -l, --min-length INT         Minimum telomeric length in bp (default: 40)
    -e, --error-rate FLOAT       Maximum error rate (default: 0.1)
    -p, --prefix STRING          Output file prefix (default: STRTelomere)

Environment Options:
    --python-version VERSION     Python version to use (default: 3.9)
    --venv-name NAME            Virtual environment name (default: telomere_env)
    --skip-env-setup           Skip environment setup (use existing)
    --force-env-rebuild        Force rebuild of virtual environment
    --system-python            Use system Python instead of virtual environment
    --minimal-modules          Load only essential modules (faster startup)

Other Options:
    -h, --help                  Show this help message
    --dry-run                   Show commands that would be executed
    --verbose                   Enable verbose output
    --test-env                  Test environment setup only (don't run analysis)

Examples:
    # Basic usage
    $0 -i /data/bam_files -o /results/telomere_output

    # With custom parameters
    $0 -i /data/bams -o /results -t 16 -c 2000 --min-length 50

    # For large datasets on HPC
    $0 -i /scratch/bams -o /scratch/results -t 32 -c 5000

EOF
}

# Function for logging
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function for error handling
error_exit() {
    echo "[ERROR] $1" >&2
    cleanup_and_exit 1
}

# Function for cleanup and exit
cleanup_and_exit() {
    local exit_code=${1:-0}
    
    log "Cleanup initiated (exit code: $exit_code)"
    
    # Clean up virtual environment unless explicitly requested to keep it
    if [[ "$KEEP_VENV" == "false" ]] && [[ -n "$VENV_PATH" ]] && [[ -d "$VENV_PATH" ]]; then
        log "Cleaning up virtual environment: $VENV_PATH"
        
        # Deactivate if active
        if [[ -n "$VIRTUAL_ENV" ]]; then
            log "Deactivating virtual environment..."
            deactivate 2>/dev/null || true
        fi
        
        # Remove virtual environment directory
        if [[ "$DRY_RUN" == "true" ]]; then
            echo "[DRY-RUN] rm -rf '$VENV_PATH'"
        else
            if rm -rf "$VENV_PATH" 2>/dev/null; then
                log "✓ Virtual environment cleaned up successfully"
            else
                log "⚠ Warning: Could not remove virtual environment (permissions?)"
            fi
        fi
    elif [[ "$KEEP_VENV" == "true" ]]; then
        log "Keeping virtual environment as requested: $VENV_PATH"
    elif [[ ! -d "$VENV_PATH" ]]; then
        log "No virtual environment to clean up"
    else
        log "Virtual environment cleanup skipped"
    fi
    
    # Clean up any temporary directories we might have created
    if [[ -n "$SLURM_JOB_ID" ]] && [[ -n "$TMPDIR" ]] && [[ -d "${TMPDIR}/telomere_analysis_$$" ]]; then
        log "Cleaning up temporary work directory..."
        rm -rf "${TMPDIR}/telomere_analysis_$$" 2>/dev/null || true
    fi
    
    log "Cleanup completed. Exiting with code $exit_code"
    exit $exit_code
}

# Set up signal traps for cleanup
trap 'log "Script interrupted by signal"; cleanup_and_exit 130' INT TERM
trap 'log "Script exiting normally"; cleanup_and_exit 0' EXIT

# Parse command line arguments
SKIP_ENV_SETUP=false
FORCE_ENV_REBUILD=false
DRY_RUN=false
VERBOSE=false
KEEP_VENV=false  # Default to always cleaning up virtual environments
SYSTEM_PYTHON=false
MINIMAL_MODULES=false
SKIP_VENV=false
TEST_ENV_ONLY=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -i|--input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -t|--threads)
            THREADS="$2"
            shift 2
            ;;
        -c|--chunk-size)
            CHUNK_SIZE="$2"
            shift 2
            ;;
        -l|--min-length)
            MIN_TELOMERE_LENGTH="$2"
            shift 2
            ;;
        -e|--error-rate)
            MAX_ERROR_RATE="$2"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="$2"
            shift 2
            ;;
        --python-version)
            PYTHON_VERSION="$2"
            shift 2
            ;;
        --venv-name)
            VENV_NAME="$2"
            shift 2
            ;;
        --skip-env-setup)
            SKIP_ENV_SETUP=true
            shift
            ;;
        --force-env-rebuild)
            FORCE_ENV_REBUILD=true
            shift
            ;;
        --system-python)
            SYSTEM_PYTHON=true
            shift
            ;;
        --minimal-modules)
            MINIMAL_MODULES=true
            shift
            ;;
        --keep-venv)
            KEEP_VENV=true
            shift
            ;;
        --cleanup-venv)
            KEEP_VENV=false
            shift
            ;;
        --dry-run)
            DRY_RUN=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --test-env)
            TEST_ENV_ONLY=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            error_exit "Unknown option: $1"
            ;;
    esac
done

# Validate required arguments
if [[ -z "$INPUT_DIR" ]]; then
    error_exit "Input directory is required. Use -i or --input-dir"
fi

if [[ -z "$OUTPUT_DIR" ]]; then
    error_exit "Output directory is required. Use -o or --output-dir"
fi

# Convert to absolute paths
INPUT_DIR=$(realpath "$INPUT_DIR" 2>/dev/null) || error_exit "Invalid input directory: $INPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR" 2>/dev/null) || OUTPUT_DIR=$(mkdir -p "$OUTPUT_DIR" && realpath "$OUTPUT_DIR")

# Validate input directory exists and contains BAM files
if [[ ! -d "$INPUT_DIR" ]]; then
    error_exit "Input directory does not exist: $INPUT_DIR"
fi

# Check for BAM files
BAM_COUNT=$(find "$INPUT_DIR" -name "*.bam" | wc -l)
if [[ $BAM_COUNT -eq 0 ]]; then
    error_exit "No BAM files found in input directory: $INPUT_DIR"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR" || error_exit "Cannot create output directory: $OUTPUT_DIR"

# Determine virtual environment path based on environment
if [[ -n "$SLURM_JOB_ID" ]]; then
    # In SLURM, use a temporary directory that's writable
    VENV_BASE_DIR="${TMPDIR:-/tmp}/slurm_venvs_${USER}"
    mkdir -p "$VENV_BASE_DIR"
    VENV_PATH="${VENV_BASE_DIR}/${VENV_NAME}_${SLURM_JOB_ID}"
    log "SLURM environment detected (Job ID: $SLURM_JOB_ID)"
    log "SLURM_SUBMIT_DIR: ${SLURM_SUBMIT_DIR:-'not set'}"
    log "TMPDIR: ${TMPDIR:-'not set'}"
else
    # Not in SLURM, use script directory
    VENV_PATH="${SCRIPT_DIR}/${VENV_NAME}"
    log "Non-SLURM environment detected"
fi

log "=== Telomere Analysis HPC Wrapper ==="
log "Script directory: $SCRIPT_DIR"
log "Script path: $SCRIPT_PATH"
log "Requirements path: $REQUIREMENTS_PATH"
log "Virtual environment path: $VENV_PATH"
log "Input directory: $INPUT_DIR"
log "Output directory: $OUTPUT_DIR"
log "BAM files found: $BAM_COUNT"
log "Threads: $THREADS"
log "Chunk size: $CHUNK_SIZE"

# Verify script exists before proceeding
if [[ ! -f "$SCRIPT_PATH" ]]; then
    log "ERROR: Python script not found at: $SCRIPT_PATH"
    log "Current working directory: $(pwd)"
    log "Contents of script directory ($SCRIPT_DIR):"
    ls -la "$SCRIPT_DIR" 2>/dev/null || log "Could not list script directory"
    error_exit "Python script not found: $SCRIPT_PATH"
fi
log "✓ Python script found: $SCRIPT_PATH"

# Function to execute or show commands
execute_cmd() {
    if [[ "$VERBOSE" == "true" ]]; then
        log "Executing: $1"
    fi
    
    if [[ "$DRY_RUN" == "true" ]]; then
        echo "[DRY-RUN] $1"
    else
        eval "$1"
    fi
}

# Environment Setup
if [[ "$SKIP_ENV_SETUP" == "false" ]]; then
    log "Setting up Python environment..."
    
    # Force system Python if requested
    if [[ "$SYSTEM_PYTHON" == "true" ]]; then
        log "Using system Python as requested (skipping virtual environment)"
        SKIP_VENV=true
    fi
    
    # Check if running on common HPC systems and load modules
    if command -v module &> /dev/null; then
        if [[ "$MINIMAL_MODULES" == "true" ]]; then
            log "Loading minimal modules for telomere analysis..."
            
            # Only load Python - most minimal approach
            for python_module in "Python/3.11.3" "Python/3.9" "python/3.9"; do
                if execute_cmd "ml load $python_module" 2>/dev/null; then
                    log "✓ Python module loaded: $python_module"
                    break
                fi
            done
        else
            log "Loading HPC modules for telomere analysis..."
            
            # Load only essential module paths
            log "Loading core module paths..."
            execute_cmd "ml use /hpcapps/lib-edda/modules/all/Core" 2>/dev/null || true
            execute_cmd "ml use /hpcapps/lib-bio/modules/all" 2>/dev/null || true
            
            # Load only essential modules for telomere analysis
            log "Loading essential modules..."
            
            # Python is the primary requirement
            PYTHON_LOADED=false
            for python_module in "Python/3.11.3" "Python/3.9" "python/3.9"; do
                if execute_cmd "ml load $python_module" 2>/dev/null; then
                    log "Successfully loaded Python module: $python_module"
                    PYTHON_LOADED=true
                    break
                fi
            done
            
            if [[ "$PYTHON_LOADED" == "false" ]]; then
                log "Warning: No Python module loaded, will try system Python"
            fi
            
            # HTSlib/SAMtools - needed for pysam BAM processing
            for hts_module in "HTSlib" "SAMtools"; do
                execute_cmd "ml load $hts_module" 2>/dev/null || log "Note: $hts_module module not available"
            done
            
            # GCC - only if needed for compilation
            execute_cmd "ml load GCC" 2>/dev/null || log "Note: GCC module not available (may not be needed)"
        fi
        
        # Show loaded modules for debugging
        if [[ "$VERBOSE" == "true" ]]; then
            log "Currently loaded modules:"
            execute_cmd "module list" 2>&1 || true
        fi
    fi
    
    # Find Python executable
    PYTHON_CMD=""
    for py_cmd in "python3" "python$PYTHON_VERSION" "python"; do
        if command -v "$py_cmd" &> /dev/null; then
            # Test if this Python can create virtual environments
            if "$py_cmd" -m venv --help >/dev/null 2>&1; then
                PYTHON_CMD="$py_cmd"
                log "Found working Python: $(which $py_cmd)"
                break
            else
                log "Python $py_cmd found but venv module not available"
            fi
        fi
    done
    
    if [[ -z "$PYTHON_CMD" ]]; then
        log "No Python with venv support found. Trying alternative approaches..."
        
        # Try using virtualenv instead of venv
        if command -v virtualenv &> /dev/null && command -v python3 &> /dev/null; then
            PYTHON_CMD="python3"
            USE_VIRTUALENV=true
            log "Will use virtualenv instead of venv"
        else
            error_exit "No suitable Python found. Please load a Python module or install Python with venv support."
        fi
    fi
    
    log "Using Python: $(which $PYTHON_CMD) ($(${PYTHON_CMD} --version 2>&1))"
    
    # Virtual environment path already determined above based on environment
    log "Virtual environment will be created at: $VENV_PATH"
    
    if [[ "$FORCE_ENV_REBUILD" == "true" ]] && [[ -d "$VENV_PATH" ]]; then
        log "Removing existing virtual environment..."
        execute_cmd "rm -rf '$VENV_PATH'"
    fi
    
    if [[ ! -d "$VENV_PATH" ]]; then
        log "Creating virtual environment: $VENV_PATH"
        
        if [[ "${USE_VIRTUALENV:-false}" == "true" ]]; then
            # Use virtualenv instead of venv
            execute_cmd "virtualenv '$VENV_PATH'" || error_exit "Failed to create virtual environment with virtualenv"
        else
            # Use standard venv
            execute_cmd "$PYTHON_CMD -m venv '$VENV_PATH'" || {
                log "venv creation failed, trying with --system-site-packages..."
                execute_cmd "$PYTHON_CMD -m venv --system-site-packages '$VENV_PATH'" || {
                    log "Standard venv failed, trying pip install --user approach..."
                    SKIP_VENV=true
                }
            }
        fi
    fi
    
    # Activate virtual environment if it was created
    if [[ "$SKIP_VENV" != "true" ]] && [[ -d "$VENV_PATH" ]]; then
        log "Activating virtual environment..."
        source "$VENV_PATH/bin/activate" || {
            log "Failed to activate virtual environment, proceeding with system Python..."
            SKIP_VENV=true
        }
    fi
    
    # Check if we're in a virtual environment or using system Python
    if [[ -n "$VIRTUAL_ENV" ]]; then
        log "Using virtual environment: $VIRTUAL_ENV"
    else
        log "Using system Python installation"
    fi
    
    # Upgrade pip if possible
    log "Upgrading pip..."
    if [[ "$SKIP_VENV" == "true" ]]; then
        execute_cmd "pip install --user --upgrade pip" || log "Warning: Could not upgrade pip"
    else
        execute_cmd "pip install --upgrade pip" || log "Warning: Could not upgrade pip"
    fi
    
    # Install requirements
    if [[ -f "$REQUIREMENTS_PATH" ]]; then
        log "Installing Python dependencies from requirements.txt..."
        if [[ "$SKIP_VENV" == "true" ]]; then
            execute_cmd "pip install --user -r '$REQUIREMENTS_PATH'" || {
                log "Requirements installation failed, trying individual packages..."
                execute_cmd "pip install --user pysam>=0.22.0" || log "Warning: Could not install pysam"
            }
        else
            execute_cmd "pip install -r '$REQUIREMENTS_PATH'" || {
                log "Requirements installation failed, trying individual packages..."
                execute_cmd "pip install pysam>=0.22.0" || log "Warning: Could not install pysam"
            }
        fi
    else
        log "Installing Python dependencies manually..."
        if [[ "$SKIP_VENV" == "true" ]]; then
            execute_cmd "pip install --user pysam>=0.22.0" || log "Warning: Could not install pysam"
        else
            execute_cmd "pip install pysam>=0.22.0" || log "Warning: Could not install pysam"
        fi
    fi
    
    # Verify pysam installation
    log "Verifying telomere analysis dependencies..."
    
    # Check if Python can import required modules
    if python -c "import pysam; print(f'✓ pysam version: {pysam.__version__}')" 2>/dev/null; then
        log "✓ pysam installation verified successfully!"
    else
        log "⚠ pysam installation verification failed, attempting installation..."
        
        # Try different installation approaches
        for install_cmd in \
            "pip install pysam --force-reinstall --no-cache-dir" \
            "pip install --user pysam --force-reinstall --no-cache-dir" \
            "pip install pysam --no-deps" \
            "pip install --user pysam --no-deps"; do
            
            log "Trying: $install_cmd"
            if execute_cmd "$install_cmd" && python -c "import pysam" 2>/dev/null; then
                log "✓ pysam installation successful with: $install_cmd"
                break
            fi
        done
        
        # Final check
        if ! python -c "import pysam" 2>/dev/null; then
            log "✗ Warning: pysam installation failed. The script may not work properly."
            log "  Telomere analysis requires pysam for BAM file processing."
            log "  You may need to:"
            log "    1. Load HTSlib module: ml load HTSlib"
            log "    2. Load development tools: ml load GCC"
            log "    3. Install pysam manually with system package manager"
        fi
    fi
    
    # Check other basic Python modules
    if python -c "import concurrent.futures, multiprocessing, csv, os, sys, subprocess, glob, re, time, argparse" 2>/dev/null; then
        log "✓ All standard Python modules available"
    else
        log "⚠ Some standard Python modules missing (unusual for Python 3.x)"
    fi
    
    log "Environment setup complete!"
    
    # If this is just a test run, exit here
    if [[ "$TEST_ENV_ONLY" == "true" ]]; then
        log "=== Environment Test Complete ==="
        log "✓ Python executable: $(which python)"
        log "✓ Virtual environment: ${VIRTUAL_ENV:-'System Python'}"
        log "✓ Script found: $SCRIPT_PATH"
        log "✓ All dependencies appear to be working"
        log "Test mode complete. Use without --test-env to run actual analysis."
        # Disable EXIT trap and cleanup normally
        trap - EXIT
        cleanup_and_exit 0
    fi
else
    log "Skipping environment setup as requested..."
fi

# Script existence already verified above

# Prepare the Python command
PYTHON_CMD_ARGS=(
    "$SCRIPT_PATH"
    "--input-dir" "$INPUT_DIR"
    "--output-dir" "$OUTPUT_DIR"
    "--threads" "$THREADS"
    "--chunk-size" "$CHUNK_SIZE"
    "--min-telomere-length" "$MIN_TELOMERE_LENGTH"
    "--max-error-rate" "$MAX_ERROR_RATE"
    "--prefix" "$PREFIX"
)

FULL_PYTHON_CMD="python ${PYTHON_CMD_ARGS[*]}"

# Create a log file for the analysis
LOG_FILE="${OUTPUT_DIR}/telomere_analysis.log"

log "Starting telomere analysis..."
log "Working directory: $(pwd)"
log "Python script: $SCRIPT_PATH"
log "Command: $FULL_PYTHON_CMD"
log "Log file: $LOG_FILE"

# Run the analysis
if [[ "$DRY_RUN" == "true" ]]; then
    echo "[DRY-RUN] $FULL_PYTHON_CMD > '$LOG_FILE' 2>&1"
else
    # Create a timestamp for the run
    echo "=== Telomere Analysis Started: $(date) ===" > "$LOG_FILE"
    echo "Command: $FULL_PYTHON_CMD" >> "$LOG_FILE"
    echo "Working directory: $(pwd)" >> "$LOG_FILE"
    echo "Environment: $(env | grep -E '^(SLURM|PBS|SGE|LSB)' || echo 'No job scheduler detected')" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
    
    # Run the analysis and capture output
    if eval "$FULL_PYTHON_CMD" >> "$LOG_FILE" 2>&1; then
        log "Analysis completed successfully!"
        log "Results saved to: $OUTPUT_DIR"
        log "Log file: $LOG_FILE"
        
        # Show quick summary
        if [[ -f "${OUTPUT_DIR}/${PREFIX}.csv" ]]; then
            log "Quick summary:"
            echo "Files processed: $(tail -n +2 "${OUTPUT_DIR}/${PREFIX}.csv" | wc -l)"
            echo "Detailed results: ${OUTPUT_DIR}/${PREFIX}_detailed.csv"
        fi
        
        # Disable EXIT trap and cleanup normally
        trap - EXIT
        cleanup_and_exit 0
    else
        error_exit "Analysis failed! Check log file: $LOG_FILE"
    fi
fi
