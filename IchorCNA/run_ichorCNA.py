"""
IchorCNA WSL Analysis Runner
Run IchorCNA analysis using WSL environment
"""

import os
import sys
from pathlib import Path

# Add the current directory to Python path
sys.path.append(str(Path(__file__).parent))

from ichorcna_wsl_pipeline import IchorCNAPipelineWSL

def run_ichorcna_wsl_analysis():
    """Run IchorCNA analysis using WSL."""
    
    # Configuration - WSL paths
    INPUT_DIR = "/mnt/d/BRCA2-misc_files/hg38_1-24"  # WSL path to your BAM files
    OUTPUT_DIR = "/home/gissu/ichorCNA"  # WSL output directory
    NORMAL_RESTARTS = "c(0.5,0.8)"  # Lower-sensitivity default (alternatives: c(0.8), c(0.5))
    
    print("🧬 IchorCNA Nanopore Analysis Pipeline (WSL)")
    print("=" * 60)
    print(f"Input directory: {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"ichorCNA normal restarts: {NORMAL_RESTARTS}")
    print("Using WSL Ubuntu environment")
    print()
    
    # Initialize pipeline
    pipeline = IchorCNAPipelineWSL(INPUT_DIR, OUTPUT_DIR, normal_restarts=NORMAL_RESTARTS)
    
    # Run the analysis
    print("Starting IchorCNA pipeline in WSL...")
    success = pipeline.run_pipeline()
    
    if success:
        print("\\n✅ Analysis completed successfully!")
        print(f"Results are available in: {OUTPUT_DIR}")
        print("\\nOutput directories:")
        print(f"  - Filtered BAMs: {OUTPUT_DIR}/filtered_bams/")
        print(f"  - Read counts: {OUTPUT_DIR}/readcounts/")
        print(f"  - IchorCNA results: {OUTPUT_DIR}/ichorCNA_results/")
        print(f"  - Logs: {OUTPUT_DIR}/logs/")
        print("\\nTo access results from Windows:")
        print(f"  \\\\wsl.localhost\\Ubuntu{OUTPUT_DIR}")
    else:
        print("\\n❌ Analysis failed. Check the logs for details.")
    
    return success

if __name__ == "__main__":
    run_ichorcna_wsl_analysis()
