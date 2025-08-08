#!/usr/bin/env python3
"""
Convert BED Files to 450K Array Format for MethAtlas Package

This script transforms methylation data from BED files into the 450K array CSV format
using the BED manifest approach, which has been shown to provide superior mapping accuracy.
The output is compatible with the MethAtlas package.
"""

import os
import sys
import glob
import csv
import time
import bisect
import traceback
from collections import Counter
import pandas as pd
import numpy as np

# Configure paths
BED_MANIFEST_FILE = "path/to/your/manifest"
BED_DIR = "/path/to/your/bed/file/directory"
OUTPUT_DIR = "/path/to/your/output/directory"
LOG_FILE = os.path.join(OUTPUT_DIR, "mapping_log.txt")

# Important parameters
MAX_DISTANCE = 20  # Maximum distance in bp to consider a match
BETA_VALUE_FIELD = 11  # Field index for beta value in BED file
USE_BEST_MATCH = True  # If True, use the closest position within MAX_DISTANCE

# Create output directory
try:
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Created output directory: {OUTPUT_DIR}")
except Exception as e:
    print(f"Error creating output directory: {e}")
    sys.exit(1)

def log(message):
    """Log message to both console and log file"""
    timestamp = time.strftime('%Y-%m-%d %H:%M:%S')
    full_message = f"[{timestamp}] {message}"
    
    # Print to console
    print(full_message, flush=True)
    
    # Write to log file
    try:
        with open(LOG_FILE, 'a') as f:
            f.write(full_message + '\n')
    except Exception as e:
        print(f"ERROR writing to log file: {e}", flush=True)

def normalize_chromosome(chrom):
    """Normalize chromosome names to handle different formats"""
    # Remove 'chr' prefix if present
    norm_chrom = str(chrom).replace('chr', '').upper()
    
    # Remove any alt/random suffixes
    if '_' in norm_chrom:
        norm_chrom = norm_chrom.split('_')[0]
    
    # Handle special chromosomes
    if norm_chrom == 'X':
        return 'X'
    elif norm_chrom == 'Y':
        return 'Y'
    elif norm_chrom in ['M', 'MT']:
        return 'MT'
    
    # Try to convert to integer for numeric chromosomes
    try:
        return str(int(norm_chrom))
    except ValueError:
        # If it's an alternative contig/scaffold that can't be parsed, return None
        if any(x in chrom.upper() for x in ['_ALT', '_RANDOM', 'HSCHR', 'UN_', 'UNLOCALIZED']):
            return None
        return norm_chrom

def convert_beta_value(beta_str):
    """Convert beta value string to float (0-1 range)"""
    try:
        # BED files appear to store beta as percentage (0-100)
        beta = float(beta_str)
        # Convert percentage to 0-1 range if needed
        if beta > 1:  # Only convert if it's in percentage format
            beta = beta / 100.0
        return beta
    except (ValueError, TypeError):
        return None

def format_illumina_id(probe_id):
    """Format probe ID to ensure it has the 'cg' prefix required by MethAtlas"""
    # If ID is numeric or doesn't start with cg, add the prefix
    if probe_id.isdigit() or not probe_id.startswith('cg'):
        return f"cg{probe_id}"
    return probe_id

def load_bed_manifest():
    """Load the Illumina 450K manifest from BED format"""
    log(f"Loading BED Illumina 450K manifest: {BED_MANIFEST_FILE}")
    
    manifest_dict = {}  # {chr: {position: probe_id}}
    manifest_count = 0
    
    try:
        # Pre-initialize dictionary for all chromosomes
        for chr_num in range(1, 23):
            manifest_dict[str(chr_num)] = {}
        manifest_dict['X'] = {}
        manifest_dict['Y'] = {}
        manifest_dict['MT'] = {}
        
        with open(BED_MANIFEST_FILE, 'r') as f:
            for line in f:
                if line.startswith('#') or line.startswith('track') or line.strip() == '':
                    continue
                
                fields = line.strip().split('\t')
                # BED format: chrom, start, end, name, score, strand, etc.
                if len(fields) < 4:  # Need at least chrom, start, end, and probe name
                    continue
                
                chr_val = fields[0]  # chromosome
                start_pos = int(fields[1])  # 0-based start position
                probe_id = fields[3]  # probe ID is typically in the name field
                
                # Format the probe ID to ensure it has the cg prefix
                formatted_probe_id = format_illumina_id(probe_id)
                
                # Normalize chromosome name
                chr_norm = normalize_chromosome(chr_val)
                
                if not chr_norm or chr_norm not in manifest_dict:
                    continue
                    
                # Store position
                manifest_dict[chr_norm][start_pos] = formatted_probe_id
                manifest_count += 1
                
                # Log progress periodically
                if manifest_count % 100000 == 0:
                    log(f"Processed {manifest_count:,} BED manifest lines...")
        
        log(f"Successfully loaded {manifest_count:,} positions from BED manifest")
        
        # Count positions per chromosome
        for chr_key, pos_dict in manifest_dict.items():
            if pos_dict:  # If not empty
                log(f"  BED Chr {chr_key}: {len(pos_dict):,} positions")
        
        return manifest_dict
    
    except Exception as e:
        log(f"ERROR loading BED manifest: {e}")
        log(traceback.format_exc())
        return None

def process_bed_file(bed_file, manifest_dict, max_distance=MAX_DISTANCE):
    """Process a BED file and map to 450K positions"""
    sample_name = os.path.basename(bed_file).replace('.bed', '')
    log(f"Processing {sample_name}...")
    
    # Results for this sample
    results = {
        'sample_name': sample_name,
        'mapped_positions': {},  # {probe_id: beta_value}
        'stats': {
            'total_positions': 0,
            'mapped_positions': 0,
            'exact_matches': 0,
            'approximate_matches': 0,
            'chromosomes': Counter(),  # Count of positions by chromosome
            'mapped_by_chr': Counter(),  # Count of mapped positions by chromosome
            'distances': Counter(),
        }
    }
    
    # Pre-compute position ranges for faster proximity mapping
    position_ranges = {}
    for chr_norm in manifest_dict:
        if not manifest_dict[chr_norm]:  # Skip empty chromosomes
            continue
            
        # Create sorted list of positions for binary search
        positions = sorted(manifest_dict[chr_norm].keys())
        position_ranges[chr_norm] = positions
    
    try:
        with open(bed_file, 'r') as f:
            # Process in batches for better performance
            batch_size = 100000
            
            # Process file in chunks
            for chunk_idx, chunk in enumerate(iter(lambda: f.readlines(batch_size), [])):
                for line_idx, line in enumerate(chunk, 1):
                    if line.startswith('#') or line.strip() == '':
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) <= BETA_VALUE_FIELD:
                        continue
                    
                    # Extract fields
                    chr_val = fields[0]
                    start_pos = int(fields[1])  # Already 0-based, no adjustment needed
                    beta_str = fields[BETA_VALUE_FIELD]
                    
                    # Normalize chromosome
                    chr_norm = normalize_chromosome(chr_val)
                    results['stats']['chromosomes'][chr_norm] += 1
                    
                    # Convert beta value
                    beta_value = convert_beta_value(beta_str)
                    if beta_value is None:
                        continue
                    
                    results['stats']['total_positions'] += 1
                    
                    # Skip if chromosome not in manifest
                    if not chr_norm or chr_norm not in manifest_dict:
                        continue
                    
                    # No need to adjust position - already in 0-based format
                    
                    # Try to find exact or approximate match
                    best_match_probe = None
                    best_match_dist = float('inf')
                    
                    # Check for exact match first (O(1) operation)
                    if start_pos in manifest_dict[chr_norm]:
                        best_match_probe = manifest_dict[chr_norm][start_pos]
                        best_match_dist = 0
                        results['stats']['exact_matches'] += 1
                    elif USE_BEST_MATCH and chr_norm in position_ranges:
                        # Use binary search to find closest position efficiently
                        positions = position_ranges[chr_norm]
                        
                        # Binary search for closest position
                        idx = bisect.bisect_left(positions, start_pos)
                        
                        # Check position to the left
                        if idx > 0:
                            left_pos = positions[idx-1]
                            left_dist = abs(start_pos - left_pos)
                            if left_dist <= max_distance and left_dist < best_match_dist:
                                best_match_probe = manifest_dict[chr_norm][left_pos]
                                best_match_dist = left_dist
                                results['stats']['approximate_matches'] += 1
                        
                        # Check position to the right
                        if idx < len(positions):
                            if positions[idx] == start_pos:  # Exact match
                                best_match_probe = manifest_dict[chr_norm][start_pos]
                                best_match_dist = 0
                                results['stats']['exact_matches'] += 1
                            elif idx < len(positions):  # Approximate match
                                right_pos = positions[idx]
                                right_dist = abs(start_pos - right_pos)
                                if right_dist <= max_distance and right_dist < best_match_dist:
                                    best_match_probe = manifest_dict[chr_norm][right_pos]
                                    best_match_dist = right_dist
                                    results['stats']['approximate_matches'] += 1
                    
                    # Store mapping if found
                    if best_match_probe is not None:
                        # If the probe already exists, we'll take the average beta value
                        if best_match_probe in results['mapped_positions']:
                            existing_beta = results['mapped_positions'][best_match_probe]
                            results['mapped_positions'][best_match_probe] = (existing_beta + beta_value) / 2
                        else:
                            results['mapped_positions'][best_match_probe] = beta_value
                            results['stats']['mapped_positions'] += 1
                            results['stats']['distances'][best_match_dist] += 1
                            results['stats']['mapped_by_chr'][chr_norm] += 1
                
                # Status update after each chunk
                if (chunk_idx + 1) % 10 == 0:
                    log(f"  Processed {results['stats']['total_positions']:,} positions...")
                    
        # Calculate mapping rate
        mapping_rate = (results['stats']['mapped_positions'] / results['stats']['total_positions'] * 100) if results['stats']['total_positions'] > 0 else 0
        log(f"  {sample_name}: {results['stats']['mapped_positions']:,} mapped of {results['stats']['total_positions']:,} positions ({mapping_rate:.2f}%)")
        log(f"  Exact matches: {results['stats']['exact_matches']:,}, Approximate matches: {results['stats']['approximate_matches']:,}")
        
        return results
    
    except Exception as e:
        log(f"ERROR processing {sample_name}: {e}")
        log(traceback.format_exc())
        return None

def export_to_csv(results, output_file, methatlas_format=True):
    """Export mapped positions to CSV format compatible with 450K array and MethAtlas"""
    log(f"Exporting data to {output_file}")
    
    try:
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            if methatlas_format:
                # Write header in MethAtlas format
                writer.writerow(['IlmnID', 'Beta_Value'])
            else:
                # Standard format
                writer.writerow(['ProbeID', 'Beta_Value'])
            
            # Write data rows
            for probe_id, beta_value in results['mapped_positions'].items():
                writer.writerow([probe_id, f"{beta_value:.6f}"])
        
        log(f"Successfully exported {len(results['mapped_positions']):,} positions to {output_file}")
        return True
    
    except Exception as e:
        log(f"ERROR exporting data: {e}")
        log(traceback.format_exc())
        return False

def export_combined_methatlas_format(all_results, output_file):
    """Export all results in a combined format suitable for MethAtlas"""
    log(f"Creating combined MethAtlas format output: {output_file}")
    
    try:
        # Create a set of all probe IDs across all samples
        all_probe_ids = set()
        for result in all_results.values():
            all_probe_ids.update(result['mapped_positions'].keys())
        
        log(f"Total unique probe IDs across all samples: {len(all_probe_ids):,}")
        
        # Create a DataFrame with probe IDs as index
        matrix_df = pd.DataFrame(index=sorted(all_probe_ids))
        matrix_df.index.name = 'IlmnID'  # Important for MethAtlas format
        
        # Add each sample as a column
        for sample_name, result in all_results.items():
            # Create a series of beta values
            beta_values = pd.Series(result['mapped_positions'], name=sample_name)
            # Add to matrix
            matrix_df[sample_name] = beta_values
        
        # Fill NaN values with NA (MethAtlas prefers NA over NaN)
        matrix_df = matrix_df.fillna('NA')
        
        # Save to CSV
        matrix_df.to_csv(output_file)
        log(f"Successfully exported combined matrix with {len(matrix_df.index)} rows and {len(matrix_df.columns)} columns")
        return True
    
    except Exception as e:
        log(f"ERROR creating combined output: {e}")
        log(traceback.format_exc())
        return False

def main():
    """Main execution function"""
    start_time = time.time()
    
    # Clear log file
    with open(LOG_FILE, 'w') as f:
        f.write("BED to 450K Array Transformation for MethAtlas\n")
        f.write("===========================================\n\n")
    
    log("Starting BED to 450K array transformation using BED manifest approach")
    log(f"BED Manifest: {BED_MANIFEST_FILE}")
    log(f"BED Directory: {BED_DIR}")
    log(f"Maximum distance allowed: {MAX_DISTANCE}bp")
    log(f"Important: Ensuring all probe IDs have 'cg' prefix for MethAtlas compatibility")
    
    # Load BED manifest
    manifest_dict = load_bed_manifest()
    if not manifest_dict:
        log("Failed to load BED manifest. Exiting.")
        return
    
    # Find all BED files
    bed_files = glob.glob(os.path.join(BED_DIR, "*.bed"))
    if not bed_files:
        log(f"No BED files found in {BED_DIR}")
        return
    
    log(f"Found {len(bed_files)} BED files")
    
    # Process each BED file
    all_results = {}
    
    for file_idx, bed_file in enumerate(bed_files, 1):
        log(f"[{file_idx}/{len(bed_files)}] Processing file: {os.path.basename(bed_file)}...")
        
        # Process with BED manifest
        result = process_bed_file(bed_file, manifest_dict)
        if result:
            all_results[result['sample_name']] = result
            
            # Export individual result in MethAtlas format
            output_file = os.path.join(OUTPUT_DIR, f"{result['sample_name']}_450k.csv")
            export_to_csv(result, output_file, methatlas_format=True)
    
    # Export combined results for MethAtlas
    if all_results:
        # Combined matrix format
        combined_output = os.path.join(OUTPUT_DIR, "methatlas_combined_450k.csv")
        export_combined_methatlas_format(all_results, combined_output)
        
    # Finish
    elapsed_time = time.time() - start_time
    log(f"Transformation completed in {elapsed_time:.2f} seconds")
    log(f"Results saved to: {OUTPUT_DIR}")
    log(f"The files are now ready to be used with the MethAtlas package")

if __name__ == "__main__":
    main()
