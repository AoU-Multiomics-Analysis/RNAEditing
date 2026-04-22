#!/usr/bin/env python

import sys
import numpy as np
from scipy.stats import rankdata, norm
from optparse import OptionParser

def inverse_normal_transform(x):
    """
    Apply inverse normal transformation to a vector.
    
    Parameters:
    x: array of values (editing levels for one site across all samples)
    
    Returns:
    array of normalized values following standard normal distribution
    """
    n = len(x)
    a = 0.5
    
    # Rank the data, convert to quantiles, then to normal distribution
    ranks = rankdata(x)
    quantiles = (ranks - a) / (n + 1.0 - 2.0*a)
    normalized = norm.ppf(quantiles)
    
    return normalized

def main(input_file, output_file):
    
    sys.stderr.write("Starting processing...\n")
    sys.stderr.write(f"Input file: {input_file}\n")
    sys.stderr.write(f"Output file: {output_file}\n\n")
    
    # Read the input file
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
    except:
        sys.stderr.write(f"Error: Cannot open input file {input_file}\n")
        return
    
    # Parse header
    header = lines[0].strip().split()
    sample_names = header[1:]  # Skip 'chrom' column
    num_samples = len(sample_names)
    
    sys.stderr.write(f"Found {num_samples} samples\n")
    
    # Storage for sites that pass filters
    sites_data = []
    sites_coords = []
    
    # Statistics
    total_sites = 0
    filtered_sex_chr = 0
    filtered_variance = 0
    sites_with_imputation = 0
    
    # Process each site
    sys.stderr.write("Processing sites...\n")
    for line in lines[1:]:
        total_sites += 1
        
        if total_sites % 1000 == 0:
            sys.stderr.write(f"  Processed {total_sites} sites...\n")
        
        fields = line.strip().split()
        chrom_full = fields[0]  # e.g., chr1:10185:10186:ENSG12234.1
        
        # Parse chromosome
        chr_parts = chrom_full.replace("chr", "").split(":")
        chromosome = chr_parts[0]
        
        # Filter 1: Remove sex chromosomes
        if chromosome in ['X', 'Y']:
            filtered_sex_chr += 1
            continue
        
        # Parse editing ratios for each sample
        ratios = fields[1:]
        editing_levels = []
        missing_indices = []
        
        for i, ratio_str in enumerate(ratios):
            num, denom = ratio_str.split('/')
            num, denom = float(num), float(denom)
            
            # Mark missing data (0/0)
            if denom < 1:
                editing_levels.append(np.nan)
                missing_indices.append(i)
            else:
                # Add pseudocount
                editing_level = (num + 0.5) / (denom + 0.5)
                editing_levels.append(editing_level)
        
        editing_levels = np.array(editing_levels)
        
        # Imputation: Replace missing values with mean of non-missing values
        num_missing = np.sum(np.isnan(editing_levels))
        if num_missing > 0:
            mean_editing = np.nanmean(editing_levels)
            editing_levels[np.isnan(editing_levels)] = mean_editing
            sites_with_imputation += 1
        
        # Filter 2: Remove low variance sites
        site_std = np.std(editing_levels)
        site_mean = np.mean(editing_levels)
        
        # Filter sites with very low mean (essentially no editing)
        if site_mean == 0:
            filtered_variance += 1
            continue
        
        # Calculate coefficient of variation
        coefficient_of_variation = site_std / site_mean
        
        # Filter sites with low coefficient of variation
        if coefficient_of_variation < 0.8:
            filtered_variance += 1
            continue
        
        # Apply inverse normal transformation (per-site)
        normalized_levels = inverse_normal_transform(editing_levels)
        
        # Store the data
        sites_coords.append(chrom_full)
        sites_data.append(normalized_levels)
    
    sys.stderr.write(f"\nFinished processing {total_sites} sites\n\n")
    
    # Print filtering statistics
    sys.stderr.write("Filtering Summary:\n")
    sys.stderr.write(f"  Total sites input: {total_sites}\n")
    sys.stderr.write(f"  Filtered (sex chromosomes): {filtered_sex_chr}\n")
    sys.stderr.write(f"  Filtered (low variance/CV): {filtered_variance}\n")
    sys.stderr.write(f"  Sites with imputed values: {sites_with_imputation}\n")
    sys.stderr.write(f"  Sites retained: {len(sites_data)}\n")
    sys.stderr.write(f"  Retention rate: {100*len(sites_data)/total_sites:.1f}%\n\n")
    
    # Convert to matrix and sort by genomic position
    sys.stderr.write("Sorting sites by genomic position...\n")
    
    # Create list of (chr, start, coordinate_string, data) for sorting
    sites_for_sorting = []
    for coord, data in zip(sites_coords, sites_data):
        parts = coord.replace("chr", "").split(":")
        chr_num = int(parts[0])
        start_pos = int(parts[1])
        sites_for_sorting.append((chr_num, start_pos, coord, data))
    
    # Sort by chromosome, then position
    sites_for_sorting.sort(key=lambda x: (x[0], x[1]))
    
    # Write output in BED format
    sys.stderr.write(f"Writing output to {output_file}...\n")
    
    with open(output_file, 'w') as out:
        # Write header
        out.write("\t".join(["#Chr", "start", "end", "ID"] + sample_names) + '\n')
        
        # Write data
        for chr_num, start_pos, coord, data in sites_for_sorting:
            parts = coord.replace("chr", "").split(":")
            chromosome = parts[0]
            start = parts[1]
            end = parts[2]
            
            # Extract gene_id if present (format: chr:start:end:gene_id)
            if len(parts) == 4:
                gene_id = parts[3]
                site_id = f"{chromosome}:{start}:{end}_{gene_id}"
            else:
                site_id = f"{chromosome}:{start}:{end}"
            
            # Format: chr, start, end, ID, sample1_value, sample2_value, ...
            data_strings = [f"{val:.6f}" for val in data]
            out.write("\t".join([chromosome, start, end, site_id] + data_strings) + '\n')
    
    sys.stderr.write("\nDone!\n")
    sys.stderr.write(f"Output saved to: {output_file}\n")

if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog -i INPUT_FILE -o OUTPUT_FILE")
    parser.add_option("-i", "--input", dest="input_file", 
                      help="Input matrix file with editing ratios", metavar="FILE")
    parser.add_option("-o", "--output", dest="output_file",
                      help="Output normalized BED file", metavar="FILE")
    
    (options, args) = parser.parse_args()
    
    # Check that both required arguments are provided
    if not options.input_file or not options.output_file:
        sys.stderr.write("Error: Both input and output files must be specified\n\n")
        sys.stderr.write("Usage:\n")
        sys.stderr.write("  python transform.py -i INPUT_FILE -o OUTPUT_FILE\n\n")
        sys.stderr.write("Options:\n")
        sys.stderr.write("  -i, --input   Input matrix file with editing ratios\n")
        sys.stderr.write("  -o, --output  Output normalized BED file\n\n")
        sys.stderr.write("Description:\n")
        sys.stderr.write("  Processes RNA editing matrix for QTL analysis:\n")
        sys.stderr.write("    1. Filters sex chromosomes (X, Y)\n")
        sys.stderr.write("    2. Filters low CV sites (CV < 0.8)\n")
        sys.stderr.write("    3. Imputes missing values (0/0) with mean editing level\n")
        sys.stderr.write("    4. Applies inverse normal transformation per-site\n")
        sys.stderr.write("    5. Outputs sorted BED format file with gene IDs\n\n")
        sys.stderr.write("Example:\n")
        sys.stderr.write("  python transform.py -i editing_matrix.txt -o normalized_output.bed\n")
        parser.print_help()
        exit(1)
    
    main(options.input_file, options.output_file)