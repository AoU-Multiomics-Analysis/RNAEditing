#!/usr/bin/env python

import sys
import numpy as np
import gzip
from scipy.stats import rankdata, norm
from optparse import OptionParser


def open_maybe_gzip(path, mode="rt"):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def ensure_gz(path):
    return path if path.endswith(".gz") else path + ".gz"


def inverse_normal_transform(x):
    n = len(x)
    a = 0.5
    ranks = rankdata(x)
    quantiles = (ranks - a) / (n + 1.0 - 2.0 * a)
    return norm.ppf(quantiles)


def compute_mad(values):
    """
    Median absolute deviation, ignoring NaNs
    """
    valid = values[~np.isnan(values)]
    if len(valid) == 0:
        return None

    med = np.median(valid)
    mad = np.median(np.abs(valid - med))
    return mad


def main(input_file, output_file):

    MAD_THRESHOLD = 0.003

    sys.stderr.write("Starting processing...\n")
    sys.stderr.write(f"Input file: {input_file}\n")
    sys.stderr.write(f"Output file: {output_file}\n\n")

    output_file = ensure_gz(output_file)

    try:
        with open_maybe_gzip(input_file, "rt") as f:
            lines = f.readlines()
    except Exception as e:
        sys.stderr.write(f"Error: Cannot open input file {input_file}: {e}\n")
        return

    header = lines[0].strip().split()
    sample_names = header[1:]
    num_samples = len(sample_names)

    sys.stderr.write(f"Found {num_samples} samples\n")

    sites_data = []
    sites_coords = []

    total_sites = 0
    filtered_low_mad = 0

    sys.stderr.write("Processing sites...\n")

    for line in lines[1:]:
        total_sites += 1

        if total_sites % 1000 == 0:
            sys.stderr.write(f"  Processed {total_sites} sites...\n")

        fields = line.strip().split()
        chrom_full = fields[0]

        chr_parts = chrom_full.replace("chr", "").split(":")
        chromosome = chr_parts[0]

        ratios = fields[1:]
        editing_levels = []

        for ratio_str in ratios:
            num, denom = ratio_str.split('/')
            num, denom = float(num), float(denom)

            if denom < 1:
                editing_levels.append(np.nan)
            else:
                val = (num + 0.5) / (denom + 0.5)
                editing_levels.append(val)

        editing_levels = np.array(editing_levels)

        # Apply MAD filter
        mad = compute_mad(editing_levels)

        if mad is None or not np.isfinite(mad):
            filtered_low_mad += 1
            continue

        if mad < MAD_THRESHOLD:
            filtered_low_mad += 1
            continue

        mask = ~np.isnan(editing_levels)
        obs = editing_levels[mask]
        
        # if nothing observed, skip (should already be handled by MAD, but safe)
        if obs.size == 0:
            continue
        
        norm_obs = inverse_normal_transform(obs)
        
        normalized_levels = np.full(editing_levels.shape, np.nan, dtype=float)
        normalized_levels[mask] = norm_obs
        
        sites_coords.append(chrom_full)
        sites_data.append(normalized_levels)

    sys.stderr.write(f"\nFinished processing {total_sites} sites\n\n")

    sys.stderr.write("Filtering Summary:\n")
    sys.stderr.write(f"  Total sites input: {total_sites}\n")
    sys.stderr.write(f"  Filtered (low MAD < {MAD_THRESHOLD}): {filtered_low_mad}\n")
    sys.stderr.write(f"  Sites retained: {len(sites_data)}\n")
    sys.stderr.write(f"  Retention rate: {100*len(sites_data)/total_sites:.1f}%\n\n")

    # -----------------------------
    # Sorting
    # -----------------------------
    sys.stderr.write("Sorting sites by genomic position...\n")

    sites_for_sorting = []
    for coord, data in zip(sites_coords, sites_data):
        parts = coord.replace("chr", "").split(":")
        chrom = parts[0]
        start_pos = int(parts[1])

        try:
            c = int(chrom)
            chr_key = (0, c) if 1 <= c <= 23 else (1, chrom)
        except ValueError:
            chr_key = (1, chrom)

        sites_for_sorting.append((chr_key, start_pos, coord, data))

    sites_for_sorting.sort(key=lambda x: (x[0], x[1]))

    # -----------------------------
    # Output
    # -----------------------------
    sys.stderr.write(f"Writing output to {output_file}...\n")

    with open_maybe_gzip(output_file, "wt") as out:
        out.write("\t".join(["#Chr", "start", "end", "ID"] + sample_names) + '\n')

        for chr_key, start_pos, coord, data in sites_for_sorting:
            parts = coord.replace("chr", "").split(":")
            chromosome = parts[0]
            start = parts[1]
            end = parts[2]

            if len(parts) == 4:
                gene_id = parts[3]
                site_id = f"{chromosome}:{start}:{end}_{gene_id}"
            else:
                site_id = f"{chromosome}:{start}:{end}"

            data_strings = [f"{val:.6f}" for val in data]
            out.write("\t".join([chromosome, start, end, site_id] + data_strings) + '\n')

    sys.stderr.write("\nDone!\n")
    sys.stderr.write(f"Output saved to: {output_file}\n")


if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog -i INPUT_FILE -o OUTPUT_FILE")
    parser.add_option("-i", "--input", dest="input_file")
    parser.add_option("-o", "--output", dest="output_file")

    (options, args) = parser.parse_args()

    if not options.input_file or not options.output_file:
        sys.stderr.write("Error: must provide input and output\n")
        sys.exit(1)

    main(options.input_file, options.output_file)