#!/usr/bin/env python

import sys
import csv
import gzip
import numpy as np
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
    """Median absolute deviation, ignoring NaNs."""
    valid = values[~np.isnan(values)]
    if valid.size == 0:
        return np.nan
    med = np.median(valid)
    mad = np.median(np.abs(valid - med))
    return mad


def main(input_file, output_file, window_size):

    MAD_THRESHOLD = 0.003

    sys.stderr.write("Starting processing...\n")
    sys.stderr.write(f"Input file: {input_file}\n")
    sys.stderr.write(f"Output file: {output_file}\n")
    sys.stderr.write(f"Window size: ±{window_size:,} bp\n\n")

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

    # retained sites for transformed BED output
    sites_data = []
    sites_coords = []

    # metadata for ALL sites (including filtered)
    # rows: site, mean_editing, mad, n_missing, kept(boolean)
    site_metadata = []

    total_sites = 0
    filtered_low_mad = 0

    sys.stderr.write("Processing sites...\n")

    for line in lines[1:]:
        total_sites += 1
        if total_sites % 1000 == 0:
            sys.stderr.write(f"  Processed {total_sites} sites...\n")

        fields = line.strip().split()
        if not fields:
            continue

        chrom_full = fields[0]
        ratios = fields[1:]

        editing_levels = []
        for ratio_str in ratios:
            num, denom = ratio_str.split('/')
            num, denom = float(num), float(denom)

            if denom < 1:
                editing_levels.append(np.nan)
            else:
                editing_levels.append((num + 0.5) / (denom + 0.5))

        editing_levels = np.array(editing_levels, dtype=float)

        n_missing = int(np.sum(np.isnan(editing_levels)))
        mean_edit = float(np.nanmean(editing_levels)) if n_missing < len(editing_levels) else np.nan
        mad = compute_mad(editing_levels)

        # decide whether we keep this site for output matrix
        keep = True
        if (not np.isfinite(mad)) or (mad < MAD_THRESHOLD):
            keep = False
            filtered_low_mad += 1

        # record metadata for ALL sites (even filtered)
        site_metadata.append((chrom_full, mean_edit, mad, n_missing, keep))

        if not keep:
            continue

        mask = ~np.isnan(editing_levels)
        obs = editing_levels[mask]
        if obs.size == 0:
            # should be rare given MAD filter, but safe
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

    sys.stderr.write(f"Writing output to {output_file}...\n")

    with open_maybe_gzip(output_file, "wt") as out:
        out.write("\t".join(["#Chr", "start", "end", "ID"] + sample_names) + "\n")

        for chr_key, start_pos, coord, data in sites_for_sorting:
            parts = coord.replace("chr", "").split(":")
            chromosome = parts[0]
            original_start = int(parts[1])
            original_end = int(parts[2])

            # Calculate windowed coordinates
            # Use the start position (editing site) as the center
            windowed_start = max(0, original_start - window_size)
            windowed_end = original_start + window_size

            # Build ID from original coordinates
            if len(parts) == 4:
                gene_id = parts[3]
                site_id = f"{chromosome}:{original_start}:{original_end}_{gene_id}"
            else:
                site_id = f"{chromosome}:{original_start}:{original_end}"

            # keep missing as NA
            data_strings = ["NA" if np.isnan(val) else f"{val:.6f}" for val in data]
            out.write("\t".join([chromosome, str(windowed_start), str(windowed_end), site_id] + data_strings) + "\n")

    sys.stderr.write("Writing site metadata to site_metadata.csv...\n")
    with open("site_metadata.csv", "w", newline="") as csvfile:
        w = csv.writer(csvfile)
        w.writerow(["site", "mean_editing", "mad", "n_missing", "kept"])
        for site, mean_edit, mad_val, n_missing, kept in site_metadata:
            # write missing mean/mad as NA for readability
            mean_out = "NA" if (mean_edit is None or not np.isfinite(mean_edit)) else mean_edit
            mad_out = "NA" if (mad_val is None or not np.isfinite(mad_val)) else mad_val
            w.writerow([site, mean_out, mad_out, n_missing, int(kept)])

    sys.stderr.write("\nDone!\n")
    sys.stderr.write(f"Output saved to: {output_file}\n")
    sys.stderr.write("Metadata saved to: site_metadata.csv\n")


if __name__ == "__main__":
    parser = OptionParser(usage="usage: %prog -i INPUT_FILE -o OUTPUT_FILE -w WINDOW_SIZE")
    parser.add_option("-i", "--input", dest="input_file",
                      help="Input file with editing ratios")
    parser.add_option("-o", "--output", dest="output_file",
                      help="Output BED file (will be gzipped if .gz extension)")
    parser.add_option("-w", "--window", dest="window_size", type="int", default=100000,
                      help="Window size in bp to extend in each direction from editing site (default: 100000)")

    (options, args) = parser.parse_args()

    if not options.input_file or not options.output_file:
        sys.stderr.write("Error: must provide input and output files\n")
        parser.print_help()
        sys.exit(1)

    main(options.input_file, options.output_file, options.window_size)