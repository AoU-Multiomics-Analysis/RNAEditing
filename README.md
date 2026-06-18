# A-to-I RNA Editing in All of Us  
This repository provides WDL workflows and scripts for quantifying A-to-I RNA editing levels from AoU data and then preparing for downstream QTL calling/fine-mapping. 

## Step 1: Quantification

The step 1 workflow quantifies RNA editing levels at a predefined set of known editing sites for aligned RNA-seq reads in a single BAM file. Note that to analyze multiple samples, you will need to submit a separate job for each sample. This workflow invokes two scripts to perform the quantification: `query_editing_levels.pl` and `parse_pileup_query.pl`. 

### Workflow: `workflows/step1_quantify.wdl`

#### Required inputs: 
Input	| Type | Description
| --- | --- | --- |
bam_file | File | Aligned RNA-seq reads (BAM format, must be indexed)
bam_index	| File | BAM index file (.bai)
reference_genome | File | Reference genome in FASTA format
reference_genome_index | File | Reference genome index (.fai)
sample_id | String | Unique sample identifier (used in output filenames)
memory | Int | Memory allocation in GB 
disk_space | Int | Disk space in GB 
num_threads | Int | Number of CPU threads 

#### Outputs:
| Output | Type | Description |
| --- | --- | --- |
| editing_counts_all | File | Gzipped TSV file named `<sample_id>.rnaediting_op.gz` containing editing measurements for all chromosomes |
| chr_files | Array[File] | 22 separate gzipped files named `<sample_id>.chr1.rnaediting_op.gz` through `<sample_id>.chr22.rnaediting_op.gz` (one per autosomal chromosome), containing editing quantifications only for sites on that chromosome. Note that sex chromosomes are excluded. |

The 'editing_counts_all' file is formatted as follows: 

```tsv
chrom    position    gene_id         coverage    editedreads    editlevel
chr1      1234567     ENSG00000123    45          12             0.267
chr1      1234890     ENSG00000456    120         88             0.733
chr2      9876543     ENSG00000789    0           0              N/A
```

`chrom`: Chromosome name. 

`position`: Genomic position (1-based). 

`gene_id`: Ensembl gene ID.

`coverage`: Number of edited reads + reference reads at that site passing quality filters. Note that it may not be equivalent to the actual total coverage at that site in cases where there are non-reference reads which do not correspond to an edited read. 

`editedreads`: Number of reads showing editing (A→G on if editing site is on + strand, T→C if editing site is on - strand).

`editlevel`: Editing frequency (edited reads / coverage), or N/A if no coverage at that site. 

---

### Script: `scripts/query_editing_levels.pl`

Runs `samtools mpileup` at predefined RNA editing sites in a BAM file and computes per-site editing levels using base-quality-filtered pileup parsing.

**Inputs:**

| Argument | Type | Description |
|----------|------|------------|
| `--bam` | File | Input BAM file (indexed) |
| `--output` | String | Output file prefix |
| `--genome` | File | Reference genome FASTA |
| `--sites` | File | BED file of known RNA editing sites |
| `--minbasequal` | Int | Minimum base quality (default: 20) |
| `--minmapqual` | Int | Minimum mapping quality (default: 250) |
| `--offset` | Int | Base quality offset (default: 33) |
| `--samtools` | String | Path to samtools executable |

**Output:**

- `<sample_id>.rnaediting_op.gz` — gzipped TSV file containing per-site editing measurements (chromosome, position, gene_id, coverage, editedreads, editlevel)


---

### Script: `scripts/parse_pileup_query.pl`

Parses `samtools mpileup` output and counts reference and nucleotide-specific read support after applying base-quality filtering. Note that this script is internally called by `query_editing_levels.pl` and should not need to be externally invoked. 

**Inputs:**

| Argument | Type | Description |
|----------|------|------------|
| `pileup line` | String | Single mpileup record |
| `minbasequal` | Int | Minimum base quality threshold |
| `offset` | Int | Base quality encoding offset |

**Output:**

- Returns a list of nucleotide counts per site:
  - reference count
  - A count
  - T count
  - C count
  - G count
 
## Step 2: Combine Sample Matrices

The step 2 workflow combines per-sample RNA editing quantifications generated in Step 1 into a single matrix suitable for downstream QTL mapping and fine-mapping analyses. Chromosome-specific sample aggregation files are also produced. This workflow invokes one script, `combine_sample_matrices.pl`, to merge editing measurements across samples while applying coverage and sample-presence filters.

### Workflow: `workflows/step2_aggregate.wdl`

#### Required inputs:

| Input | Type | Description |
| --- | --- | --- |
| chr_fofns | Array[File] | Array of file-of-filenames (FOFN) files, one per chromosome. Each FOFN contains paths to the corresponding chromosome-specific editing quantification files generated in Step 1. |
| output_prefix | String | Prefix used for all output matrix filenames. |
| memory | Int | Memory allocation in GB |
| disk_space | Int | Disk space in GB |
| num_threads | Int | Number of CPU threads |

#### Outputs:

| Output | Type | Description |
| --- | --- | --- |
| per_chr_matrices | Array[File] | Chromosome-specific editing matrices named `<output_prefix>.chr1.gz` through `<output_prefix>.chr22.gz` |
| combined_matrix | File | Genome-wide merged matrix named `<output_prefix>.all_chr.gz` |

The chromosome-specific and combined matrices are formatted as follows:

```text
chrom sample_1 sample_2 sample_3 ...
chr1:1234566:1234567:ENSG00000123 12/45 8/41 0/0 ...
chr1:1234889:1234890:ENSG00000456 88/120 74/110 91/125 ...
chr2:9876542:9876543:ENSG00000789 0/0 15/52 22/61 ...
```

The first column contains a unique site identifier in the format:

```text
chrom:start:end:gene_id
```

where:

- `chrom`: Chromosome name.
- `start`: 0-based genomic start coordinate.
- `end`: 1-based genomic end coordinate.
- `gene_id`: Ensembl gene identifier associated with the editing site. If site falls within multiple gene annotations, their Ensembl gene identifiers are concatenated with a semicolon. 

Each subsequent column corresponds to a sample and contains editing counts in the format:

```text
edited_reads/coverage
```

where:

- `edited_reads`: Number of reads supporting editing at the site.
- `coverage`: Total filtered coverage at the site.

Sites not meeting coverage thresholds in a given sample are represented as:

```text
0/0
```

Only sites passing the minimum sample-presence threshold are retained in the final matrix.

---

### Script: `scripts/combine_sample_matrices.pl`

Combines per-sample RNA editing quantification files into a single sample-by-site matrix. The script applies coverage filtering, retains sites observed in a user-specified fraction of samples, and outputs a compressed matrix suitable for downstream association analyses.

**Inputs:**

| Argument | Type | Description |
|----------|------|------------|
| `--input`, `-i` | Directory | Directory containing per-sample `.rnaediting_op` or `.rnaediting_op.gz` files |
| `--output`, `-o` | String | Output matrix filename |
| `--mincov`, `-c` | Int | Minimum coverage required for a site within a sample to be included (default: 20) |
| `--minsamps`, `-s` | Float | Minimum fraction of samples required to pass coverage filtering at a site for the site to be retained (default: 0.9) |
| `--help`, `-h` | Flag | Display usage information |

**Output:**

- `<output>.gz` — gzipped matrix containing editing counts for all retained sites across all samples.

**Filtering behavior:**

1. A site is considered present in a sample only if:

   ```text
   coverage >= mincov
   ```

2. A site is retained in the final matrix only if it is present in at least these many samples:

   ```text
   ceil(minsamps × number_of_samples)
   ```

3. Samples failing the coverage threshold for a retained site are represented as:

   ```text
   0/0
   ```

**Output format:**

```text
chrom sample_1 sample_2 sample_3 ...
chr1:1234566:1234567:ENSG00000123 12/45 8/41 0/0 ...
```

where each sample column contains:

```text
edited_reads/coverage
```

The sample columns are sorted alphabetically to ensure consistent ordering across chromosome-specific matrices.

## Step 3: Normalization and Principal Component Analysis

The step 3 workflow converts editing count matrices into normalized phenotype matrices suitable for QTL mapping. This workflow first filters low-variance editing sites, imputes missing values, and applies inverse normal transformation using `transform.py`. It then computes phenotype principal components (PCs) using `calculate_pcs.R` for use as covariates in downstream association analyses.

### Workflow: `workflows/step3_normalize.wdl`

#### Required inputs:

| Input | Type | Description |
| --- | --- | --- |
| matrix_file | File | Combined editing matrix generated in Step 2 |
| sample_list | File | TSV file containing samples to retain. Must contain a `research_id` column matching sample IDs in the matrix. |
| output_prefix | String | Prefix used for output filenames |
| mad_threshold | Float | Minimum median absolute deviation (MAD) required for a site to be retained |
| memory | Int | Memory allocation in GB |
| disk_space | Int | Disk space in GB |
| num_threads | Int | Number of CPU threads |

#### Outputs:

| Output | Type | Description |
| --- | --- | --- |
| normalized_bed | File | BED-format phenotype matrix containing inverse-normal transformed editing levels |
| site_metadata | File | CSV file containing summary statistics and filtering information for all sites |
| phenotype_pcs | File | Tab-delimited file containing phenotype principal components |

The normalized BED file is formatted as follows:

```text
chr     start     end     phenotype_id                    sample_1    sample_2    sample_3
chr1    1234566   1234567 chr1:1234566:1234567_ENSG...   -0.8432     0.1211      1.3448
chr1    1234889   1234890 chr1:1234889:1234890_ENSG...    0.5517    -1.2450      0.3028
```

where:

- `chr`: Chromosome name.
- `start`: Genomic start coordinate (0-based).
- `end`: Genomic end coordinate (1-based).
- `phenotype_id`: Unique site identifier.
- Remaining columns contain inverse-normal transformed editing phenotypes for each sample.

The phenotype PC file is formatted as follows:

```text
ID     PC1      PC2      PC3 ...
sample1
sample2
sample3
```

where each row corresponds to a sample and each column contains a phenotype principal component.

---

### Script: `scripts/transform.py`

Converts editing count matrices into normalized phenotype matrices for QTL mapping. The script filters low-variance editing sites, calculates editing proportions, imputes missing values, applies inverse normal transformation, and outputs a BED-format phenotype matrix. A metadata file containing filtering statistics for all sites is also generated.

**Inputs:**

| Argument | Type | Description |
|----------|------|------------|
| `--input`, `-i` | File | Editing matrix generated in Step 2 |
| `--output`, `-o` | String | Output BED filename |
| `--mad_threshold`, `-m` | Float | Minimum median absolute deviation required for a site to be retained (default: 0.003) |
| `--sample_list` | File | TSV file containing a `research_id` column specifying samples to retain |

**Outputs:**

- `<output>.gz` — BED-format matrix containing normalized editing phenotypes.
- `site_metadata.csv` — Metadata file containing summary statistics and filtering information for all sites.

The metadata file is formatted as follows:

```text
site,mean_editing,mad,n_missing,kept
chr1:1234566:1234567:ENSG00000123,0.27,0.011,14,1
chr1:1234889:1234890:ENSG00000456,0.73,0.001,3,0
```

where:

- `site`: Site identifier.
- `mean_editing`: Mean editing proportion across observed samples. Note that a pseudocount is applied here, so that editing levels are computed as: (edited_reads + 0.5) / (coverage + 0.5)
- `mad`: Median absolute deviation.
- `n_missing`: Number of samples with insufficient coverage (i.e. less than `mincov` from Step 2) at that editing site.
- `kept`: Indicator of whether the site passed filtering (`1`) or was removed (`0`).

---

### Script: `scripts/calculate_pcs.R`

Computes phenotype principal components from the normalized BED matrix generated by `transform.py`. The script performs principal component analysis (PCA) and automatically determines the number of informative PCs using the Gavish-Donoho method implemented in PCAtools.

**Inputs:**

| Argument | Type | Description |
|----------|------|------------|
| `--bed_file` | File | Normalized BED file generated by `transform.py` |
| `--output_prefix` | String | Prefix used for output filenames |

**Output:**

- `<output_prefix>_phenotype_PCs.tsv` — Tab-delimited file containing phenotype principal components for all samples.

**Output format:**

```text
ID     PC1      PC2      PC3 ...
sample1
sample2
sample3
```

where:

- `ID` is the sample identifier.
- `PC1`, `PC2`, etc. are phenotype principal components that can be used as covariates in downstream QTL analyses.

## Acknowledgments

This repository was built upon code developed in:

* Li, Q., Gloudemans, M.J., Geisinger, J.M. et al. RNA editing underlies genetic risk of common inflammatory diseases. Nature 608, 569–577 (2022). https://doi.org/10.1038/s41586-022-05052-x

* Original repository: https://github.com/vargasliqin/GTEx_edQTL
