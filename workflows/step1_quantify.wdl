version 1.0

task quantify_editing_single_sample {
  input {
    File bam_file
    File bam_index
    File reference_genome
    File reference_genome_index
    String sample_id

    Int memory
    Int disk_space
    Int num_threads
  }

  command <<<
    set -euo pipefail

    ln -s ~{bam_file} input.bam
    ln -s ~{bam_index} input.bam.bai

    ln -s ~{reference_genome} reference.fasta
    ln -s ~{reference_genome_index} reference.fasta.fai

    PERLLIB=/opt/scripts perl /opt/scripts/query_editing_levels.pl \
      --bam input.bam \
      --output ~{sample_id}.rnaediting_op.gz \
      --genome reference.fasta \
      --sites /opt/scripts/augmented_editing_sites_REDIportal_plus_original_annotated_autosomes.bed.gz
  >>>

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/rnaediting/quantification:main"
    memory: "~{memory}GB"
    disks: "local-disk ~{disk_space} HDD"
    cpu: "~{num_threads}"
  }

  output {
    File editing_counts_all = "~{sample_id}.rnaediting_op.gz"
  }
}

task split_editing_by_chr {
  input {
    File editing_counts_all
    String sample_id

    Int memory
    Int disk_space
    Int num_threads
  }

  command <<<
    set -euo pipefail
    mkdir -p per_chr

    python - <<'PY'
import gzip, os

inp = "~{editing_counts_all}"
outdir = "per_chr"
os.makedirs(outdir, exist_ok=True)

outs = {}
header = None

def out_handle(name):
    if name not in outs:
        fn = os.path.join(outdir, f"~{sample_id}.{name}.rnaediting_op.gz")
        h = gzip.open(fn, "wt")
        h.write(header)
        outs[name] = h
    return outs[name]

with gzip.open(inp, "rt") as f:
    header = f.readline()

    # pre-create outputs: chr1..chr22
    for i in range(1, 23):
        out_handle(f"chr{i}")

    for line in f:
        if not line.strip():
            continue
        chrname = line.split("\t", 1)[0]

        # since sites list is autosomes only, keep only chr1..chr22
        if chrname.startswith("chr"):
            rest = chrname[3:]
            if rest.isdigit():
                n = int(rest)
                if 1 <= n <= 22:
                    out_handle(f"chr{n}").write(line)

for h in outs.values():
    h.close()
PY
  >>>

  runtime {
    docker: "ghcr.io/aou-multiomics-analysis/rnaediting/quantification:main"
    memory: "~{memory}GB"
    disks: "local-disk ~{disk_space} HDD"
    cpu: "~{num_threads}"
  }

  output {
    Array[File] editing_counts_chr = glob("per_chr/~{sample_id}.chr*.rnaediting_op.gz")
  }
}

workflow quantify_editing_single_sample_workflow {
  input {
    File bam_file
    File bam_index
    String sample_id

    File reference_genome
    File reference_genome_index

    Int memory
    Int disk_space
    Int num_threads
  }

  call quantify_editing_single_sample {
    input:
      bam_file = bam_file,
      bam_index = bam_index,
      sample_id = sample_id,
      reference_genome = reference_genome,
      reference_genome_index = reference_genome_index,
      memory = memory,
      disk_space = disk_space,
      num_threads = num_threads
  }

  call split_editing_by_chr {
    input:
      editing_counts_all = quantify_editing_single_sample.editing_counts_all,
      sample_id = sample_id,
      memory = memory,
      disk_space = disk_space,
      num_threads = num_threads
  }

  output {
    File editing_counts_all = quantify_editing_single_sample.editing_counts_all
    Array[File] editing_counts_chr = split_editing_by_chr.editing_counts_chr
  }
}