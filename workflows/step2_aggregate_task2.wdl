version 1.0

task merge_chromosome_matrices {
  input {
    Array[File] chr_matrices
    String output_file

    Int memory
    Int disk_space
    Int num_threads
  }

  command <<<
    set -euo pipefail

    # First file: keep header
    zcat ~{chr_matrices[0]} > merged.txt

    for f in ~{sep=' ' chr_matrices}; do
      if [[ "$f" != "~{chr_matrices[0]}" ]]; then
        zcat "$f" | tail -n +2 >> merged.txt
      fi
    done

    gzip merged.txt
    mv merged.txt.gz ~{output_file}.gz
  >>>

  runtime {
    docker: "ubuntu:22.04"
    memory: "~{memory}GB"
    disks: "local-disk ~{disk_space} HDD"
    cpu: "~{num_threads}"
  }

  output {
    File merged_matrix = "~{output_file}.gz"
  }
}

workflow merge_all_chr_workflow {
  input {
    Array[File] chr_matrices
    String output_file

    Int memory
    Int disk_space
    Int num_threads
  }

  call merge_chromosome_matrices {
    input:
      chr_matrices = chr_matrices,
      output_file = output_file,
      memory = memory,
      disk_space = disk_space,
      num_threads = num_threads
  }

  output {
    File merged_matrix = merge_chromosome_matrices.merged_matrix
  }
}