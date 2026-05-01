version 1.0

task aggregate_samples {
    input {
        File EditingCountsFOFN 
        String output_file
        
        Int min_coverage = 20
        Int min_samples = 450
        
        Int memory
        Int disk_space
        Int num_threads
    }
    
    command <<<
    set -euo pipefail

    mkdir -p localized

    export GSUTIL_PARALLEL_PROCESS_COUNT=32
    export GSUTIL_PARALLEL_THREAD_COUNT=8

    awk '{print $1}' ~{EditingCountsFOFN} | grep -v '^$' > file_paths.txt

    echo "=== file_paths.txt count ==="
    wc -l file_paths.txt

    xargs -a file_paths.txt -n 100 sh -c 'gsutil -m cp "$@" localized/' sh

    echo "=== localized file count ==="
    find localized -maxdepth 1 -type f | wc -l

    # Aggregate
    perl /opt/scripts/combine_sample_matrices.pl \
        --input localized \
        --output ~{output_file}.gz \
        --mincov ~{min_coverage} \
        --minsamps ~{min_samples}
    >>>
        
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/rnaediting/aggregation:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }
    
    output {
        File aggregated_matrix = "~{output_file}.gz"
    }
}

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

workflow aggregate_workflow {
    input {
        # One FOFN per chromosome (must match WDL1 outputs)
        Array[File] chr_fofns
        Array[String] chr_names   # ["chr1", ..., "chr22"]

        String output_prefix
        
        Int min_coverage = 20
        Int min_samples = 450
        
        Int memory
        Int disk_space
        Int num_threads
    }

    scatter (i in range(length(chr_fofns))) {

        call aggregate_samples as aggregate_per_chr {
            input:
                EditingCountsFOFN = chr_fofns[i],
                output_file = output_prefix + "." + chr_names[i],
                min_coverage = min_coverage,
                min_samples = min_samples,
                memory = memory,
                disk_space = disk_space,
                num_threads = num_threads
        }
    }

    call merge_chromosome_matrices {
        input:
            chr_matrices = aggregate_per_chr.aggregated_matrix,
            output_file = output_prefix + ".all_chr",
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    output {
        Array[File] per_chr_matrices = aggregate_per_chr.aggregated_matrix
        File combined_matrix = merge_chromosome_matrices.merged_matrix
    }
}