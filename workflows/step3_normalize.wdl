version 1.0

task normalize_transform {
    input {
        File input_matrix
        File sample_list
        String output_file
        Int window_size = 100000

        Int memory
        Int disk_space
        Int num_threads
    }

    command <<<
        python /opt/scripts/transform.py \
            --input ~{input_matrix} \
            --output ~{output_file} \
            --window ~{window_size} \
            --sample_list ~{sample_list}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/rnaediting/quantification:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }

    output {
        File normalized_bed = "~{output_file}"
        File site_metadata = "site_metadata.csv"
    }
}

task calculate_pcs {
    input {
        File bed_file
        String output_prefix

        Int memory
        Int disk_space
        Int num_threads
    }

    command <<<
        Rscript /opt/scripts/calculate_pcs.R \
            --bed_file ~{bed_file} \
            --output_prefix ~{output_prefix}
    >>>

    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/rnaediting/quantification:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }

    output {
        File phenotype_pcs = "~{output_prefix}_phenotype_PCs.tsv"
    }
}

workflow normalize_workflow {
    input {
        File matrix_file
        File sample_list
        String output_prefix
        Int window_size = 100000

        Int memory
        Int disk_space
        Int num_threads
    }

    call normalize_transform {
        input:
            input_matrix = matrix_file,
            sample_list = sample_list,
            output_file = "~{output_prefix}.bed.gz",
            window_size = window_size,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    call calculate_pcs {
        input:
            bed_file = normalize_transform.normalized_bed,
            output_prefix = output_prefix,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }

    output {
        File normalized_bed = normalize_transform.normalized_bed
        File site_metadata = normalize_transform.site_metadata
        File phenotype_pcs = calculate_pcs.phenotype_pcs
    }
}
