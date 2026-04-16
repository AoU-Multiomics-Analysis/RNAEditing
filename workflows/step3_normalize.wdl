version 1.0

task normalize_transform {
    input {
        File input_matrix
        
        String output_file
        
        Int memory
        Int disk_space
        Int num_threads
    }
    
    command <<<
        python /opt/scripts/transform.py \
            --input ~{input_matrix} \
            --output ~{output_file}
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/rnaediting:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }
    
    output {
        File normalized_bed = "~{output_file}"
    }
}

workflow normalize_workflow {
    input {
        File matrix_file
        
        String output_file
        
        Int memory
        Int disk_space
        Int num_threads
    }
    
    call normalize_transform {
        input:
            input_matrix = matrix_file,
            output_file = output_file,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }
    
    output {
        File normalized_bed = normalize_transform.normalized_bed
    }
}