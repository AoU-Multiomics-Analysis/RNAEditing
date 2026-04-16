version 1.0

task aggregate_samples {
    input {
        #Array[File] individual_count_files
        File EditingCountsFOFN 
        String output_file
        
        Int min_coverage = 20
        Int min_samples = 450
        
        Int memory
        Int disk_space
        Int num_threads
    }
    
    command <<<
    mkdir -p input_files
    mkdir -p localized

    export GSUTIL_PARALLEL_PROCESS_COUNT=32
    export GSUTIL_PARALLEL_THREAD_COUNT=8

    awk '{print $1}' ~{EditingCountsFOFN} | grep -v '^$' > file_paths.txt

    echo "=== file_paths.txt count ==="
    wc -l file_paths.txt

    xargs -a file_paths.txt -n 100 sh -c 'gsutil -m cp "$@" localized/' sh

    echo "=== localized file count ==="
    find localized -maxdepth 1 -type f | wc -l

    find "$(pwd)/localized" -maxdepth 1 -type f | sort > filelist.txt

    perl /opt/scripts/combine_sample_matrices.pl \
        --input localized \
        --output ~{output_file} \
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
        File aggregated_matrix = "~{output_file}"
    }
}

workflow aggregate_workflow {
    input {
        #Array[File] count_files
        File EditingCountsFOFN 
 
        String output_file
        
        Int min_coverage = 20
        Int min_samples = 450
        
        Int memory
        Int disk_space
        Int num_threads
    }
    
    call aggregate_samples {
        input:
            EditingCountsFOFN = EditingCountsFOFN,
            output_file = output_file,
            min_coverage = min_coverage,
            min_samples = min_samples,
            memory = memory,
            disk_space = disk_space,
            num_threads = num_threads
    }
    
    output {
        File combined_matrix = aggregate_samples.aggregated_matrix
    }
}
