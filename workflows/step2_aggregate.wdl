version 1.0

task aggregate_samples {
    input {
        File EditingCountsFOFN 
        String output_file
        
        Int min_coverage = 20
        Int min_samples = 900
        
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

    echo "Validating sample consistency across chromosome matrices..."
    
    # Get reference header from first file
    REFERENCE_HEADER=$(zcat ~{chr_matrices[0]} | head -1)
    REFERENCE_FILE="~{chr_matrices[0]}"
    REFERENCE_COUNT=$(echo "$REFERENCE_HEADER" | wc -w)
    
    echo "  File 0 ($REFERENCE_FILE): $REFERENCE_COUNT columns"
    
    # Check all other files against reference
    MATRICES=(~{sep=' ' chr_matrices})
    ALL_MATCH=true
    
    for i in "${!MATRICES[@]}"; do
        if [ $i -gt 0 ]; then
            CURRENT_FILE="${MATRICES[$i]}"
            CURRENT_HEADER=$(zcat "$CURRENT_FILE" | head -1)
            CURRENT_COUNT=$(echo "$CURRENT_HEADER" | wc -w)
            
            if [ "$CURRENT_HEADER" != "$REFERENCE_HEADER" ]; then
                echo ""
                echo "ERROR: Sample mismatch detected!"
                echo "  Reference file: $REFERENCE_FILE"
                echo "  Current file:   $CURRENT_FILE"
                echo "  Reference: $REFERENCE_COUNT columns"
                echo "  Current:   $CURRENT_COUNT columns"
                echo ""
                echo "This will cause column misalignment in the merged matrix!"
                ALL_MATCH=false
                exit 1
            else
                echo "  File $i ($CURRENT_FILE): $CURRENT_COUNT columns ✓"
            fi
        fi
    done
    
    if [ "$ALL_MATCH" = true ]; then
        echo ""
        echo "✓ All ${#MATRICES[@]} files have identical headers ($REFERENCE_COUNT columns)"
    fi

    # Validation passed - now merge
    echo ""
    echo "Merging chromosome matrices..."
    
    # First file: keep header
    zcat ~{chr_matrices[0]} > merged.txt
    echo "  Added ~{chr_matrices[0]} (with header)"

    # Remaining files: skip header
    for i in "${!MATRICES[@]}"; do
        if [ $i -gt 0 ]; then
            echo "  Adding ${MATRICES[$i]} (data only)..."
            zcat "${MATRICES[$i]}" | tail -n +2 >> merged.txt
        fi
    done

    echo ""
    echo "Compressing output..."
    gzip merged.txt
    mv merged.txt.gz ~{output_file}.gz
    
    echo "Merge complete"
    echo "Output file: ~{output_file}.gz"
    ls -lh ~{output_file}.gz
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
        Array[String] chr_names = [ "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22"]

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
