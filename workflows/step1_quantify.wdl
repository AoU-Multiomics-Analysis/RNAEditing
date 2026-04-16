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
        ln -s ~{bam_file} input.bam
        ln -s ~{bam_index} input.bam.bai
        
        ln -s ~{reference_genome} reference.fasta
        ln -s ~{reference_genome_index} reference.fasta.fai
        
        env PERLLIB=/opt/scripts \
        perl /opt/scripts/query_editing_levels.pl \
            --bam input.bam \
            --output ~{sample_id}.rnaediting_op \
            --genome reference.fasta \
            --sites /opt/scripts/All.AG.stranded.annovar.Hg38_multianno.AnnoAlu.AnnoRep.NR.bed.gz
    >>>
    
    runtime {
        docker: "ghcr.io/aou-multiomics-analysis/rnaediting:main"
        memory: "~{memory}GB"
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_threads}"
    }
    
    output {
        File editing_counts = "~{sample_id}.rnaediting_op"
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
    
    output {
        File editing_counts = quantify_editing_single_sample.editing_counts
    }
}