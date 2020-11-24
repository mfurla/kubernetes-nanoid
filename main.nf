#!/usr/bin/env nextflow
/*
    PROCESSES
*/

/*
    Reads extraction
*/
if(params.readExtraction){
    process readExtraction{
        
        publishDir ".", mode: 'copy'
    
        input:
        output:
    
        script:
            """
                cd /workspace/ieo4032/nanoid
                Rscript read.extraction.R
            """
    }
}

/*
    Reads extraction - rewrite better with proper Naxtflow functionalities
*/
if(params.minimap2){
    process minimap2 {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        minimap2 -ax splice -k14 --secondary=no -I1G -t${task.cpus} GRCm38.primary_assembly.genome.fa ${params.FASTQ}_DNA/*.fastq > data/minimap.sam
        samtools view data/minimap.sam -bhq20 -t GRCm38.primary_assembly.genome.fa.fai -F 2324 | samtools sort -o data/minimap.bam
        samtools index data/minimap.bam data/minimap.bam.bai        
    """  
    }
}

/*
    Alignment extraction
*/
if(params.alignmentExtraction){
    process alignmentExtraction {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript alignment.extraction.R
    """  
    }
}


