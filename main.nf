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
    Reads extraction
*/
if(params.minimap2){
    process minimap2 {
      publishDir "/workspace/ieo4032/nanoid/data", mode: 'copy'
      input:
      output:
    
    script:
    """
        cp -r ${params.FASTQ} ${params.FASTQ}_DATA
        for i in ${params.FASTQ}_DNA/*fastq; do
            awk '{ if (NR%4 == 2) {gsub(/U/,"T",$1); print $1} else print }' $i > ${i%.fastq}.U2T.fastq;
            rm $i
        done
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


