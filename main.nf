#!/usr/bin/env nextflow
/*
    PROCESSES
*/
/*
    Reads extraction
*/
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
