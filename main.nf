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
            cd nanoid
            Rscript read.extraction.R
        """
}
