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
        minimap2 -ax splice -k14 --secondary=no -I1G -t${task.cpus} Mus_musculus.GRCm38.dna.primary_assembly.fa ${params.FASTQ}_DNA/*.fastq > data/minimap.sam
        samtools view data/minimap.sam -bhq20 -t Mus_musculus.GRCm38.dna.primary_assembly.fa.fai -F 2324 | samtools sort -o data/minimap.bam
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

/*
    Sequencing summary
*/
if(params.sequencingSummary){
    process sequencingSummary {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript sequencing.summary.R
    """  
    }
}

/*
    Alignment Reconstruction
*/
if(params.alignmentReconstruction){
    process alignmentReconstruction {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript alignment.reconstruction.R
    """  
    }
}

/*
    5-mer Alignment Reconstruction
*/
if(params.fivemerAlignmentReconstruction){
    process fivemerAlignmentReconstruction {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript five.mer.alignment.reconstruction.R
    """  
    }
}

/*
    Trace Model
*/
if(params.traceModel){
    process traceModel {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript trace.model.R
    """  
    }
}

/*
    Trace Model Add On
*/
if(params.traceModelAddOn){
    process traceModelAddOn {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript trace.model.add.on.R
    """  
    }
}

/*
    Raw Signal
*/
if(params.rawSignal){
    process rawSignal {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.R
    """  
    }
}

/*
    Raw Signal 5-mers
*/
if(params.rawSignalFiveMers){
    process rawSignalFiveMers {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.five.mers.R
    """  
    }
}

/*
    Raw Signal 5-mers Add On
*/
if(params.rawSignalFiveMersAddOn){
    process rawSignalFiveMersAddOn {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.five.mers.add.on.R
    """  
    }
}

/*
    Raw Signal 5-mers Add On
*/
if(params.dataFormatting){
    process dataFormatting {
      input:
      output:
    
    script:
    """
        cd /workspace/ieo4032/nanoid
        Rscript data.formatting.R
    """  
    }
}