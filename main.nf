#!/usr/bin/env nextflow
/*
    PROCESSES
*/

/*
    Minimap2 - rewrite better with proper Naxtflow functionalities
*/

    process minimap2 {
      publishDir ".", mode: 'copy'
    
      input:

      output:
      val "minimap2" into minimap2Flag
    
    script:
    if(params.minimap2='true')
    """
        cd /workspace/ieo4032/nanoid
        minimap2 -ax splice -k14 --secondary=no -I1G -t${task.cpus} Mus_musculus.GRCm38.dna.primary_assembly.fa ${params.FASTQ}_DNA/*.fastq > data/minimap.sam
        samtools view data/minimap.sam -bhq20 -t Mus_musculus.GRCm38.dna.primary_assembly.fa.fai -F 2324 | samtools sort -o data/minimap.bam
        samtools index data/minimap.bam data/minimap.bam.bai        
    """
    else
    """
        echo "Skipped"
    """
    }

/*
    Reads extraction
*/

    process readExtraction{
      publishDir ".", mode: 'copy'
    
      input:
      val minimap2Flag from minimap2Flag

      output:
      val "readExtraction" into readExtractionFlag
    
    script:
    if(params.readExtraction='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript read.extraction.R
    """
    else
    """
        echo "Skipped"
    """
    }

/*
    Alignment extraction
*/

    process alignmentExtraction {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag

      output:
      val "alignmentExtraction" into alignmentExtractionFlag
    
    script:
    if(params.alignmentExtraction='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript alignment.extraction.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Sequencing summary
*/

    process sequencingSummary {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag

      output:
      val "sequencingSummary" into sequencingSummaryFlag
    
    script:
    if(params.sequencingSummary='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript sequencing.summary.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Alignment Reconstruction
*/

    process alignmentReconstruction {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag

      output:
      val "alignmentReconstruction" into alignmentReconstructionFlag
    
    script:
    if(params.alignmentReconstruction='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript alignment.reconstruction.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    5-mer Alignment Reconstruction
*/

    process fivemerAlignmentReconstruction {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag
      val alignmentReconstructionFlag from alignmentReconstructionFlag

      output:
      val "fivemerAlignmentReconstruction" into fivemerAlignmentReconstructionFlag
    
    script:
    if(params.fivemerAlignmentReconstruction='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript five.mer.alignment.reconstruction.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Trace Model
*/

    process traceModel {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag
      val alignmentReconstructionFlag from alignmentReconstructionFlag
      val fivemerAlignmentReconstructionFlag from fivemerAlignmentReconstructionFlag

      output:
      val "traceModel" into traceModelFlag
    
    script:
    if(params.traceModel='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript trace.model.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Trace Model Add On
*/

    process traceModelAddOn {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag
      val alignmentReconstructionFlag from alignmentReconstructionFlag
      val fivemerAlignmentReconstructionFlag from fivemerAlignmentReconstructionFlag
      val traceModelFlag from traceModelFlag

      output:
      val "traceModelAddOn" into traceModelAddOnFlag
    
    script:
    if(params.traceModelAddOn='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript trace.model.add.on.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Raw Signal
*/

    process rawSignal {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag
      val alignmentReconstructionFlag from alignmentReconstructionFlag
      val fivemerAlignmentReconstructionFlag from fivemerAlignmentReconstructionFlag
      val traceModelFlag from traceModelFlag
      val traceModelAddOnFlag from traceModelAddOnFlag

      output:
      val "rawSignal" into rawSignalFlag
    
    script:
    if(params.rawSignal='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Raw Signal 5-mers
*/

    process rawSignalFiveMers {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag
      val alignmentReconstructionFlag from alignmentReconstructionFlag
      val fivemerAlignmentReconstructionFlag from fivemerAlignmentReconstructionFlag
      val traceModelFlag from traceModelFlag
      val traceModelAddOnFlag from traceModelAddOnFlag
      val rawSignalFlag from rawSignalFlag

      output:
      val "rawSignalFiveMers" into rawSignalFiveMersFlag
    
    script:
    if(params.rawSignalFiveMers='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.five.mers.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Raw Signal 5-mers Add On
*/

    process rawSignalFiveMersAddOn {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag
      val alignmentReconstructionFlag from alignmentReconstructionFlag
      val fivemerAlignmentReconstructionFlag from fivemerAlignmentReconstructionFlag
      val traceModelFlag from traceModelFlag
      val traceModelAddOnFlag from traceModelAddOnFlag
      val rawSignalFlag from rawSignalFlag
      val rawSignalFiveMersFlag from rawSignalFiveMersFlag

      output:
      val "rawSignalFiveMersAddOn" into rawSignalFiveMersAddOnFlag
    
    script:
    if(params.rawSignalFiveMersAddOn='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.five.mers.add.on.R
    """  
    else
    """
        echo "Skipped"
    """
    }

/*
    Raw Signal 5-mers Add On
*/

    process dataFormatting {
      publishDir ".", mode: 'copy'

      input:
      val minimap2Flag from minimap2Flag
      val readExtractionFlag from readExtractionFlag
      val alignmentExtractionFlag from alignmentExtractionFlag
      val sequencingSummaryFlag from sequencingSummaryFlag
      val alignmentReconstructionFlag from alignmentReconstructionFlag
      val fivemerAlignmentReconstructionFlag from fivemerAlignmentReconstructionFlag
      val traceModelFlag from traceModelFlag
      val traceModelAddOnFlag from traceModelAddOnFlag
      val rawSignalFlag from rawSignalFlag
      val rawSignalFiveMersFlag from rawSignalFiveMersFlag
      val rawSignalFiveMersAddOnFlag from rawSignalFiveMersAddOnFlag
        
      output:
      val "dataFormatting" into dataFormattingFlag
    
    script:
    if(params.dataFormatting='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript data.formatting.R
    """  
    else
    """
        echo "Skipped"
    """
    }