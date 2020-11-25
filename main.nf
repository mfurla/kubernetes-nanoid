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
    if(${paramsparams.minimap2}='true')
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
}

/*
    Reads extraction
*/

    process readExtraction{
      publishDir ".", mode: 'copy'
    
      input:

      output:
      val "readExtraction" into readExtractionFlag
    
    script:
    if(${paramsparams.readExtraction}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript read.extraction.R
    """
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Alignment extraction
*/

    process alignmentExtraction {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "alignmentExtraction" into alignmentExtractionFlag
    
    script:
    if(${paramsparams.alignmentExtraction}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript alignment.extraction.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Sequencing summary
*/

    process sequencingSummary {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "sequencingSummary" into sequencingSummaryFlag
    
    script:
    if(${paramsparams.sequencingSummary}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript sequencing.summary.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Alignment Reconstruction
*/

    process alignmentReconstruction {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "alignmentReconstruction" into alignmentReconstructionFlag
    
    script:
    if(${paramsparams.alignmentReconstruction}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript alignment.reconstruction.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    5-mer Alignment Reconstruction
*/

    process fivemerAlignmentReconstruction {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "fivemerAlignmentReconstruction" into fivemerAlignmentReconstructionFlag
    
    script:
    if(${paramsparams.fivemerAlignmentReconstruction}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript five.mer.alignment.reconstruction.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Trace Model
*/

    process traceModel {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "traceModel" into traceModelFlag
    
    script:
    if(${paramsparams.traceModel}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript trace.model.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Trace Model Add On
*/

    process traceModelAddOn {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "traceModelAddOn" into traceModelAddOnFlag
    
    script:
    if(${paramsparams.traceModelAddOn}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript trace.model.add.on.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Raw Signal
*/

    process rawSignal {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "rawSignal" into rawSignalFlag
    
    script:
    if(${paramsparams.rawSignal}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Raw Signal 5-mers
*/

    process rawSignalFiveMers {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "rawSignalFiveMers" into rawSignalFiveMersFlag
    
    script:
    if(${paramsparams.rawSignalFiveMers}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.five.mers.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}

/*
    Raw Signal 5-mers Add On
*/

    process rawSignalFiveMersAddOn {
      publishDir ".", mode: 'copy'

      input:

      output:
      val "rawSignalFiveMersAddOn" into rawSignalFiveMersAddOnFlag
    
    script:
    if(${paramsparams.rawSignalFiveMersAddOn}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript raw.signal.five.mers.add.on.R
    """  
    else
    """
        echo "Skipped"
    """
    }
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
    if(${paramsparams.dataFormatting}='true')
    """
        cd /workspace/ieo4032/nanoid
        Rscript data.formatting.R
    """  
    else
    """
        echo "Skipped"
    """
    }
}