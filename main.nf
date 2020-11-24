#!/usr/bin/env nextflow

/*
    CHANNELS
*/

/*
        Fast5 directory
*/
if( params.fast5 ){
    Channel
        .fromPath(params.fast5, checkIfExists:true)
        .into{FAST5pathGuppy}
}
else {
    exit 1, "No fast5 files."
}

/*
    PROCESSES
*/
/*
        Guppy basecalling
*/

process guppy {
    
    publishDir ".", mode: 'copy'

    input:
        file 'fast5path' from FAST5pathGuppy

    output:
        set val("${sample}"), file("guppy") into guppy_outputs_pycoqc, guppy_outputs_minimap, guppy_outputs_nanopolish

    script:
        def gpu_opts = ""
        if (params.GPU == "true") {
            gpu_opts = "-x 'cuda:0' --gpu_runners_per_device ${params.guppy_runners_per_device} --chunks_per_runner ${params.guppy_chunks_per_runner} --chunk_size ${params.guppy_chunk_size}"
        }
        """
            guppy_basecaller -i fast5path -s guppy --fast5_out true ${gpu_opts} --recursive --num_callers ${task.cpus} --disable_pings --reverse_sequence true --u_substitution true --trim_strategy rna --flowcell ${params.flowcell} --kit ${params.kit}
        """
}
