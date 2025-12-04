#!/usr/bin/env nextflow
nextflow.enable.dsl=2
//nextflow.preview.output = true

// Modules
include { make_input } from "./modules/make_input.nf"
include { run_stai } from "./modules/run_stai.nf"
include { make_output } from "./modules/make_output.nf"


// Workflow
workflow {
    spdata = Channel.fromPath("${params.adata_sp}")
    scdata = Channel.fromPath("${params.adata_sc}")

    input = make_input(spdata, scdata)
    imputed = run_stai(input.batches.flatten())
    results = make_output(imputed.collect())

    //publish:
    //results >> 'out'

    emit:
    results_ch = results

}

// Output
//output {
//    out { path "./imputed" }
//}

