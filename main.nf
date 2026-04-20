#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { HCVPIPE } from './workflows/hcvpipe'

workflow {
    HCVPIPE()
}
