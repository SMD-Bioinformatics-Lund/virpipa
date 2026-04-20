#!/bin/bash
# Run the virpipa Nextflow pipeline on Hopper HPC
# Usage: bash runme.sh [resume]

set -e

# Load required modules
source /etc/profile
module load Java/23.0.2 nextflow/25.10.0 apptainer

# Change to pipeline directory
cd /fs1/jonas/src/virpipa

# Run nextflow
if [ "$1" = "resume" ]; then
    nextflow run main.nf -profile hpc,apptainer --input assets/test_samplesheet.csv --ref_dir refgenomes -resume
else
    nextflow run main.nf -profile hpc,apptainer --input assets/test_samplesheet.csv --ref_dir refgenomes
fi
