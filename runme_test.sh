#!/bin/bash

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: bash runme_test.sh <module_name> [resume]" >&2
    exit 1
fi

source /etc/profile
module load Java/23.0.2 nextflow/25.10.0 apptainer

cd /fs1/jonas/src/virpipa
git pull

if [[ "${2:-}" == "resume" ]]; then
    nextflow run test_module.nf -profile hpc,apptainer --module "$1" -resume
else
    nextflow run test_module.nf -profile hpc,apptainer --module "$1"
fi
