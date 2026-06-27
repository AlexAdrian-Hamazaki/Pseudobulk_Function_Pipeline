#!/bin/bash

# Check arguments
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <num_bootstraps> <publishDir> [max_parallel]"
    exit 1
fi

num_bootstraps="$1"
publishDir="$2"
max_parallel="${3:-4}"   # default parallel jobs = 4

run_bootstrap () {
    iteration=$1
    publishDir=$2

    echo "~~~~~~~~~~~~Running $iteration~~~~~~~~~~~~~"

    nextflow run bulkPipe.nf \
        --publish "$publishDir/$iteration" \
        -entry run_on_all_Gtex
}

export -f run_bootstrap

seq 1 "$num_bootstraps" | \
    xargs -n1 -P "$max_parallel" -I{} bash -c 'run_bootstrap "$@"' _ {} "$publishDir"

# Run using this command
# ./bootstrap_pipe.sh 100 data/bootstrapped maxforks