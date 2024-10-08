#!/bin/bash

# Check if num_bootstraps and publishDir are provided
if [ -z "$1" ] || [ -z "$2" ]; then
    echo "Usage: $0 <num_bootstraps> <publishDir>"
    exit 1
fi

num_bootstraps="$1"
publishDir="$2"

# Iterate through a while loop
iteration=1

while [ "$iteration" -le "$num_bootstraps" ]; do
	echo "\n~~~~~~~~~~~~Running $iteration~~~~~~~~~~~~~"
    # Modify the "nextflow.config" file
    sed -i "s|^params.publish.*|params.publish = \"$publishDir/$iteration\"|" nextflow.config

    # Run the "nextflow" command
    nextflow main.nf

    # Increment the iteration counter
    ((iteration++))
done


# Run using this command
# ./bootstrap_pipe.sh 100 data/bootstrapped