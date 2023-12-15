#!/usr/bin/env python3
import os
import sys
import pandas as pd


def main():
    bootstrap = sys.argv[1]
    organism_part_path = sys.argv[2]
    organism_part = sys.argv[3]
    sample_size = int(sys.argv[4])
    
    # Load adata
    df = pd.read_csv(organism_part_path)
    
    # Generate a sample with replacement
    sub_sampled_data = df.sample(n=sample_size, replace=True)
    # save subsample
    sub_sampled_data.to_csv(f"{organism_part}_{bootstrap}.csv.gz", index=False)
    


if __name__ == '__main__':
    main()
    
    
