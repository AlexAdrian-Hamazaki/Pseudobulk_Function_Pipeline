#!/usr/bin/env python3


import sys
import scanpy as sc
import sys
import pandas as pd

def main():
    
    sample_size_ch = sys.argv[1]
    variance_ch = sys.argv[2]
    num_simulations_ch = sys.argv[3]
    lo_csv_paths = sys.argv[4:]
    
    print(lo_csv_paths)
    
    # load each csv path in a list
    loCSVs = [pd.read_csv(csv_path, index_col=0) for csv_path in lo_csv_paths]

    
    # Concatenate the dataframes
    simulated_bulk_dataset = pd.concat(loCSVs, axis = 1)
    
    print(simulated_bulk_dataset.shape)
    print(f"simulated_ss${sample_size_ch}_var${variance_ch}_nsim${num_simulations_ch}.csv")
    
    simulated_bulk_dataset.to_csv(f"simulated_ss{sample_size_ch}_var{variance_ch}_nsim{num_simulations_ch}.csv")
    


if __name__ == "__main__":
    main()