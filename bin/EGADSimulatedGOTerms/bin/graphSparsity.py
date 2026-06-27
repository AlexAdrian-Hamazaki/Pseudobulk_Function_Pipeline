#!/usr/bin/env python3

import os
import sys
import pandas as pd
from functools import reduce
import plotly.express as px


def main():
    #
    ct_profile_path = sys.argv[1]
    organism_part_name = sys.argv[2]
    
    # load ct profile
    df = pd.read_csv(ct_profile_path,index_col=0)
    
    # Count the number of zeros in each row
    zeros_in_row = (df == 0).sum(axis=1)
    # Calc percent of zeros    
    percent_zeros = zeros_in_row/df.shape[1]
    
    # Write a bar chart
    fig = px.bar(percent_zeros, x=percent_zeros.index, y=percent_zeros.values, labels={'x': 'Cell Type Profile', 'y': 'Percent Zeros'}, title=f'Sparsity of CT Profiles from {organism_part_name}')
    fig.write_image(f"{organism_part_name}_sparsity_graph.png")

if __name__ == "__main__":
    main()