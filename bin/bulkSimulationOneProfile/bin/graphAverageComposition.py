#!/usr/bin/env python3

import os
import sys
import pandas as pd
from functools import reduce
import re
import json
import plotly.express as px


def main():
    lo_composition_file_paths = sys.argv[1:]
    
    # Load all the pandas dataframes
    lodfs = [pd.read_csv(file_path, index_col=0) for file_path in lo_composition_file_paths]
    
    # Get the num of each simulation
    lo_sim_ids = [path.split("_")[2] for path in lo_composition_file_paths]
    
    
    # Concat all the dfs and change column names
    merged_df = pd.concat(lodfs, axis = 1)
    merged_df.columns = lo_sim_ids
    
    
    
    # Melt df
    melted = merged_df.reset_index().melt(id_vars=['index'])
    melted = melted.sort_values("variable", ascending=True)
    
    # Extract the organism part substring
    op_substring = lo_composition_file_paths[0].split("_")[3]
    
    # Plot
    fig = px.line(melted, x = 'variable', y = 'value', color = 'index')
    fig.update_layout(title = f"Average Percent Composition of {op_substring} Across Variance Levels")

    # Update x-axis and y-axis labels
    fig.update_xaxes(title_text='Variance Level')
    fig.update_yaxes(title_text='Average Percent Composition')
    fig.update_layout(legend_title_text='Cell Types')


    
    # Save the figure as a PNG file
    fig.write_image("avg_composition_graph.png", format="png")
    
    
if __name__ == "__main__":
    main()