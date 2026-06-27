#!/usr/bin/env python3

import os
import sys
import pandas as pd
from functools import reduce
import re
import json
import plotly.express as px


def main():
    composition_file_path = sys.argv[1]
        
    # Load the composition dataframe
    df  = pd.read_csv(composition_file_path, index_col=0)
    
    # Get the substraing for what level of variance this is
    var_level = composition_file_path.split("_")[2]
    
    # Melt df
    melted = df.reset_index().melt(id_vars=['index', ])
    melted = melted.sort_values("variable")

    # Make graph
    fig = px.line(melted, x = 'variable', y = 'value', color = 'index')

    fig.update_layout(title = f"Percent Composition Across Simulations with {var_level} Variance")

    # Update x-axis and y-axis labels
    fig.update_xaxes(title_text="Simulation", showticklabels=False)
    fig.update_yaxes(title_text='Percent Composition')
    fig.update_layout(legend_title_text='Cell Types (CTs)')

    
    fig.write_image("simulation_consistancy.png")
    
    
if __name__ == "__main__":
    main()
    
    
