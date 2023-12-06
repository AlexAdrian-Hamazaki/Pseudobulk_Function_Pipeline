#!/usr/bin/env python3

import pandas as pd
import os
import functools
import sys
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt



def main():
    lo_melted_df_path = sys.argv[1:]

    for melted_df_path in lo_melted_df_path:
        # Get the organism part we are using
        op = melted_df_path.split("_")[0]    
        
        # load melted df
        melted_df = pd.read_csv(melted_df_path, index_col=0)
                
        # Make boxplot of melted df
        graph_AUCs_boxplot(melted_df, op)
        
        # Make avg lineplot of melted df
        graph_AUCs_lineplot(melted_df=melted_df, op = op)
    
def graph_AUCs_boxplot(melted_df, op):
    """Graphs AUCS. 
    Make one panel per organism part. 
    Args:
        lo_processed_melted_dfs (list): _description_
    """
    # Defining a color palette
    custom_palette = {'exp': 'red', 'contrl': 'blue'}
    # Setting the custom color palette
    sns.set_palette([custom_palette[val] for val in melted_df['type'].unique()])# Set plot labels and title
    
    # Creating a boxplot using seaborn
    sns.boxplot(data=melted_df, x='variable', y='value', hue='type')
    # Set plot labels and title
    plt.xlabel('Variance')
    plt.ylabel('AUC')
    plt.title('AUC across Variance Levels')
    plt.savefig(f"{op}_EGAD_avg_AUC_boxplot.png")
    plt.figure().clear()
    

def graph_AUCs_lineplot(melted_df, op):
    plt.figure().clear()

    print(melted_df.head())
    # get the average auc performance for each variance level between exp and contrl
    melted_mean_df = melted_df.groupby(['variable', 'type']).mean().reset_index()
    
    # Defining a color palette
    custom_palette = {'exp': 'red', 'contrl': 'blue'}
    # Setting the custom color palette
    sns.set_palette([custom_palette[val] for val in melted_mean_df['type'].unique()])# Set plot labels and title
    # Creating a line plot using seaborn
    sns.lineplot(data=melted_mean_df, x='variable', y='value', hue='type')
    plt.xlabel('Variable')
    plt.ylabel('Value')
    plt.title('Line Plot Grouped by Type')
    # Show the plot
    print(f"{op}_EGAD_avg_AUC_lineplot.png")
    plt.savefig(f"{op}_EGAD_avg_AUC_lineplot.png")
    plt.figure().clear()

    
    
if __name__ == "__main__":
	main()