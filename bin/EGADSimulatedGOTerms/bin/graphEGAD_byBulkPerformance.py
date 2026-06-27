#!/usr/bin/env python3

import pandas as pd
import sys
import plotly.graph_objects as go
import seaborn as sns
import matplotlib.pyplot as plt



def main():
    organism_part = sys.argv[1]
    melted_df_path = sys.argv[2]
    bulk_df_path = sys.argv[3]
    
    # load melted df
    melted_df = pd.read_csv(melted_df_path, index_col=0)
    # load bulk df
    bulk_df = pd.read_csv(bulk_df_path)
    # Load ctaffiliation_df
    ctaffiliation_df = pd.read_csv(ctaffiliation_path).fillna("NAN")
    
    # Add the CTAffiliation info to the melted df
    melted_df = pd.merge(left=melted_df, right=ctaffiliation_df, left_on="index", right_on="ID")
    
    # Get the organism part we are using
    op = melted_df_path.split("_")[0]
    # Make boxplot of melted df
    graph_AUCs_boxplot(melted_df, op)
    
    # Make avg lineplot of melted df
    graph_AUCs_lineplot(melted_df, op = op)
    
def graph_AUCs_boxplot(melted_df, op):
    """Graphs AUCS. 
    Make one panel per organism part. 
    Args:
        lo_processed_melted_dfs (list): _description_
    """
    # Defining a color palette
    custom_palette = {True: 'red', False: 'blue', 'NAN':'grey'}
    # Setting the custom color palette
    g = sns.FacetGrid(melted_df, col = 'type', col_order=["contrl", "exp"])
    g.map(sns.boxplot, 'variable', 'value', 'CTAffiliation', palette=custom_palette)
    # Set plot labels and title
    g.set_axis_labels('Compositional Variance Levels', 'AUC')
    # Create a custom legend
    g.add_legend()
    g.savefig(f"{op}_EGAD_avg_AUC_boxplot_CTAffiliation.png")
    
    plt.figure().clear()

def graph_AUCs_lineplot(melted_df, op):
    
    # get the average auc performance for each variance level between exp and contrl
    melted_mean_df = melted_df.groupby(['variable', 'type']).mean().reset_index()
    # Defining a color palette
    custom_palette = {True: 'red', False: 'blue', 'NAN':'grey'}

    
    # Creating a line plot using seaborn
    g = sns.FacetGrid(melted_df, col = 'type', col_order=["contrl", "exp"])
    g.map(sns.lineplot, 'variable', 'value', 'CTAffiliation', palette = custom_palette)
    # Set plot labels and title
    g.set_axis_labels('Compositional Variance Levels', 'AUC')
    # Create a custom legend
    g.add_legend()
    # Show the plot
    print(f"{op}_EGAD_avg_AUC_lineplot.png")
    plt.savefig(f"{op}_EGAD_avg_AUC_lineplot_CTAffiliation.png")
    
    
if __name__ == "__main__":
	main()