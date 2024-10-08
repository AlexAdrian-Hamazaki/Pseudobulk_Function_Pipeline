#!/usr/bin/env python3

import sys
import pandas as pd

def main():
	path_to_all_ct_pros = sys.argv[1]
	
	ct_pros = pd.read_csv(path_to_all_ct_pros, index_col=0)
 
	# Fix index names
	ct_pros.index = ct_pros.index.str.replace(" ", "_")

	for index, expression in ct_pros.iterrows():
		ct_name = index

		# Get a baseline control for the ct I'm dropping out
		ct_pros_control, ct_pros_experiment = get_cont_and_exp_profiles(ct_pros, ct_name)
		# Save Save the control experiment 
		ct_pros_control.to_csv(f"{ct_name}_removed_control.csv")
		ct_pros_experiment.to_csv(f"{ct_name}_removed_experiment.csv")

def get_cont_and_exp_profiles(ct_pros, ct_name):
	"""Get a control background of ct profiles for my profile of interest
	 by removing ct profiles that strongly correlate to my profile of interest

	Args:
		ct_pros (_type_): _description_
		ct_name (_type_): _description_
	return:
		ct_pros_control (pd.Dataframe): Control dataframe that contains only ct profiles that don't strongly correlate with our profile of interest. also contains our profile of interest
		ct_pros_experiment (pd.Dataframe): Same as control dataframe, but with our ct profile of interest dropped out
	"""
	# Correlate
	cor=ct_pros.T.corr()
	# make a mask for what CTs I nwat to keep
	filter = cor.loc[:,ct_name] <= 0.7
	# Keep only some cts
	kept_cts = cor.loc[filter].index
	# keep only those cts,
 	# the ct_pros without our ct is technically our dropout experiment so return that as well
	ct_pros_experiment =ct_pros[ct_pros.index.isin(kept_cts)]
	# Add back in the ct of interest
	ct_pro_of_interest = ct_pros[ct_pros.index==ct_name]
	ct_pros_control = pd.concat([ct_pros_experiment,ct_pro_of_interest], axis = 0)
	# the ct_pros without our ct is technically our dropout experiment so return that as well
	return ct_pros_control, ct_pros_experiment

if __name__ == "__main__":
	main()