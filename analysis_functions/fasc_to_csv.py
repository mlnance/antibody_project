#!/usr/bin/python
__author__ = "morganlnance"



import argparse

parser = argparse.ArgumentParser(description="Use Rosetta to calculate RMSD between a native pose and a directory of structures")
parser.add_argument("fasc_file", type=str, help="the path to the .fasc file with the relevant data.")
parser.add_argument("resulting_filename", type=str, help="what do you want the resulting csv file to be called? This program will add the .csv extension for you")
parser.add_argument("utility_dir", type=str, help="the path to the utility functions directory.")
input_args = parser.parse_args()




#################
#### IMPORTS ####
#################

import os, sys, csv
try:
    import pandas
except:
    pass

# try to add the utility directory
try:
    sys.path.append( input_args.utility_dir )
except:
    print "It seems like you gave me an incorrect path to the utility directory, exiting"
    sys.exit()    
from util import read_fasc_file, get_score_term_from_fasc_data_dict



#########################
#### DATA EXTRACTION ####
#########################

# metrics data gets held in these lists in a dictionary
filename = []
total_score = []
delta_total_score = []
rmsd = []
fa_atr = []
fa_rep = []
fa_sol = []
fa_intra_rep = []
fa_elec = []
pro_close = []
hbond_sr_bb = []
hbond_lr_bb = []
hbond_bb_sc = []
hbond_sc = []
dslf_fa13 = []
atom_pair_constraint = []
rama = []
omega = []
fa_dun = []
p_aa_pp = []
yhh_planarity = []
ref = []
sugar_bb = []
glycan_rmsd = []
Fc_glycan_rmsd = []
pseudo_interface_energy = []
delta_pseudo_interface_energy = []
delta_biggest_score_diff_by_scoretype = []
biggest_score_diff_scoretype = []
std_interface_interaction_energy = []
delta_std_interface_interaction_energy = []
hbonds = []
delta_hbonds = []
glycan_to_protein_contacts = []
delta_glycan_to_protein_contacts = []
Fc_glycan_to_protein_Fnat_recovered_contacts = []
Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A = []
Fc_glycan_to_FcR_glycan_contacts = []
delta_Fc_glycan_to_FcR_glycan_contacts = []
Fc_glycan_to_FcR_glycan_Fnat_recovered_contacts = []
interface_res_contacts_8_A = []
delta_interface_res_contacts_8_A = []
interface_sasa = []
delta_interface_sasa = []
GlcNAc_to_its_Phe_contacts_5A = []
Fc_glycan_to_FcR_glycan_Fnat_tot_contacts_recovered_10A = []
Fc_glycan_sasa_contributed = []
delta_Fc_glycan_sasa_contributed = []
MonteCarlo_acceptance_rate = []

metrics = { "filename":filename, "total_score":total_score, "delta_total_score":delta_total_score, "rmsd":rmsd, "fa_atr":fa_atr, "fa_rep":fa_rep, "fa_sol":fa_sol, "fa_intra_rep":fa_intra_rep, "fa_elec":fa_elec, "pro_close":pro_close, "hbond_sr_bb":hbond_sr_bb, "hbond_lr_bb":hbond_lr_bb, "hbond_bb_sc":hbond_bb_sc, "hbond_sc":hbond_sc, "dslf_fa13":dslf_fa13, "atom_pair_constraint":atom_pair_constraint, "rama":rama, "omega":omega, "fa_dun":fa_dun, "p_aa_pp":p_aa_pp, "yhh_planarity":yhh_planarity, "ref":ref, "sugar_bb":sugar_bb, "glycan_rmsd":glycan_rmsd, "Fc_glycan_rmsd":Fc_glycan_rmsd, "pseudo_interface_energy":pseudo_interface_energy, "delta_pseudo_interface_energy":delta_pseudo_interface_energy, "std_interface_interaction_energy":std_interface_interaction_energy, "delta_biggest_score_diff_by_scoretype":delta_biggest_score_diff_by_scoretype, "biggest_score_diff_scoretype":biggest_score_diff_scoretype, "delta_std_interface_interaction_energy":delta_std_interface_interaction_energy, "hbonds":hbonds, "delta_hbonds":delta_hbonds, "glycan_to_protein_contacts":glycan_to_protein_contacts, "delta_glycan_to_protein_contacts":delta_glycan_to_protein_contacts, "Fc_glycan_to_protein_Fnat_recovered_contacts":Fc_glycan_to_protein_Fnat_recovered_contacts, "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A, "Fc_glycan_to_FcR_glycan_contacts":Fc_glycan_to_FcR_glycan_contacts, "delta_Fc_glycan_to_FcR_glycan_contacts":delta_Fc_glycan_to_FcR_glycan_contacts, "Fc_glycan_to_FcR_glycan_Fnat_recovered_contacts":Fc_glycan_to_FcR_glycan_Fnat_recovered_contacts, "interface_res_contacts_8_A":interface_res_contacts_8_A, "delta_interface_res_contacts_8_A":delta_interface_res_contacts_8_A, "interface_sasa":interface_sasa, "delta_interface_sasa":delta_interface_sasa, "Fc_glycan_to_FcR_glycan_Fnat_tot_contacts_recovered_10A":Fc_glycan_to_FcR_glycan_Fnat_tot_contacts_recovered_10A, "GlcNAc_to_its_Phe_contacts_5A":GlcNAc_to_its_Phe_contacts_5A, "Fc_glycan_sasa_contributed":Fc_glycan_sasa_contributed, "delta_Fc_glycan_sasa_contributed":delta_Fc_glycan_sasa_contributed, "MonteCarlo_acceptance_rate":MonteCarlo_acceptance_rate }



## check the .fasc file
working_dir = os.getcwd() + '/'
if not os.path.isfile( input_args.fasc_file ):
    print "It seems like you gave me an incorrect path to the .fasc file, exiting"
    sys.exit()
else:
    fasc_file = input_args.fasc_file

# get the fasc_data_dict
fasc_data_dict = read_fasc_file( fasc_file )

# add the name of the metric to the top of each list
for metric_name in metrics:
    metrics[ metric_name ].append( metric_name )

# loop over each decoy and pull out desired data
used_metric_names = [ "filename" ]
unused_metric_names = []
for decoy in fasc_data_dict.keys():
    metrics[ "filename" ].append( decoy )

    decoy_metrics = fasc_data_dict[ decoy ]
    for metric_name in metrics.keys():
        if metric_name not in unused_metric_names:
            if metric_name is not "filename":
                try:
                    metrics[ metric_name ].append( decoy_metrics[ metric_name ] )
                    if metric_name not in used_metric_names:
                        used_metric_names.append( metric_name )
                except:
                    unused_metric_names.append( metric_name )

# check out the passed result filename for .csv extension
dump_filename = input_args.resulting_filename
if not dump_filename.endswith( ".csv" ):
    dump_filename += ".csv"

# write out the .csv data file
rows = zip( *[ metrics[ metric ] for metric in used_metric_names ] )
with open( dump_filename, "wb" ) as f:
    writer = csv.writer( f )
    writer.writerows( rows )
