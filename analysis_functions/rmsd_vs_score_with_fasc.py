#!/usr/bin/python
__author__ = "morganlnance"


'''
run rmsd_vs_score_with_fasc.py pdb_copies_dont_touch/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb structures_from_louis/3ay4_Fc_FcgRIIIa/glycan_sampling_just_25_LCM/base_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII_removed_Fc_sugar/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII_removed_Fc_sugar_glycosylated_then_just_25_LCM.fasc test_fasc_out.csv ~/antibody_project/send_to_louis/project_utility_files/
'''


import argparse

parser = argparse.ArgumentParser(description="Use Rosetta to calculate RMSD between a native pose and a directory of structures")
parser.add_argument("fasc_file", type=str, help="the path to the .fasc file with the relevant data.")
parser.add_argument("metric_name", type=str, help="what is the name of the score metric in the .fasc file you want to compare against glycan_rmsd?")
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

# try to add the utility directory to the path to get pose_metrics_util
try:
    sys.path.append( input_args.utility_dir )
except:
    print "It seems like you gave me an incorrect path to the utility directory, exiting"
    sys.exit()    
from util import read_fasc_file, get_score_term_from_fasc_data_dict



#########################
#### DATA EXTRACTION ####
#########################

# check the .fasc file
working_dir = os.getcwd() + '/'
if not os.path.isfile( input_args.fasc_file ):
    print "It seems like you gave me an incorrect path to the .fasc file, exiting"
    sys.exit()
else:
    fasc_file = input_args.fasc_file

# get the fasc_data_dict
fasc_data_dict = read_fasc_file( fasc_file )

'''
# get the ligand_rmsd data from the .fasc file
glycan_rmsd_data = get_score_term_from_fasc_data_dict( fasc_data_dict, "glycan_rmsd" )
#glycan_rmsd = glycan_rmsd_data.values()

# get the total_score data from the .fasc file
tot_score_data = get_score_term_from_fasc_data_dict( fasc_data_dict, "total_score" )
#tot_score = tot_score_data.values()
'''

# check that the passed metric exists in the .fasc file
try:
    metric = get_score_term_from_fasc_data_dict( fasc_data_dict, input_args.metric_name )
except:
    print "\nIt appears that '%s' is not a valid metric in the given .fasc file. Please check your input.\n" %input_args.metric_name
    sys.exit()

# prepare a .csv data file
#header = [ "pdb_names", "pseudo_interface_energy", "glycan_rmsd" ]
#header = [ "pdb_names", "total_score", "glycan_rmsd" ]
header = [ "pdb_names", input_args.metric_name, "glycan_rmsd" ]
df_data = []
df_data.append( header )

# loop over each decoy and pull out desired data
for decoy_num in fasc_data_dict.decoy_nums:
    try:
        # get the pdb names
        pdb_name_temp = fasc_data_dict[ decoy_num ][ "filename" ]
        pdb_name = '_'.join( pdb_name_temp.split('/')[-1].split('_')[-4:] )

        # get the glycan_rmsd
        glycan_rmsd = fasc_data_dict[ decoy_num ][ "glycan_rmsd" ]

        # get the comparing metric
        #pseudo_interface_energy = fasc_data_dict[ decoy_num ][ "pseudo_interface_energy" ]
        #score = fasc_data_dict[ decoy_num ][ "total_score" ]
        metric = fasc_data_dict[ decoy_num ][ input_args.metric_name ]

        # add the data together to append to the data frame
        #data = [ pdb_name, pseudo_interface_energy, glycan_rmsd ]
        #data = [ pdb_name, score, glycan_rmsd ]
        data = [ pdb_name, metric, glycan_rmsd ]
        df_data.append( data )
        
    # skip the entry if there is an issue
    except:
        pass

# check out the passed result filename for .csv extension
filename = input_args.resulting_filename
if not filename.endswith( ".csv" ):
    filename += ".csv"

# write out the .csv data file
with open( filename, "wb" ) as f:
    writer = csv.writer( f )
    writer.writerows( df_data )
