#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 20 } )
import numpy as np
import pandas as pd



# fill in appropriate filenames and paths here
path_to_native_3ay4 = "/Users/Research/antibody_project/send_to_louis/project_structs/native_crystal_struct_3ay4_Fc_FcgRIII.pdb"
path_to_low_E_double_pack_min_native = "/Users/Research/antibody_project/send_to_louis/project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb"
path_to_low_E_double_pack_min_native_no_Fc_glycan = "/Users/Research/antibody_project/send_to_louis/project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII_removed_Fc_sugar.pdb"



## LCM 50 data
# pseudo_interface_energy no ramp
path_to_pseudo_interface_energy_LCM_50_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_LCM_50_no_ramp.csv"
pseudo_interface_energy_LCM_50_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_LCM_50_no_ramp )
pseudo_interface_energy_LCM_50_no_ramp_data = pseudo_interface_energy_LCM_50_no_ramp_data.sort( "glycan_rmsd" )

# pseudo_interface_energy with ramp
path_to_pseudo_interface_energy_LCM_50_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_LCM_50_with_ramp.csv"
pseudo_interface_energy_LCM_50_with_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_LCM_50_with_ramp )
pseudo_interface_energy_LCM_50_with_ramp_data = pseudo_interface_energy_LCM_50_with_ramp_data.sort( "pseudo_interface_energy" )

# pseudo_interface_energy no ramp
path_to_pseudo_interface_energy_LCM_50_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_LCM_50_no_ramp.csv"
pseudo_interface_energy_LCM_50_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_LCM_50_no_ramp )
pseudo_interface_energy_LCM_50_no_ramp_data = pseudo_interface_energy_LCM_50_no_ramp_data.sort( "glycan_rmsd" )

# total_score no ramp
path_to_total_score_LCM_50_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_LCM_50_no_ramp.csv"
total_score_LCM_50_no_ramp_data = pd.read_csv( path_to_total_score_LCM_50_no_ramp )
total_score_LCM_50_no_ramp_data = total_score_LCM_50_no_ramp_data.sort( "total_score" )

# total_score with ramp
path_to_total_score_LCM_50_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_LCM_50_with_ramp.csv"
total_score_LCM_50_with_ramp_data = pd.read_csv( path_to_total_score_LCM_50_with_ramp )
total_score_LCM_50_with_ramp_data = total_score_LCM_50_with_ramp_data.sort( "total_score" )



## LCM 100 data
# pseudo_interface_energy no ramp
path_to_pseudo_interface_energy_LCM_100_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_LCM_100_no_ramp.csv"
pseudo_interface_energy_LCM_100_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_LCM_100_no_ramp )
pseudo_interface_energy_LCM_100_no_ramp_data = pseudo_interface_energy_LCM_100_no_ramp_data.sort( "glycan_rmsd" )

# total_score no ramp
path_to_total_score_LCM_100_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_LCM_100_no_ramp.csv"
total_score_LCM_100_no_ramp_data = pd.read_csv( path_to_total_score_LCM_100_no_ramp )
total_score_LCM_100_no_ramp_data = total_score_LCM_100_no_ramp_data.sort( "total_score" )

# pseudo_interface_energy random reset no ramp
path_to_pseudo_interface_energy_LCM_100_random_reset_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_LCM_100_random_reset_no_ramp.csv"
pseudo_interface_energy_LCM_100_random_reset_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_LCM_100_random_reset_no_ramp )
pseudo_interface_energy_LCM_100_random_reset_no_ramp_data = pseudo_interface_energy_LCM_100_random_reset_no_ramp_data.sort( "glycan_rmsd" )

# total_score random reset no ramp
path_to_total_score_LCM_100_random_reset_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_LCM_100_random_reset_no_ramp.csv"
total_score_LCM_100_random_reset_no_ramp_data = pd.read_csv( path_to_total_score_LCM_100_random_reset_no_ramp )
total_score_LCM_100_random_reset_no_ramp_data = total_score_LCM_100_random_reset_no_ramp_data.sort( "total_score" )



## GRM 50 data
# pseudo_interface_energy
path_to_pseudo_interface_energy_GRM_50 = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_GRM_50.csv"
pseudo_interface_energy_GRM_50_data = pd.read_csv( path_to_pseudo_interface_energy_GRM_50 )
pseudo_interface_energy_GRM_50_data = pseudo_interface_energy_GRM_50_data.sort( "pseudo_interface_energy", ascending = False )

# total_score
path_to_total_score_GRM_50 = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_GRM_50.csv"
total_score_GRM_50_data = pd.read_csv( path_to_total_score_GRM_50 )
total_score_GRM_50_data = total_score_GRM_50_data.sort( "total_score", ascending = False )




## SmallMover 50 data
# pseudo_interface_energy no ramp
path_to_pseudo_interface_energy_small_moves_50_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_small_moves_50_no_ramp.csv"
pseudo_interface_energy_small_moves_50_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_small_moves_50_no_ramp )
pseudo_interface_energy_small_moves_50_no_ramp_data = pseudo_interface_energy_small_moves_50_no_ramp_data.sort( "glycan_rmsd", ascending=False )

# pseudo_interface_energy with ramp
path_to_pseudo_interface_energy_small_moves_50_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_small_moves_50_with_ramp.csv"
pseudo_interface_energy_small_moves_50_with_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_small_moves_50_with_ramp )
pseudo_interface_energy_small_moves_50_with_ramp_data = pseudo_interface_energy_small_moves_50_with_ramp_data.sort( "pseudo_interface_energy" )

# total_score no ramp
path_to_total_score_small_moves_50_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_small_moves_50_no_ramp.csv"
total_score_small_moves_50_no_ramp_data = pd.read_csv( path_to_total_score_small_moves_50_no_ramp )
total_score_small_moves_50_no_ramp_data = total_score_small_moves_50_no_ramp_data.sort( "total_score" )

# total_score with ramp
path_to_total_score_small_moves_50_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_small_moves_50_with_ramp.csv"
total_score_small_moves_50_with_ramp_data = pd.read_csv( path_to_total_score_small_moves_50_with_ramp )
total_score_small_moves_50_with_ramp_data = total_score_small_moves_50_with_ramp_data.sort( "total_score" )




## SugarSmallMover 50 data
# pseudo_interface_energy random reset constraints on Gal with ramp (normal) and angle_multiplier of 2
path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_am2_with_ramp.csv"
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp )
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data.sort( "glycan_rmsd" )

# pseudo_interface_energy random reset constraints on Gal with ramp (normal), angle_multiplier of 2, and double hbond score term
path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_am2_double_hbond_with_ramp.csv"
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp )
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data.sort( "glycan_rmsd" )

# pseudo_interface_energy random reset constraints on Gal with ramp (normal)
path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp.csv"
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp )
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data.sort( "glycan_rmsd" )

# pseudo_interface_energy random reset constraints on Gal no ramp (normal)
path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_no_ramp.csv"
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp )
pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data.sort( "pseudo_interface_energy", ascending = False )

# pseudo_interface_energy random reset no omega no ramp
path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_no_omega_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_sugar_small_moves_50_random_reset_no_omega_no_ramp.csv"
pseudo_interface_energy_sugar_small_moves_50_random_reset_no_omega_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_no_omega_no_ramp )
pseudo_interface_energy_sugar_small_moves_50_random_reset_no_omega_no_ramp_data = pseudo_interface_energy_sugar_small_moves_50_random_reset_no_omega_no_ramp_data.sort( "pseudo_interface_energy", ascending = False )

# pseudo_interface_energy random reset with omega no ramp
path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_with_omega_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_sugar_small_moves_50_random_reset_with_omega_no_ramp.csv"
pseudo_interface_energy_sugar_small_moves_50_random_reset_with_omega_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_sugar_small_moves_50_random_reset_with_omega_no_ramp )
pseudo_interface_energy_sugar_small_moves_50_random_reset_with_omega_no_ramp_data = pseudo_interface_energy_sugar_small_moves_50_random_reset_with_omega_no_ramp_data.sort( "glycan_rmsd", ascending=False )

# total_score random reset constraints on Gal with ramp (normal) and angle_multiplier of 2
path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_am2_with_ramp.csv"
total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data = pd.read_csv( path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp )
total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data = total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data.sort( "glycan_rmsd" )

# total_score random reset constraints on Gal with ramp (normal), angle_multiplier of 2, and double hbond score term
path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_am2_double_hbond_with_ramp.csv"
total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data = pd.read_csv( path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp )
total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data = total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data.sort( "glycan_rmsd" )

# total_score random reset constraints on Gal with ramp (normal)
path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp.csv"
total_score_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data = pd.read_csv( path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp )
total_score_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data = total_score_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data.sort( "total_score", ascending = False )

# total_score random reset constraints on Gal no ramp (normal)
path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_no_ramp.csv"
total_score_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data = pd.read_csv( path_to_total_score_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp )
total_score_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data = total_score_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data.sort( "total_score", ascending = False )

# total_score random reset no omega no ramp
path_to_total_score_sugar_small_moves_50_random_reset_no_omega_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_sugar_small_moves_50_random_reset_no_omega_no_ramp.csv"
total_score_sugar_small_moves_50_random_reset_no_omega_no_ramp_data = pd.read_csv( path_to_total_score_sugar_small_moves_50_random_reset_no_omega_no_ramp )
total_score_sugar_small_moves_50_random_reset_no_omega_no_ramp_data = total_score_sugar_small_moves_50_random_reset_no_omega_no_ramp_data.sort( "total_score" )

# total_score random reset with omega no ramp
path_to_total_score_sugar_small_moves_50_random_reset_with_omega_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_sugar_small_moves_50_random_reset_with_omega_no_ramp.csv"
total_score_sugar_small_moves_50_random_reset_with_omega_no_ramp_data = pd.read_csv( path_to_total_score_sugar_small_moves_50_random_reset_with_omega_no_ramp )
total_score_sugar_small_moves_50_random_reset_with_omega_no_ramp_data = total_score_sugar_small_moves_50_random_reset_with_omega_no_ramp_data.sort( "total_score" )

# glycan_to_protein_contacts_Frecovered random reset constraints on Gal with ramp (normal)
path_to_glycan_to_protein_contacts_Frecovered_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/glycan_to_protein_contacts_Frecovered_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp.csv"
glycan_to_protein_contacts_Frecovered_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp_data = pd.read_csv( path_to_glycan_to_protein_contacts_Frecovered_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp )
glycan_to_protein_contacts_Frecovered_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp_data = glycan_to_protein_contacts_Frecovered_vs_glycan_rmsd_sugar_small_moves_50_random_reset_Gal_cst_with_ramp_data.sort( "glycan_to_protein_contacts_Frecovered" )




## RMSD vs Score pseudo_interface_energy
# LCM
LCM_50_rmsd_list_pseudo_interface_energy_data_no_ramp = pseudo_interface_energy_LCM_50_no_ramp_data[ "glycan_rmsd" ]
LCM_50_score_list_pseudo_interface_energy_data_no_ramp = pseudo_interface_energy_LCM_50_no_ramp_data[ "pseudo_interface_energy" ]
LCM_50_rmsd_list_pseudo_interface_energy_data_with_ramp = pseudo_interface_energy_LCM_50_with_ramp_data[ "glycan_rmsd" ]
LCM_50_score_list_pseudo_interface_energy_data_with_ramp = pseudo_interface_energy_LCM_50_with_ramp_data[ "pseudo_interface_energy" ]
LCM_100_rmsd_list_pseudo_interface_energy_data_no_ramp = pseudo_interface_energy_LCM_100_no_ramp_data[ "glycan_rmsd" ]
LCM_100_score_list_pseudo_interface_energy_data_no_ramp = pseudo_interface_energy_LCM_100_no_ramp_data[ "pseudo_interface_energy" ]
LCM_100_rmsd_list_pseudo_interface_energy_data_random_reset_no_ramp = pseudo_interface_energy_LCM_100_random_reset_no_ramp_data[ "glycan_rmsd" ]
LCM_100_score_list_pseudo_interface_energy_data_random_reset_no_ramp = pseudo_interface_energy_LCM_100_random_reset_no_ramp_data[ "pseudo_interface_energy" ]


# GRM
GRM_50_rmsd_list_pseudo_interface_energy_data = pseudo_interface_energy_GRM_50_data[ "glycan_rmsd" ]
GRM_50_score_list_pseudo_interface_energy_data = pseudo_interface_energy_GRM_50_data[ "pseudo_interface_energy" ]


# SmallMover
small_moves_50_rmsd_list_pseudo_interface_energy_data_no_ramp = pseudo_interface_energy_small_moves_50_no_ramp_data[ "glycan_rmsd" ]
small_moves_50_score_list_pseudo_interface_energy_data_no_ramp = pseudo_interface_energy_small_moves_50_no_ramp_data[ "pseudo_interface_energy" ]
small_moves_50_rmsd_list_pseudo_interface_energy_data_with_ramp = pseudo_interface_energy_small_moves_50_with_ramp_data[ "glycan_rmsd" ]
small_moves_50_score_list_pseudo_interface_energy_data_with_ramp = pseudo_interface_energy_small_moves_50_with_ramp_data[ "pseudo_interface_energy" ]


# SugarSmallMover
sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_with_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_with_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data[ "pseudo_interface_energy" ]
sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data[ "pseudo_interface_energy" ]
sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_with_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_with_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data[ "pseudo_interface_energy" ]
sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_no_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_no_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data[ "pseudo_interface_energy" ]
sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_no_omega_no_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_no_omega_no_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_no_omega_no_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_no_omega_no_ramp_data[ "pseudo_interface_energy" ]
sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_with_omega_no_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_with_omega_no_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_with_omega_no_ramp = pseudo_interface_energy_sugar_small_moves_50_random_reset_with_omega_no_ramp_data[ "pseudo_interface_energy" ]



## RMSD vs Score total_score
# LCM
LCM_50_rmsd_list_total_score_data_no_ramp = total_score_LCM_50_no_ramp_data[ "glycan_rmsd" ]
LCM_50_score_list_total_score_data_no_ramp = total_score_LCM_50_no_ramp_data[ "total_score" ]
LCM_50_rmsd_list_total_score_data_with_ramp = total_score_LCM_50_with_ramp_data[ "glycan_rmsd" ]
LCM_50_score_list_total_score_data_with_ramp = total_score_LCM_50_with_ramp_data[ "total_score" ]
LCM_100_rmsd_list_total_score_data_no_ramp = total_score_LCM_100_no_ramp_data[ "glycan_rmsd" ]
LCM_100_score_list_total_score_data_no_ramp = total_score_LCM_100_no_ramp_data[ "total_score" ]
LCM_100_rmsd_list_total_score_data_random_reset_no_ramp = total_score_LCM_100_random_reset_no_ramp_data[ "glycan_rmsd" ]
LCM_100_score_list_total_score_data_random_reset_no_ramp = total_score_LCM_100_random_reset_no_ramp_data[ "total_score" ]


# GRM
GRM_50_rmsd_list_total_score_data = total_score_GRM_50_data[ "glycan_rmsd" ]
GRM_50_score_list_total_score_data = total_score_GRM_50_data[ "total_score" ]


# SmallMover
small_moves_50_rmsd_list_total_score_data_no_ramp = total_score_small_moves_50_no_ramp_data[ "glycan_rmsd" ]
small_moves_50_score_list_total_score_data_no_ramp = total_score_small_moves_50_no_ramp_data[ "total_score" ]
small_moves_50_rmsd_list_total_score_data_with_ramp = total_score_small_moves_50_with_ramp_data[ "glycan_rmsd" ]
small_moves_50_score_list_total_score_data_with_ramp = total_score_small_moves_50_with_ramp_data[ "total_score" ]


# SugarSmallMover
sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_am2_with_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_am2_with_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_with_ramp_data[ "total_score" ]
sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_am2_double_hbond_with_ramp_data[ "total_score" ]
sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_with_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_with_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_with_ramp_data[ "total_score" ]
sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_no_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_no_ramp = total_score_sugar_small_moves_50_random_reset_Gal_constraints_no_ramp_data[ "total_score" ]
sugar_small_moves_50_rmsd_list_total_score_data_random_reset_no_omega_no_ramp = total_score_sugar_small_moves_50_random_reset_no_omega_no_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_total_score_data_random_reset_no_omega_no_ramp = total_score_sugar_small_moves_50_random_reset_no_omega_no_ramp_data[ "total_score" ]
sugar_small_moves_50_rmsd_list_total_score_data_random_reset_with_omega_no_ramp = total_score_sugar_small_moves_50_random_reset_with_omega_no_ramp_data[ "glycan_rmsd" ]
sugar_small_moves_50_score_list_total_score_data_random_reset_with_omega_no_ramp = total_score_sugar_small_moves_50_random_reset_with_omega_no_ramp_data[ "total_score" ]





################
#### PLOT 1 ####
################
# LCM and SmallMoves 50 no ramp and fa_intra_rep scorefunction

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( LCM_50_rmsd_list_pseudo_interface_energy_data_no_ramp, LCM_50_score_list_pseudo_interface_energy_data_no_ramp )
sub1.set_title( "LCM-50" )
sub1.set_xlabel( "glycan RMSD" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [-8, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( small_moves_50_rmsd_list_pseudo_interface_energy_data_no_ramp, small_moves_50_score_list_pseudo_interface_energy_data_no_ramp )
sub2.set_title( "SmallMoves-50" )
sub2.set_xlabel( "glycan RMSD" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "pseudo_interface_energy" )
sub2.set_ylim( [-8, 0 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( LCM_50_rmsd_list_total_score_data_no_ramp, LCM_50_score_list_total_score_data_no_ramp )
sub3.set_title( "LCM-50" )
sub3.set_xlabel( "glycan RMSD" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "total_score" )
sub3.set_ylim( [ -470, -420 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( small_moves_50_rmsd_list_total_score_data_no_ramp, small_moves_50_score_list_total_score_data_no_ramp )
sub4.set_title( "SmallMoves-50" )
sub4.set_xlabel( "glycan RMSD" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ -470, -420 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing LCM-50 and SmallMoves-50 with no ramp protocols with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 2 ####
################
# LCM and SmallMoves 50 with ramp and fa_intra_rep scorefunction

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( LCM_50_rmsd_list_pseudo_interface_energy_data_with_ramp, LCM_50_score_list_pseudo_interface_energy_data_with_ramp )
sub1.set_title( "LCM-50" )
sub1.set_xlabel( "glycan RMSD" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [-8, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( small_moves_50_rmsd_list_pseudo_interface_energy_data_with_ramp, small_moves_50_score_list_pseudo_interface_energy_data_with_ramp )
sub2.set_title( "SmallMoves-50" )
sub2.set_xlabel( "glycan RMSD" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "pseudo_interface_energy" )
sub2.set_ylim( [-8, 0 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( LCM_50_rmsd_list_total_score_data_with_ramp, LCM_50_score_list_total_score_data_with_ramp )
sub3.set_title( "LCM-50" )
sub3.set_xlabel( "glycan RMSD" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "total_score" )
sub3.set_ylim( [ -480, -430 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( small_moves_50_rmsd_list_total_score_data_with_ramp, small_moves_50_score_list_total_score_data_with_ramp )
sub4.set_title( "SmallMoves-50" )
sub4.set_xlabel( "glycan RMSD" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [-480, -430 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SmallMoves-50 and LCM-50 with REVERSE ramp protocols with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 3 ####
################
# LCM 50 vs 100 with fa_intra_rep scorefunction no ramp

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( LCM_50_rmsd_list_pseudo_interface_energy_data_no_ramp, LCM_50_score_list_pseudo_interface_energy_data_no_ramp )
sub1.set_title( "LCM-50", fontsize = 26 )
sub1.set_xlabel( "glycan RMSD", fontsize = 24 )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy", fontsize = 24 )
sub1.set_ylim( [ -7, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( LCM_100_rmsd_list_pseudo_interface_energy_data_no_ramp, LCM_100_score_list_pseudo_interface_energy_data_no_ramp )
sub2.set_title( "LCM-100", fontsize = 26 )
sub2.set_xlabel( "glycan RMSD", fontsize = 24 )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "pseudo_interface_energy", fontsize = 24 )
sub2.set_ylim( [ -7, 0 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( LCM_50_rmsd_list_total_score_data_no_ramp, LCM_50_score_list_total_score_data_no_ramp )
sub3.set_title( "LCM-50", fontsize = 26 )
sub3.set_xlabel( "glycan RMSD", fontsize = 24 )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "total_score", fontsize = 24 )
sub3.set_ylim( [ -475, -430 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( LCM_100_rmsd_list_total_score_data_no_ramp, LCM_100_score_list_total_score_data_no_ramp )
sub4.set_title( "LCM-100", fontsize = 26 )
sub4.set_xlabel( "glycan RMSD", fontsize = 24 )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score", fontsize = 24 )
sub4.set_ylim( [ -475, -430 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing LCM-50 and LCM-100 with no ramp protocols with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 4 ####
################
# GRM 50 and fa_intra_rep scorefunction

fig = plt.figure(figsize=(38, 12))
sub1 = plt.subplot(1, 2, 1)
sub1.scatter( GRM_50_rmsd_list_pseudo_interface_energy_data, GRM_50_score_list_pseudo_interface_energy_data )
sub1.set_title( "GRM-50", fontsize = 28 )
sub1.set_xlabel( "glycan RMSD", fontsize = 24 )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy", fontsize = 24 )
sub1.set_ylim( [-6, 0 ] )

sub2 = plt.subplot(1, 2, 2)
sub2.scatter( GRM_50_rmsd_list_total_score_data, GRM_50_score_list_total_score_data )
sub2.set_title( "GRM-50", fontsize = 28 )
sub2.set_xlabel( "glycan RMSD", fontsize = 24 )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score", fontsize = 24 )
sub2.set_ylim( [ -500, -350 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing GRM-50 protocol with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 34 )
plt.subplots_adjust(top=0.83)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 5 ####
################
# SugarSmallMoves 50 and random reset with or without omega torsion

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_no_omega_no_ramp, sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_no_omega_no_ramp )
sub1.set_title( "SugarSmallMoves-50 no omega" )
sub1.set_xlabel( "glycan RMSD" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [-7, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_with_omega_no_ramp, sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_with_omega_no_ramp )
sub2.set_title( "SugarSmallMoves-50 with omega" )
sub2.set_xlabel( "glycan RMSD" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "pseudo_interface_energy" )
sub2.set_ylim( [-7, 0 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( sugar_small_moves_50_rmsd_list_total_score_data_random_reset_no_omega_no_ramp, sugar_small_moves_50_score_list_total_score_data_random_reset_no_omega_no_ramp )
sub3.set_title( "SugarSmallMoves-50 no omega" )
sub3.set_xlabel( "glycan RMSD" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "total_score" )
sub3.set_ylim( [ -500, 0 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( sugar_small_moves_50_rmsd_list_total_score_data_random_reset_with_omega_no_ramp, sugar_small_moves_50_score_list_total_score_data_random_reset_with_omega_no_ramp )
sub4.set_title( "SugarSmallMoves-50 with omega" )
sub4.set_xlabel( "glycan RMSD" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [-500, 0 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves with random reset, no ramp, and with-or-without omega protocol with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 6 ####
################
# LCM-100 with and without random reset

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( LCM_100_rmsd_list_pseudo_interface_energy_data_random_reset_no_ramp, LCM_100_score_list_pseudo_interface_energy_data_random_reset_no_ramp )
sub1.set_title( "LCM-100 with random reset" )
sub1.set_xlabel( "glycan RMSD" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [-7, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( LCM_100_rmsd_list_pseudo_interface_energy_data_no_ramp, LCM_100_score_list_pseudo_interface_energy_data_no_ramp )
sub2.set_title( "LCM-100 no random reset" )
sub2.set_xlabel( "glycan RMSD" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "pseudo_interface_energy" )
sub2.set_ylim( [-7, 0 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( LCM_100_rmsd_list_total_score_data_random_reset_no_ramp, LCM_100_score_list_total_score_data_random_reset_no_ramp )
sub3.set_title( "LCM-100 with random reset" )
sub3.set_xlabel( "glycan RMSD" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "total_score" )
sub3.set_ylim( [ -480, -360 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( LCM_100_rmsd_list_total_score_data_no_ramp, LCM_100_score_list_total_score_data_no_ramp )
sub4.set_title( "LCM-100 no random reset" )
sub4.set_xlabel( "glycan RMSD" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [-470, -440 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing LCM-100 with no ramp, and with-or-without random reset protocol with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 7 ####
################
# SugarSmallMoves-50 random reset with (normal) ramp and constraints on Gal

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_with_ramp, sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_with_ramp )
sub1.set_title( "SugarSmallMoves-50 with ramp" )
sub1.set_xlabel( "glycan RMSD" )
#sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [-7, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_no_ramp, sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_no_ramp )
sub2.set_title( "SugarSmallMoves-50 no ramp" )
sub2.set_xlabel( "glycan RMSD" )
#sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "pseudo_interface_energy" )
#sub2.set_ylim( [-7, 0 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_with_ramp, sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_with_ramp )
sub3.set_title( "SugarSmallMoves-50 with ramp" )
sub3.set_xlabel( "glycan RMSD" )
#sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "total_score" )
sub3.set_ylim( [ -500, 0 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_no_ramp, sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_no_ramp )
sub4.set_title( "SugarSmallMoves-50 no ramp" )
sub4.set_xlabel( "glycan RMSD" )
#sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [-500, 0 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 with random reset, constraints on Gal, and with-or-without normal ramp protocol with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 8 ####
################
# SugarSmallMoves-50 random reset with (normal) ramp constraints on Gal and angle_multiplier of 2

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_with_ramp, sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_with_ramp )
sub1.set_title( "SugarSmallMoves-50 normal hbond weight" )
sub1.set_xlabel( "glycan RMSD" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [-8, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( sugar_small_moves_50_rmsd_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp, sugar_small_moves_50_score_list_pseudo_interface_energy_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp )
sub2.set_title( "SugarSmallMoves-50 *2 hbond weight" )
sub2.set_xlabel( "glycan RMSD" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "pseudo_interface_energy" )
sub2.set_ylim( [-15, 0 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_am2_with_ramp, sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_am2_with_ramp )
sub3.set_title( "SugarSmallMoves-50 normal hbond weight" )
sub3.set_xlabel( "glycan RMSD" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "total_score" )
sub3.set_ylim( [ -500, 0 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( sugar_small_moves_50_rmsd_list_total_score_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp, sugar_small_moves_50_score_list_total_score_data_random_reset_Gal_constraints_am2_double_hbond_with_ramp )
sub4.set_title( "SugarSmallMoves-50 *2 hbond weight" )
sub4.set_xlabel( "glycan RMSD" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [-900, 0 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 with random reset, constraints on Gal, normal ramp, angle_multiplier of 2, and with hbond * 1 or * 2 protocol with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 26 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )
