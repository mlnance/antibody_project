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
# pseduo_interface_energy no ramp
path_to_pseudo_interface_energy_LCM_50_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_LCM_50_no_ramp.csv"
pseudo_interface_energy_LCM_50_no_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_LCM_50_no_ramp )
pseudo_interface_energy_LCM_50_no_ramp_data = pseudo_interface_energy_LCM_50_no_ramp_data.sort( "pseudo_interface_energy" )

# pseduo_interface_energy with ramp
path_to_pseudo_interface_energy_LCM_50_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/pseudo_interface_energy_vs_glycan_rmsd_LCM_50_with_ramp.csv"
pseudo_interface_energy_LCM_50_with_ramp_data = pd.read_csv( path_to_pseudo_interface_energy_LCM_50_with_ramp )
pseudo_interface_energy_LCM_50_with_ramp_data = pseudo_interface_energy_LCM_50_with_ramp_data.sort( "pseudo_interface_energy" )

# total_score no ramp
path_to_total_score_LCM_50_no_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_LCM_50_with_ramp.csv"
total_score_LCM_50_no_ramp_data = pd.read_csv( path_to_total_score_LCM_50_no_ramp )
total_score_LCM_50_no_ramp_data = total_score_LCM_50_no_ramp_data.sort( "total_score" )

# total_core with ramp
path_to_total_score_LCM_50_with_ramp = "/Users/Research/pyrosetta_dir/metrics_vs_glycan_rmsd_data/total_score_vs_glycan_rmsd_LCM_50_no_ramp.csv"
total_score_LCM_50_with_ramp_data = pd.read_csv( path_to_total_score_LCM_50_with_ramp )
total_score_LCM_50_with_ramp_data = total_score_LCM_50_with_ramp_data.sort( "total_score" )



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
pseudo_interface_energy_small_moves_50_no_ramp_data = pseudo_interface_energy_small_moves_50_no_ramp_data.sort( "pseudo_interface_energy" )

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





## RMSD vs Score pseudo_interface_energy
# LCM
LCM_rmsd_list_pseudo_interface_data_no_ramp = pseudo_interface_energy_LCM_50_no_ramp_data[ "glycan_rmsd" ]
LCM_score_list_pseudo_interface_data_no_ramp = pseudo_interface_energy_LCM_50_no_ramp_data[ "pseudo_interface_energy" ]
LCM_rmsd_list_pseudo_interface_data_with_ramp = pseudo_interface_energy_LCM_50_with_ramp_data[ "glycan_rmsd" ]
LCM_score_list_pseudo_interface_data_with_ramp = pseudo_interface_energy_LCM_50_with_ramp_data[ "pseudo_interface_energy" ]


# GRM
GRM_rmsd_list_pseudo_interface_data = pseudo_interface_energy_GRM_50_data[ "glycan_rmsd" ]
GRM_score_list_pseudo_interface_data = pseudo_interface_energy_GRM_50_data[ "pseudo_interface_energy" ]


# SmallMover
small_moves_rmsd_list_pseudo_interface_data_no_ramp = pseudo_interface_energy_small_moves_50_no_ramp_data[ "glycan_rmsd" ]
small_moves_score_list_pseudo_interface_data_no_ramp = pseudo_interface_energy_small_moves_50_no_ramp_data[ "pseudo_interface_energy" ]
small_moves_rmsd_list_pseudo_interface_data_with_ramp = pseudo_interface_energy_small_moves_50_with_ramp_data[ "glycan_rmsd" ]
small_moves_score_list_pseudo_interface_data_with_ramp = pseudo_interface_energy_small_moves_50_with_ramp_data[ "pseudo_interface_energy" ]



## RMSD vs Score total_score
# LCM
LCM_rmsd_list_total_score_data_no_ramp = total_score_LCM_50_no_ramp_data[ "glycan_rmsd" ]
LCM_score_list_total_score_data_no_ramp = total_score_LCM_50_no_ramp_data[ "total_score" ]
LCM_rmsd_list_total_score_data_with_ramp = total_score_LCM_50_with_ramp_data[ "glycan_rmsd" ]
LCM_score_list_total_score_data_with_ramp = total_score_LCM_50_with_ramp_data[ "total_score" ]


# GRM
GRM_rmsd_list_total_score_data = total_score_GRM_50_data[ "glycan_rmsd" ]
GRM_score_list_total_score_data = total_score_GRM_50_data[ "total_score" ]


# SmallMover
small_moves_rmsd_list_total_score_data_no_ramp = total_score_small_moves_50_no_ramp_data[ "glycan_rmsd" ]
small_moves_score_list_total_score_data_no_ramp = total_score_small_moves_50_no_ramp_data[ "total_score" ]
small_moves_rmsd_list_total_score_data_with_ramp = total_score_small_moves_50_with_ramp_data[ "glycan_rmsd" ]
small_moves_score_list_total_score_data_with_ramp = total_score_small_moves_50_with_ramp_data[ "total_score" ]




################
#### PLOT 1 ####
################

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( LCM_rmsd_list_pseudo_interface_data_no_ramp, LCM_score_list_pseudo_interface_data_no_ramp )
sub1.set_title( "3ay4-type glycosylation with LCM-50 no ramp" )
sub1.set_xlabel( "glycan RMSD" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [-8, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( LCM_rmsd_list_total_score_data_no_ramp, LCM_score_list_total_score_data_no_ramp )
sub2.set_title( "3ay4-type glycosylation with LCM-50 no ramp" )
sub2.set_xlabel( "glycan RMSD" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ -470, -430 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( small_moves_rmsd_list_pseudo_interface_data_no_ramp, small_moves_score_list_pseudo_interface_data_no_ramp )
sub3.set_title( "3ay4-type glycosylation with SmallMoves-50 no ramp" )
sub3.set_xlabel( "glycan RMSD" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [-10, 0 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( small_moves_rmsd_list_total_score_data_no_ramp, small_moves_score_list_total_score_data_no_ramp )
sub4.set_title( "3ay4-type glycosylation with SmallMoves-50 no ramp" )
sub4.set_xlabel( "glycan RMSD" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ -500, -400 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing no ramp protocols with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )




################
#### PLOT 2 ####
################

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( LCM_rmsd_list_pseudo_interface_data_with_ramp, LCM_score_list_pseudo_interface_data_with_ramp )
sub1.set_title( "3ay4-type glycosylation with LCM-50 with ramp" )
sub1.set_xlabel( "glycan RMSD" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [-10, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( LCM_rmsd_list_total_score_data_with_ramp, LCM_score_list_total_score_data_with_ramp )
sub2.set_title( "3ay4-type glycosylation with LCM-50 with ramp" )
sub2.set_xlabel( "glycan RMSD" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ -470, -430 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( small_moves_rmsd_list_pseudo_interface_data_with_ramp, small_moves_score_list_pseudo_interface_data_with_ramp )
sub3.set_title( "3ay4-type glycosylation with SmallMoves-50 with ramp" )
sub3.set_xlabel( "glycan RMSD" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [-8, 0 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( small_moves_rmsd_list_total_score_data_with_ramp, small_moves_score_list_total_score_data_with_ramp )
sub4.set_title( "3ay4-type glycosylation with SmallMoves-50 with ramp" )
sub4.set_xlabel( "glycan RMSD" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [-480, -430 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing with ramp protocols with metrics vs glycan rmsd"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
