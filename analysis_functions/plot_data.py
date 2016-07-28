#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 20 } )
import numpy as np
from scipy.stats import linregress
import pandas as pd



############################
#### Base Pack/Min Pose ####
############################

path_to_double_pack_and_min_only_of_native = "/Users/Research/pyrosetta_dir/metric_data/double_pack_and_min_only_of_native.csv"
double_pack_and_min_only_of_native_data = pd.read_csv( path_to_double_pack_and_min_only_of_native )

fig = plt.figure(figsize=(17, 8))
plt.scatter( double_pack_and_min_only_of_native_data[ "rmsd" ], double_pack_and_min_only_of_native_data[ "total_score" ] )
plt.xlabel( "rmsd" )
plt.xlim( [ 0.5, 1.1 ] )
plt.ylabel( "total_score" )
plt.ylim( [-1070, -1010 ] )

# save the plot
plt.tight_layout()
plot_title = "Double pack and min of native 3ay4 into low energy structure"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )




#################################
#### SugarSmallMoves-50 data ####
#################################

###########
#### 1 ####
###########
## random reset tight Gal cst angle multiplier of 3 and 5 moves per turn
path_to_random_reset_tight_Gal_constraints_am3_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/sugar_small_moves_50_random_reset_tight_Gal_constraints_am3_5_mpt.csv"
random_reset_tight_Gal_constraints_am3_5_mpt_data = pd.read_csv( path_to_random_reset_tight_Gal_constraints_am3_5_mpt )
random_reset_tight_Gal_constraints_am3_5_mpt_data = random_reset_tight_Gal_constraints_am3_5_mpt_data.sort( "glycan_rmsd" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] )
#sub1.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "atom_pair_constraint" ] )
sub1.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub1.set_xlim( [ 80, 20 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ -10, 0 ] )
#sub1.set_ylim( [ -7, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub2.set_xlim( [ 80, 20 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ -430, -380 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] )
#sub3.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "atom_pair_constraint" ] )
sub3.set_xlabel( "glycan_rmsd" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ 10, 30 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] )
sub4.set_xlabel( "glycan_rmsd" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ -430, -380 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 with random reset, ramp, tight Gal constraints, angle multiplier 3, 5 moves per turn"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





###########
#### 2 ####
###########
## random reset tight Gal and its GlcNAc cst, angle multiplier of 3, triple hbond, and 5 moves per turn
path_to_random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/sugar_small_moves_50_random_reset_tight_Gal_and_its_GlcNAc_am3_triple_hbond_5_mpt.csv"
random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data = pd.read_csv( path_to_random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt )
random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data = random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data.sort( "glycan_rmsd" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "pseudo_interface_energy" ] )
#sub1.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "atom_pair_constraint" ] )
sub1.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub1.set_xlim( [ 90, 20 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ -30, 0 ] )
#sub1.set_ylim( [ 0, 20 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub2.set_xlim( [ 90, 20 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ -1165, -1000 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "pseudo_interface_energy" ] )
#sub3.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "atom_pair_constraint" ] )
sub3.set_xlabel( "glycan_rmsd" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ -30, 0 ] )
#sub3.set_ylim( [ 0, 20 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "total_score" ] )
sub4.set_xlabel( "glycan_rmsd" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ -1165, -1000 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 with random reset, ramp, tight Gal and its GlcNAc constraints, angle multiplier 3, 5 moves per turn"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





###########
#### 3 ####
###########
## native reset, tight Gal and its GlcNAc, am1, 5 mpt, with ramp
path_to_native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/sugar_small_moves_50_native_reset_tight_Gal_and_its_GlcNAc_am1_5_mpt_with_ramp.csv"
native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data = pd.read_csv( path_to_native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp )
native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data = native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data.sort( "glycan_rmsd" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
#sub1.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
sub1.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub1.set_xlim( [ 100, 80 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ -30, 0 ] )
#sub1.set_ylim( [ 0, 20 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub2.set_xlim( [ 100, 80 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ -465, -440 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "glycan_rmsd" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
#sub3.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "glycan_rmsd" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
sub3.set_xlabel( "glycan_rmsd" )
sub3.set_xlim( [ 1.7, 2.4 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ -30, 0 ] )
#sub3.set_ylim( [ 0, 20 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "glycan_rmsd" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "glycan_rmsd" )
sub4.set_xlim( [ 1.7, 2.4 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ -465, -440 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 with native reset, ramp, tight Gal and its GlcNAc constraints, angle multiplier 1, 5 moves per turn"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





########################################################
#### SugarSmallMoves-50 data using native and reset ####
########################################################

###########
#### 1 ####
###########
# random reset, Gal and its GlcNAc tight, triple hbond sf, am 3, 5 moves per trial, with ramp data
path_to_using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt.csv"
using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data = pd.read_csv( path_to_using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt )
using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data = using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data.sort( "Fc_glycan_rmsd" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "delta_pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "delta_pseudo_interface_energy" )
sub1.set_ylim( [ 0, 50 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "delta_total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "delta_total_score" )
sub2.set_ylim( [ 0, 200 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 on the native Pose with ramp, tight Gal and its GlcNAc constraints, triple hbond sf weights, angle multiplier 3, and 5 moves per turn"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 2 ####
###########
# random reset, Gal and its GlcNAc tight, triple hbond sf, am 1, 5 moves per trial, with ramp data
path_to_using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt.csv"
using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data = pd.read_csv( path_to_using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt )
using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data = using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data.sort( "Fc_glycan_rmsd" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "delta_pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "delta_pseudo_interface_energy" )
sub1.set_ylim( [ 0, 50 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "delta_total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "delta_total_score" )
sub2.set_ylim( [ 0, 200 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 on the native Pose with ramp, tight Gal and its GlcNAc constraints, triple hbond sf weights, angle multiplier 1, and 5 moves per turn"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





##############################################################
#### SugarSmallMoves-100 data using native Fc branch only ####
##############################################################
###########
#### 1 ####
###########

# using native, random reset of Fc glycan branch only, Gal and its GlcNAc high and tight constraints, am3, 3 mpt, with ramp
path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials.csv"
using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data = pd.read_csv( path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_pseudo_interface_energy <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_pseudo_interface_energy, 80 ) ]
sub1.scatter( data["Fc_glycan_rmsd"], data["delta_pseudo_interface_energy"] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "delta_pseudo_interface_energy" )
#sub1.set_ylim( [ 0, 50 ] )

sub2 = plt.subplot(2, 2, 2)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_total_score <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_total_score, 80 ) ]
sub2.scatter( data["Fc_glycan_rmsd"], data["delta_total_score"] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "delta_total_score" )
#sub2.set_ylim( [ 0, 200 ] )

sub3 = plt.subplot(2, 2, 3)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_pseudo_interface_energy <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_pseudo_interface_energy, 80 ) ]
sub3.scatter( data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], data["delta_pseudo_interface_energy"] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "delta_pseudo_interface_energy" )
#sub3.set_ylim( [ 0, 200 ] )

sub4 = plt.subplot(2, 2, 4)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_total_score <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_total_score, 80 ) ]
sub4.scatter( data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], data["delta_total_score"] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "delta_total_score" )
#sub4.set_ylim( [ 0, 200 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-100 on the native with random reset of Fc glycan branch only, Gal and its GlcNAc high and tight cst, am2, 5mpt, 100 trials, with ramp"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 2 ####
###########
# using native, random reset of Fc glycan branch only, Gal and its GlcNAc high and tight constraints, triple hbond sf, am3, 3 mpt, with ramp
path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials.csv"
using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data = pd.read_csv( path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_pseudo_interface_energy <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_pseudo_interface_energy, 80 ) ]
sub1.scatter( data["Fc_glycan_rmsd"], data["delta_pseudo_interface_energy"] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "delta_pseudo_interface_energy" )
#sub1.set_ylim( [ 0, 50 ] )

sub2 = plt.subplot(2, 2, 2)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_total_score <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_total_score, 80 ) ]
sub2.scatter( data["Fc_glycan_rmsd"], data["delta_total_score"] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "delta_total_score" )
#sub2.set_ylim( [ 0, 200 ] )

sub3 = plt.subplot(2, 2, 3)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_pseudo_interface_energy <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_pseudo_interface_energy, 80 ) ]
sub3.scatter( data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], data["delta_pseudo_interface_energy"] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "delta_pseudo_interface_energy" )
#sub3.set_ylim( [ 0, 200 ] )

sub4 = plt.subplot(2, 2, 4)
data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_total_score <= np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.delta_total_score, 80 ) ]
sub4.scatter( data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], data["delta_total_score"] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "delta_total_score" )
#sub4.set_ylim( [ 0, 200 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-100 on the native with random reset of Fc glycan branch only, Gal and its GlcNAc high and tight cst, triple hbond sf, am2, 5mpt, 100 trials, with ramp"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





#######################################################
#### SugarSmallMoves-50 data using native no reset ####
#######################################################

###########
#### 1 ####
###########
# Gal and its GlcNAc high and tight, am 1, 3 moves per trial, with ramp data
path_to_using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt = "/Users/Research/pyrosetta_dir/metric_data/using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt.csv"
using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data = pd.read_csv( path_to_using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt )
using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data = using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data.sort( "delta_total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
data = using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data[ using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data.delta_pseudo_interface_energy <= np.percentile( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data.delta_pseudo_interface_energy, 80 ) ]
sub1.scatter( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data[ "delta_pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 1 ] )
sub1.set_ylabel( "delta_pseudo_interface_energy" )
#sub1.set_ylim( [ 0, 50 ] )

sub2 = plt.subplot(2, 2, 2)
data = using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data[ using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data.delta_total_score <= np.percentile( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data.delta_total_score, 80 ) ]
sub2.scatter( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_data[ "delta_total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 1 ] )
sub2.set_ylabel( "delta_total_score" )
#sub2.set_ylim( [ 0, 200 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 on the native Pose with ramp, high and tight Gal and its GlcNAc constraints, angle multiplier 1, and 3 moves per turn"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )
