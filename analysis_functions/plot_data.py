#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 20 } )
import numpy as np
from scipy.stats import linregress
import pandas as pd
from collections import Counter
from math import floor, ceil


def rounddown( x ):
    return int( floor( x / 5 ) * 5 )

def roundup( x ):
    return int( ceil( x / 5 ) * 5 )



############################
#### Base Pack/Min Pose ####
############################
###########
#### 1 ####
###########
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



###########
#### 2 ####
###########
## using fa_intra_rep 0.440 sf
path_to_double_pack_and_min_only_of_native_using_fa_intra_rep_sf = "/Users/Research/pyrosetta_dir/metric_data/double_pack_and_min_only_of_native_using_fa_intra_rep_sf.csv"
double_pack_and_min_only_of_native_using_fa_intra_rep_sf_data = pd.read_csv( path_to_double_pack_and_min_only_of_native_using_fa_intra_rep_sf )

fig = plt.figure(figsize=(17, 8))
plt.scatter( double_pack_and_min_only_of_native_using_fa_intra_rep_sf_data[ "rmsd" ], double_pack_and_min_only_of_native_using_fa_intra_rep_sf_data[ "total_score" ] )
plt.xlabel( "rmsd" )
#plt.xlim( [ 0.5, 1.1 ] )
plt.ylabel( "total_score" )
#plt.ylim( [-1070, -1010 ] )

# save the plot
plt.tight_layout()
plot_title = "Double pack and min of native 3ay4 into low energy structure using fa_intra_rep sf"
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
ymin = floor( min( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_constraints_am3_5_mpt_data.pseudo_interface_energy, 80 ) )
sub1.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_constraints_am3_5_mpt_data.total_score, 80 ) )
sub2.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] )
sub2.set_xlabel( "glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_constraints_am3_5_mpt_data.pseudo_interface_energy, 80 ) )
sub3.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_constraints_am3_5_mpt_data.total_score, 80 ) )
sub4.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 with random reset, ramp, tight Gal constraints, angle multiplier 3, 5 mpt"
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
ymin = floor( min( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data.pseudo_interface_energy, 80 ) )
sub1.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "total_score" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data.total_score, 80 ) )
sub2.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "total_score" ] )
sub2.set_xlabel( "glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data.pseudo_interface_energy, 80 ) )
sub3.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub3.set_xlim( [ 90, 20 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "total_score" ] ) )
ymax = ceil( np.percentile( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data.total_score, 80 ) )
sub4.scatter( random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_and_its_GlcNAc_constraints_am3_triple_hbond_5_mpt_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub4.set_xlim( [ 90, 20 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 with random reset, ramp, tight Gal and its GlcNAc constraints, triple hbond, angle multiplier 3, 5 mpt"
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
ymin = floor( min( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub1.set_xlim( [ 100, 0 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub2.set_xlim( [ 100, 0 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "glycan_rmsd" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "glycan_rmsd" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "glycan_rmsd" ], native_reset_tight_Gal_and_its_GlcNAc_constraints_am1_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "glycan_rmsd" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 with native reset, ramp, tight Gal and its GlcNAc constraints, angle multiplier 1, 5 mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





############################################
#### SSM-50 data using native and reset ####
############################################
###########
#### 1 ####
###########
# random reset, Gal and its GlcNAc tight, triple hbond sf, am 3, 5 moves per trial, with ramp data
path_to_using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt.csv"
using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data = pd.read_csv( path_to_using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt )
using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data = using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data.sort( "Fc_glycan_rmsd" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data.total_score, 80 ) )
sub2.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am3_5_mpt_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose with ramp, random reset, tight Gal and its GlcNAc constraints, triple hbond sf weights, angle multiplier 3, and 5 mpt"
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
ymin = floor( min( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data.total_score, 80 ) )
sub2.scatter( using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_random_reset_Gal_and_its_GlcNAc_tight_triple_hbond_am1_5_mpt_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose with ramp, random reset, tight Gal and its GlcNAc constraints, triple hbond sf weights, angle multiplier 1, and 5 mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





##################################################
#### SSM-100 data using native Fc branch only ####
##################################################
###########
#### 1 ####
###########

# using native, random reset of Fc glycan branch only, Gal and its GlcNAc high and tight constraints, am3, 3 mpt, with ramp
path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials.csv"
using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data = pd.read_csv( path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ "delta_pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["delta_pseudo_interface_energy"] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "delta_pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ "delta_total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_total_score, 80 ) )
sub2.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["delta_total_score"] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "delta_total_score" )
sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ "delta_pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["delta_pseudo_interface_energy"] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "delta_pseudo_interface_energy" )
sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data[ "delta_total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data.delta_total_score, 80 ) )
sub4.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_am2_3mpt_100trials_data["delta_total_score"] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "delta_total_score" )
sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-100 on native with random reset of Fc glycan branch only, Gal and its GlcNAc high and tight cst, am2, 3mpt, 100 trials, with ramp"
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
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["pseudo_interface_energy"] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.total_score, 80 ) )
sub2.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["total_score"] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["pseudo_interface_energy"] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data.total_score, 80 ) )
sub4.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_triple_hbond_am2_3mpt_100trials_data["total_score"] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-100 on native with random reset of Fc glycan branch only, Gal and its GlcNAc high and tight cst, triple hbond sf, am2, 3mpt, 100 trials, with ramp"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 3 ####
###########
# using native, random reset of Fc glycan branch only, Gal and its GlcNAc high and tight constraints, double fa_atr and half fa_rep sf, am3, 3 mpt, with ramp
path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials.csv"
using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data = pd.read_csv( path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials )
using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data.sort( "delta_total_score" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["pseudo_interface_energy"] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data.total_score, 80 ) )
sub2.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["total_score"] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["pseudo_interface_energy"] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data.total_score, 80 ) )
sub4.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data["total_score"] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-100 on native with random reset of Fc glycan branch only, Gal and its GlcNAc high and tight cst, double fa_atr and half fa_rep sf, am2, 3mpt, 100 trials, with ramp"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 4 ####
###########
# using native, random reset of Fc glycan branch only, Gal and its GlcNAc high and tight constraints, double fa_elec sf, am3, 3 mpt, with ramp
path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials = "/Users/Research/pyrosetta_dir/metric_data/using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials.csv"
using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data = pd.read_csv( path_to_using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials )
using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data.sort( "total_score" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["pseudo_interface_energy"] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data.total_score, 80 ) )
sub2.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["Fc_glycan_rmsd"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["total_score"] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["pseudo_interface_energy"] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data.total_score, 80 ) )
sub4.scatter( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"], using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_elec_am2_3mpt_100trials_data["total_score"] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-100 on native with random reset of Fc glycan branch only, Gal and its GlcNAc high and tight cst, double fa_elec sf, am2, 3mpt, 100 trials, with ramp"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





####################
#### histograms ####
####################

#using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data = using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data.sort( "biggest_score_diff_by_scoretype" )
#data = list( using_native_random_reset_only_branch_glycan_Gal_and_its_GlcNAc_high_and_tight_double_fa_atr_half_fa_rep_am2_3mpt_100trials_data.biggest_score_diff_scoretype )
#counted_data = Counter( data )





############################################################
#### SSM-50 data using Fc glycan branch native no reset ####
############################################################
###########
#### 1 ####
###########
# Fc branch only, Gal and its GlcNAc high and tight, am 1, 3 moves per trial, with ramp data
path_to_using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp.csv"
using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data = pd.read_csv( path_to_using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp )
using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data = using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose Fc branch only with ramp, high and tight Gal and its GlcNAc constraints, angle multiplier 1, and 3 moves per turn"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )





##########################################################
#### SSM-50 data using full Fc glycan native no reset ####
##########################################################
###########
#### 1 ####
###########
# full glycan (except core GlcNAc), am1, 3mpt, high and tight Gal and GlcNAc cst, no reset, with ramp
path_to_using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp.csv"
using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data = pd.read_csv( path_to_using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp )
using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data = using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_high_and_tight_am1_3_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose no reset, with ramp, high and tight Gal and GlcNAc cst, am1, 3mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 2 ####
###########
# full glycan (except core GlcNAc), am1, 3mpt, tight Gal and GlcNAc cst, no reset, with ramp
path_to_using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp.csv"
using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data = pd.read_csv( path_to_using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp )
using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data = using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose no reset, with ramp, tight Gal and GlcNAc cst, am1, 3mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 3 ####
###########
# full glycan (except core GlcNAc), am1, 3mpt, Gal and GlcNAc cst, no reset, with ramp
path_to_using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp.csv"
using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data = pd.read_csv( path_to_using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp )
using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data = using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_tight_am1_3_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose no reset, with ramp, Gal and GlcNAc cst, am1, 3mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 4 ####
###########
# full glycan (except core GlcNAc), am2, 5mpt, tight Gal and GlcNAc cst, no reset, with ramp
path_to_using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp.csv"
using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp )
using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data = using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose no reset, with ramp, tight Gal and GlcNAc cst, am2, 5mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 5 ####
###########
# full glycan (except core GlcNAc), am2, 5mpt, Gal and GlcNAc cst, no reset, with ramp
path_to_using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp.csv"
using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp )
using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data = using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose no reset, with ramp, Gal and GlcNAc cst, am2, 5mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 6 ####
###########
# full glycan (except core GlcNAc), am1, 3mpt, NO ramp or reset
path_to_using_native_am1_3_mpt = "/Users/Research/pyrosetta_dir/metric_data/using_native_am1_3_mpt.csv"
using_native_am1_3_mpt_data = pd.read_csv( path_to_using_native_am1_3_mpt )
using_native_am1_3_mpt_data = using_native_am1_3_mpt_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_am1_3_mpt_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_am1_3_mpt_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_am1_3_mpt_data[ "Fc_glycan_rmsd" ], using_native_am1_3_mpt_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 10 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_am1_3_mpt_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_am1_3_mpt_data.total_score, 80 ) )
sub2.scatter( using_native_am1_3_mpt_data[ "Fc_glycan_rmsd" ], using_native_am1_3_mpt_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 10 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 glycan comparing SSM-50 on native Pose no ramp or reset, angle multiplier 1, and 3 mpt"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )
