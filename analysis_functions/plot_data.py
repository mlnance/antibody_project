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
plt.title( "Base Pack/Min Pose" )
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

#######
### Fc_glycan_to_protein_Fnat_recovered_contacts
#######

# pseudo_interface_energy vs Fc_glycan_to_protein_Fnat_recovered_contacts random reset tight Gal cst angle multiplier of 3 and 5 moves per turn
path_to_random_reset_tight_Gal_constraints_am3_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/sugar_small_moves_50_random_reset_tight_Gal_constraints_am3_5_mpt.csv"
random_reset_tight_Gal_constraints_am3_5_mpt_data = pd.read_csv( path_to_random_reset_tight_Gal_constraints_am3_5_mpt )
random_reset_tight_Gal_constraints_am3_5_mpt_data = random_reset_tight_Gal_constraints_am3_5_mpt_data.sort( "Fc_glycan_to_protein_Fnat_recovered_contacts" )

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
sub1.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] )
#sub1.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "atom_pair_constraint" ] )
sub1.set_title( "tight Gal only" )
sub1.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub1.set_xlim( [ 20, 80 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
sub1.set_ylim( [ -7, 0 ] )

sub2 = plt.subplot(2, 2, 2)
sub2.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "Fc_glycan_to_protein_Fnat_recovered_contacts" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] )
sub2.set_title( "tight Gal only" )
sub2.set_xlabel( "Fc_glycan_to_protein_Fnat_recovered_contacts" )
sub2.set_xlim( [ 20, 80 ] )
sub2.set_ylabel( "total_score" )
sub2.set_ylim( [ -430, -380 ] )

sub3 = plt.subplot(2, 2, 3)
sub3.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "pseudo_interface_energy" ] )
#sub3.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "atom_pair_constraint" ] )
sub3.set_title( "tight Gal only" )
sub3.set_xlabel( "glycan_rmsd" )
sub3.set_xlim( [ 0, 10 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
sub3.set_ylim( [ -7, 0 ] )

sub4 = plt.subplot(2, 2, 4)
sub4.scatter( random_reset_tight_Gal_constraints_am3_5_mpt_data[ "glycan_rmsd" ], random_reset_tight_Gal_constraints_am3_5_mpt_data[ "total_score" ] )
sub4.set_title( "tight Gal only" )
sub4.set_xlabel( "glycan_rmsd" )
sub4.set_xlim( [ 0, 10 ] )
sub4.set_ylabel( "total_score" )
sub4.set_ylim( [ -430, -380 ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4-type glycosylation comparing SugarSmallMoves-50 with random reset, ramp, angle multiplier 3, 5 moves per turn, and either tight Gal or tight Gal and its GlcNAc"
plt.suptitle( plot_title, fontsize = 30 )
plt.subplots_adjust(top=0.87)
plt.savefig( plot_title, dpi=120, transparent=True )
