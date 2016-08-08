#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 20 } )
import numpy as np
from scipy.stats import linregress, normaltest
import pandas as pd
from collections import Counter
from math import floor, ceil, log10



def normalize_data( data ):
    new_data = [ ( x - min( data ) / ( max( data ) - min( data ) ) ) for x in data ]

    return new_data



####################################################
#### SSM-50 data using full Fc glycan LCM reset ####
####################################################
###########
#### 1 ####
###########
# am2_5_mpt_half_fa_sol_glycan_to_protein_atoms, with ramp
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

m, b, r, p, std_err = linregress( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
print r**2

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, half fa_sol sf, glycan to protein atoms cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()



###########
#### 2 ####
###########
# am2_5_mpt_half_fa_sol_glycan_to_protein_atoms_tight, with ramp
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp )

using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data.sort( "atom_pair_constraint" )
#print using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ ["filename", "Fc_glycan_rmsd", "atom_pair_constraint" ] ]

m, b, r, p, std_err = linregress( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
print r**2

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, half fa_sol sf, glycan to protein atoms tight cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()



###########
#### 3 ####
###########
# am2_5_mpt_half_fa_sol_glycan_to_protein_atoms_tightest, with ramp
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp )

using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.sort( "atom_pair_constraint" )
#print using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ ["filename", "Fc_glycan_rmsd", "atom_pair_constraint" ] ]

m, b, r, p, std_err = linregress( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
print r**2

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, half fa_sol sf, glycan to protein atoms tightest cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()



###########
#### 4 ####
###########
# fa_intra_rep native, am2_5_mpt_half_fa_sol_glycan_to_protein_atoms_tightest, with ramp
path_to_using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp.csv"
using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp )

using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.sort( "atom_pair_constraint" )
#print using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ ["filename", "Fc_glycan_rmsd", "atom_pair_constraint" ] ]

m, b, r, p, std_err = linregress( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
print r**2

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_fa_intra_rep_native_full_glycan_LCM_reset_half_fa_sol_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using fa_intra_rep native and SSM-50 on Fc glycan and LCM reset, half fa_sol sf, glycan to protein atoms tightest cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()



###########
#### 5 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest, with ramp
path_to_using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp )

using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.sort( "atom_pair_constraint" )
#print using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ ["filename", "Fc_glycan_rmsd", "atom_pair_constraint" ] ]

metrics = list( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.columns.values )
r_squared_to_metric_dict = {}
for metric in metrics:
    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ metric ] )
        #if pval >= 0.05:
        m, b, r, p, std_err = linregress( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ metric ] )
        r_squared_to_metric_dict[ r**2 ] = metric

r_squared_keys = r_squared_to_metric_dict.keys()
r_squared_keys.sort( reverse = True )
for r_squared in r_squared_keys:
    print r_squared, r_squared_to_metric_dict[ r_squared ]
print

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan and LCM reset, two glycan to protein atoms tightest cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()



###########
#### 6 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest, with ramp
path_to_using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp )

using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data.sort( "atom_pair_constraint" )
#print using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ ["filename", "Fc_glycan_rmsd", "atom_pair_constraint" ] ]

'''
metrics = list( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data.columns.values )
r_squared_to_metric_dict = {}
for metric in metrics:
    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ metric ] )
        #if pval >= 0.05:
        m, b, r, p, std_err = linregress( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ metric ] )
        r_squared_to_metric_dict[ r**2 ] = metric
r_squared_keys = r_squared_to_metric_dict.keys()
r_squared_keys.sort( reverse = True )
for r_squared in r_squared_keys:
    print r_squared, r_squared_to_metric_dict[ r_squared ]
print
'''

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_two_glycan_to_protein_atoms_high_and_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan and LCM reset, two glycan to protein atoms high and tightest cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()
