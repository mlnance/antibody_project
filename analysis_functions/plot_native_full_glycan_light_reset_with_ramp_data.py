#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 20 } )
import numpy as np
from scipy.stats import linregress
import pandas as pd
from collections import Counter
from math import floor, ceil



######################################################
#### SSM-50 data using full Fc glycan light reset ####
######################################################
###########
#### 1 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 10-20)
path_to_using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (10-20), with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_10_20_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 2 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 5-10), Man3 to F241 and GlcNAc6 to F243 cst
path_to_using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (5-10), Man3 to F241 and GlcNAc6 to F241 cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_5_10_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 3 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 10-15), Man3 to F241 and GlcNAc6 to F243 cst
path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.sort( "delta_interface_sasa" )
#print using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ ["filename", "Fc_glycan_rmsd", "delta_interface_sasa" ] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "fa_sol" ] )
#plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (10-15), Man3 to F241 and GlcNAc6 to F241 cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 4 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 10-15), Man3 to F241 and GlcNAc6 to F243 cst, double fa_elec sf
path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.sort( "atom_pair_constraint" )

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (10-15), Man3 to F241 and GlcNAc6 to F241 cst, with ramp, double fa_elec sf, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_double_fa_elec_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 5 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 10-15), Man3 to F241 and GlcNAc6 to F243 cst, half fa_sol sf
path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.sort( "atom_pair_constraint" )

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (10-15), Man3 to F241 and GlcNAc6 to F241 cst, with ramp, half fa_sol sf, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 6 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 10-15), Man3 to F241 and GlcNAc6 to F243 cst, third fa_sol sf
path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy" )
print_these = using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ (using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd <= 7.5 ) & (using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd >= 6 ) ]
#print_these = print_these.sort( "pseudo_interface_energy" )
#print print_these[ ["filename", "Fc_glycan_rmsd", "pseudo_interface_energy"] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "fa_rep" ] )
#plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (10-15), Man3 to F241 and GlcNAc6 to F241 cst, with ramp, third fa_sol sf, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_10_15_Man3_to_F241_and_GlcNAc6_to_F243_third_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 7 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 20-25), Man3 to F241 and GlcNAc6 to F243 cst, half fa_sol sf
path_to_using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_20_25_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy" )
print_these = using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ (using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd <= 7.5 ) & (using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd >= 6 ) ]
#print_these = print_these.sort( "pseudo_interface_energy" )
#print print_these[ ["filename", "Fc_glycan_rmsd", "pseudo_interface_energy"] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "fa_rep" ] )
#plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (20-25), Man3 to F241 and GlcNAc6 to F241 cst, with ramp, half fa_sol sf, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_20_25_Man3_to_F241_and_GlcNAc6_to_F243_half_fa_sol_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 8 ####
###########
# am2, 5mpt, with ramp, light reset (+/- 10-15), two_glycan_to_protein_atoms_tightest
path_to_using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_20_25_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy" )
print_these = using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ (using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd <= 7.5 ) & (using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd >= 6 ) ]
#print_these = print_these.sort( "pseudo_interface_energy" )
#print print_these[ ["filename", "Fc_glycan_rmsd", "pseudo_interface_energy"] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd_after_packmin" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset (10-15), two glycan to protein atoms tightest cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_10_15_two_glycan_to_protein_atoms_tightest_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"
