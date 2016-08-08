#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 20 } )
import numpy as np
from scipy.stats import linregress
import pandas as pd
from collections import Counter
from math import floor, ceil



####################################################
#### SSM-50 data using full Fc glycan LCM reset ####
####################################################
###########
#### 1 ####
###########
# am2, 5mpt, with ramp
path_to_using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 2 ####
###########
# am2, 5mpt, with ramp, Gal 2A flexbility cst
path_to_using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Gal cst 2A flexibility, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 3 ####
###########
# am2, 5mpt, with ramp, Gal tight 0.5A flexbility cst
path_to_using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Gal tight cst with 0,5A flexibility, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 4 ####
###########
# am2, 5mpt, with ramp, Gal and its GlcNAc 2A flexbility cst
path_to_using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Gal and its GlcNAc cst 2A flexibility, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"



###########
#### 5 ####
###########
# am2, 5mpt, with ramp, tight Gal and its GlcNAc 0.5A flexbility cst
path_to_using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Gal and its GlcNAc tight cst with 0,5A flexibility, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"





####################################################
#### SSM-50 data using full Fc glycan LCM reset ####
####################################################
###########
#### 1 ####
###########
# am1, 7mpt, with ramp
path_to_using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp )
using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, am1, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 2 ####
###########
# am1, 7mpt, with ramp, Gal 2A flexbility cst
path_to_using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Gal cst 2A flexibility, am1, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 3 ####
###########
# am1, 7mpt, with ramp, Gal and its GlcNAc tight 0.5A flexbility cst
path_to_using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Gal and its GlcNAc tight cst 0,5A flexibility, am1, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"





####################################################
#### SSM-50 data using full Fc glycan LCM reset ####
####################################################
###########
#### 1 ####
###########
# am2, 5mpt, with ramp, Man3_to_F241_and_GlcNAc6_to_F243 1.0A flexibility cst
path_to_using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "delta_Fc_glycan_sasa_contributed" ] )
#plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Man3 to F241 and GlcNAc6 to F243 cst, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 2 ####
###########
# am2, 5mpt, with ramp, double hbond sf, Man3_to_F241_and_GlcNAc6_to_F243 high (0.5 stdev) 1.0A flexibility cst
path_to_using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, double hbond sf, Man3 to F241 and GlcNAc6 to F243 high cst, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_double_hbond_Man3_to_F241_and_GlcNAc6_to_F243_high_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 3 ####
###########
# am2, 5mpt, with ramp, half fa_sol sf, Man3_to_F241_and_GlcNAc6_to_F243 1.0A flexibility cst
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Man3 to F241 and GlcNAc6 to F243 cst, half fa_sol sf, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 4 ####
###########
# am2, 5mpt, with ramp, double fa_atr and half fa_sol sf, Man3_to_F241_and_GlcNAc6_to_F243 1.0A flexibility and Gal to chainA/B tight (0.5A flexibility) cst
path_to_using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Man3 to F241 and GlcNAc6 to F243 cst and Gal to protein tight, double fa_atr and half fa_sol sf, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_double_fa_atr_half_fa_sol_native_3ay4_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 5 ####
###########
# am2, 5mpt, with ramp, half fa_sol sf, Man3_to_F241_and_GlcNAc6_to_F243 1.0A flexibility and Gal to chainA/B tight (0.5A flexibility) cst
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Man3 to F241 and GlcNAc6 to F243 cst and Gal to protein tight, half fa_sol sf, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_chainA_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 6 ####
###########
# am2, 5mpt, with ramp, Man3_to_F241_and_GlcNAc6_to_F243 1.0A flexibility and Gal to Pro244-Pro245 (1.2A flexibility) cst
path_to_using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, Man3 to F241 and GlcNAc6 to F243 cst and Gal to Pro244, Pro245, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 7 ####
###########
# am2, 5mpt, with ramp, half fa_sol sf, Man3_to_F241_and_GlcNAc6_to_F243 1.0A flexibility and Gal to Pro244-Pro245 (1.2A flexibility) cst
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.sort( "Fc_glycan_rmsd")
#print using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ ["filename","Fc_glycan_rmsd","atom_pair_constraint", "pseudo_interface_energy"] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, half fa_sol sf, Man3 to F241 and GlcNAc6 to F243 cst and Gal to Pro244, Pro245, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 8 ####
###########
# am2, 5mpt, with ramp, half fa_sol sf, Man3_to_F241_and_GlcNAc6_to_F243 high and tight 0.5A flexibility and Gal to Pro244-Pro245 high and tight (0.5A flexibility) cst
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.sort( "Fc_glycan_rmsd")
#print using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ ["filename","Fc_glycan_rmsd","atom_pair_constraint" ] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "atom_pair_constraint" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, half fa_sol sf, Man3 to F241 and GlcNAc6 to F243 cst and Gal to Pro244, Pro245 high and tight, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Pro_245_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



###########
#### 9 ####
###########
# am2, 5mpt, with ramp, half fa_sol sf, Man3_to_F241_and_GlcNAc6_to_F243 high and tight 0.5A flexibility and Gal to Pro244-Glu258 high and tight (0.5A flexibility) cst
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.sort( "Fc_glycan_rmsd")
#print using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ ["filename","Fc_glycan_rmsd","atom_pair_constraint" ] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd_after_packmin" ] - using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd_after_reset" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, half fa_sol sf, Man3 to F241 and GlcNAc6 to F243 cst and Gal to Pro244, Glu258 high and tight, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"



############
#### 10 ####
############
# am2, 5mpt, with ramp, half fa_sol double atom_pair_constraint sf, Man3_to_F241_and_GlcNAc6_to_F243 high and tight 0.5A flexibility and Gal to Pro244-Glu258 high and tight (0.5A flexibility) cst
path_to_using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp )
using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data = using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.sort( "Fc_glycan_rmsd")
#print using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ ["filename","Fc_glycan_rmsd","atom_pair_constraint" ] ]

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
#sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd_after_packmin" ] - using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd_after_reset" ] )
plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and LCM reset, with ramp, half fa_sol double atom_pair_constraint sf, Man3 to F241 and GlcNAc6 to F243 cst and Gal to Pro244, Glu258 high and tight, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_half_fa_sol_double_apc_Man3_to_F241_and_GlcNAc6_to_F243_Gal_to_Pro_244_Glu_258_high_and_tight_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"
