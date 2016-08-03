#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 20 } )
import numpy as np
from scipy.stats import linregress
import pandas as pd
from collections import Counter
from math import floor, ceil



###################################################
#### SSM-50 data using full Fc glycan no reset ####
###################################################
###########
#### 1 ####
###########
# am2, 7mpt, with ramp
path_to_using_native_full_glycan_no_reset_am2_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_am2_7_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_am2_7_mpt_with_ramp )
using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data = using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
sc = plt.scatter( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "total_score" ] )
plt.colorbar(sc)
ymin = floor( min( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 2 ] )
plt.ylabel( "pseudo_interface_energy" )
#plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
sc = plt.scatter( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.colorbar(sc)
ymin = floor( min( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 2 ] )
plt.ylabel( "total_score" )
#plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
sc = plt.scatter( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "total_score" ] )
plt.colorbar(sc)
ymin = floor( min( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 60 ] )
plt.ylabel( "pseudo_interface_energy" )
#plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
sc = plt.scatter( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "total_score" ], c=using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.colorbar(sc)
ymin = floor( min( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 60 ] )
plt.ylabel( "total_score" )
#plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, with ramp, am2, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )

#print "mean:", np.mean( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.total_score ), "\t(am2, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_no_reset_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 7mpt, ramp)"



###########
#### 2 ####
###########
# am2, 7mpt, with ramp, Gal only cst with 2A flexibility
path_to_using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.sort( "total_score")
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.total_score ), "\t(am2, 7mpt, Gal 2A, ramp)"
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 7mpt, Gal 2A, ramp)"

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, Gal only cst with 2A flexibility, with ramp, am2, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 3 ####
###########
# am2, 7mpt, with ramp, tight Gal only cst with 0.5A flexibility
path_to_using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.sort( "total_score")
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.total_score ), "\t(am2, 7mpt, Gal tight 0.5A, ramp)"
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 7mpt, Gal tight 0.5A, ramp)"

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_tight_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, tight Gal only cst with 0,5A flexibility, with ramp, am2, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 4 ####
###########
# am2, 7mpt, with ramp, Gal and its GlcNAc cst with 2A flexibility
path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.sort( "total_score")
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.total_score ), "\t(am2, 7mpt, Gal and its GlcNAc 2A, ramp)"
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 7mpt, Gal and its GlcNAc 2A, ramp)"

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, Gal and its GlcNAc cst with 2A flexibility, with ramp, am2, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 5 ####
###########
# am2, 7mpt, with ramp, Gal and its GlcNAc_tight tight cst with 0.5A flexibility
path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.sort( "total_score")
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.total_score ), "\t(am2, 7mpt, Gal and its GlcNAc 0.5A, ramp)"
#print "mean:", np.mean( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 7mpt, Gal and its GlcNAc 0.5A, ramp)"

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am2_7_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, Gal and its GlcNAc tight cst with 0,5A flexibility, with ramp, am2, 7 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )





###################################################
#### SSM-50 data using full Fc glycan no reset ####
###################################################
###########
#### 1 ####
###########
# am3, 5mpt, with ramp
path_to_using_native_full_glycan_no_reset_am3_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_am3_5_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_am3_5_mpt_with_ramp )
using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data = using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, with ramp, am3, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 2 ####
###########
# am3, 5mpt, with ramp, Gal only cst with 2A flexibility
path_to_using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, Gal only cst with 2A flexibility, with ramp, am3, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 3 ####
###########
# am3, 5mpt, with ramp, tight Gal only cst with 0.5A flexibility
path_to_using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_tight_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, tight Gal only cst with 0,5A flexibility, with ramp, am3, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 4 ####
###########
# am3, 5mpt, with ramp, Gal and its GlcNAc cst with 2A flexibility
path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, Gal and its GlcNAc cst with 2A flexibility, with ramp, am3, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )



###########
#### 5 ####
###########
# am3, 5mpt, with ramp, Gal and its GlcNAc_tight tight cst with 0.5A flexibility
path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp )
using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data = using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_Gal_and_its_GlcNAc_tight_am3_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, Gal and its GlcNAc tight cst with 0,5A flexibility, with ramp, am3, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )





###################################################
#### SSM-50 data using full Fc glycan no reset ####
###################################################
###########
#### 1 ####
###########
# am1, 5mpt, with ramp
path_to_using_native_full_glycan_no_reset_am1_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_no_reset_am1_5_mpt_with_ramp.csv"
using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_no_reset_am1_5_mpt_with_ramp )
using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data = using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data.sort( "total_score")

fig = plt.figure(figsize=(30, 15))
sub1 = plt.subplot(2, 2, 1)
ymin = floor( min( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub1.scatter( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub1.set_xlabel( "Fc_glycan_rmsd" )
sub1.set_xlim( [ 0, 2 ] )
sub1.set_ylabel( "pseudo_interface_energy" )
#sub1.set_ylim( [ ymin, ymax ] )

sub2 = plt.subplot(2, 2, 2)
ymin = floor( min( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data.total_score, 80 ) )
sub2.scatter( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "total_score" ] )
sub2.set_xlabel( "Fc_glycan_rmsd" )
sub2.set_xlim( [ 0, 2 ] )
sub2.set_ylabel( "total_score" )
#sub2.set_ylim( [ ymin, ymax ] )

sub3 = plt.subplot(2, 2, 3)
ymin = floor( min( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sub3.scatter( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
sub3.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub3.set_xlim( [ 100, 0 ] )
sub3.set_ylabel( "pseudo_interface_energy" )
#sub3.set_ylim( [ ymin, ymax ] )

sub4 = plt.subplot(2, 2, 4)
ymin = floor( min( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data.total_score, 80 ) )
sub4.scatter( using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_no_reset_am1_5_mpt_with_ramp_data[ "total_score" ] )
sub4.set_xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
sub4.set_xlim( [ 100, 0 ] )
sub4.set_ylabel( "total_score" )
#sub4.set_ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and no reset, with ramp, am1, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )





######################################################
#### SSM-50 data using full Fc glycan light reset ####
######################################################
###########
#### 1 ####
###########
# am2, 5mpt, with ramp, light reset
path_to_using_native_full_glycan_light_reset_am2_5_mpt_with_ramp = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_light_reset_am2_5_mpt_with_ramp.csv"
using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data = pd.read_csv( path_to_using_native_full_glycan_light_reset_am2_5_mpt_with_ramp )
using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data = using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.sort( "pseudo_interface_energy")

fig = plt.figure(figsize=(30, 15))
plt.subplot( 221 )
ymin = floor( min( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 222 )
ymin = floor( min( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ 0, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 223 )
ymin = floor( min( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.pseudo_interface_energy, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 224 )
ymin = floor( min( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "total_score" ] ) )
ymax = ceil( np.percentile( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.total_score, 80 ) )
sc = plt.scatter( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 100, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using native and SSM-50 on Fc glycan and light reset, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )

#print "mean:", np.mean( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.total_score ), "\t(am2, 5mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_light_reset_am2_5_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am2, 5mpt, ramp)"





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

#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.total_score ), "\t(am1, 7mpt, ramp)"
#print "mean:", np.mean( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "median:", np.median( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "stdev:", np.std( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "var:", np.var( using_native_full_glycan_LCM_reset_Gal_and_its_GlcNAc_tight_am1_7_mpt_with_ramp_data.Fc_glycan_rmsd ), "\t(am1, 7mpt, ramp)"
