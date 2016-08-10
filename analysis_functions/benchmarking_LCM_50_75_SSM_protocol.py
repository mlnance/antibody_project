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


def get_r_of_line_of_best_fit( data, metric, in_type = 1 ):
    """
    Return m, b, r, and p corresponding from linregress( data[metric], type( data[metric] ) )
    Default is a linear line of best fit of x to y
    Can do things like log10
    :param data: DataFrame
    :param metric: metric of interest
    :param in_type: math function
    """
    try:
        m, b, r, p, std_err = linregress( data[ "Fc_glycan_rmsd" ], [ in_type(y) for y in data[ metric ] ] )
    except:
        m, b, r, p, std_err = linregress( data[ "Fc_glycan_rmsd" ], [ y * in_type for y in data[ metric ] ] )

    return r


def print_r_squared_data( data, r_squared_dict, name ):
    r_squared_keys = r_squared_dict.keys()
    r_squared_keys.sort( reverse = True )
    
    for r_squared in r_squared_keys:
        print r_squared, r_squared_dict[ r_squared ], name


def print_other_data( data ):
    print "MC min:", min( data[ "MonteCarlo_acceptance_rate" ] ), "max:", max( data[ "MonteCarlo_acceptance_rate" ] ), "mean:", np.mean( data[ "MonteCarlo_acceptance_rate" ] ), "median:", np.median( data[ "MonteCarlo_acceptance_rate" ] )

    rmsd_data = data[ data["Fc_glycan_rmsd"] <= 2 ]
    rmsd_count = len( rmsd_data["Fc_glycan_rmsd"] )
    print "% Fc_glycan_rmsd <= 2:", round( ( float(rmsd_count) / float(len( data["Fc_glycan_rmsd"])) ) * 100, 2 )

    Fnat_data = data[ data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"] >= 80 ]
    Fnat_count = len( Fnat_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"] )
    print "% Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A >= 80:", round( ( float(Fnat_count) / float(len( data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"])) ) * 100, 2 )




low_E_native_data = pd.read_csv( "/Users/Research/pyrosetta_dir/metric_data/low_E_native_data.csv" )
low_E_native_fa_intra_rep_data = pd.read_csv( "/Users/Research/pyrosetta_dir/metric_data/low_E_native_fa_intra_rep_data.csv" )
low_E_native_fa_intra_rep_pseudo_interface_energy = -20.6597224902  # number 365


####################################################
#### SSM-50 data using full Fc glycan LCM reset ####
####################################################
###########
#### 1 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest, with ramp, fa_intra_rep low E native
path_to_using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_50_SSM_benchmark_am2_5_mpt_fa_intra_rep.csv"
using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data = pd.read_csv( path_to_using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep )
path_to_using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep = "/Users/Research/pyrosetta_dir/metric_data/using_native_full_glycan_50_SSM_no_LCM_reset_benchmark_am2_5_mpt_fa_intra_rep.csv"
using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data = pd.read_csv( path_to_using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep )

metrics = list( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
ymin = floor( min( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "pseudo_interface_energy" ] ) - 5 )
ymax = ceil( np.percentile( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data.pseudo_interface_energy, 80 ) + 5 )
x = 0.0
y = low_E_native_fa_intra_rep_pseudo_interface_energy # number 365
if floor(y) < ymin:
    ymin = floor(y - 5)
sc = plt.scatter( x, y, marker='o', s=24, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "pseudo_interface_energy" ], marker='v', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 322 )
ymin = floor( min( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "total_score" ] ) - 5 )
ymax = ceil( np.percentile( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data.total_score, 80 ) + 5 )
x = 0.0
y = min( low_E_native_fa_intra_rep_data[ "total_score" ] )
if floor(y) < ymin:
    ymin = floor(y - 5)
sc = plt.scatter( x, y, marker='o', s=24, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "total_score" ], marker='v', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "total_score" ], c=using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "atom_pair_constraint" ] )
#plt.colorbar(sc)
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 323 )
ymin = floor( min( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "pseudo_interface_energy" ] ) - 5 )
ymax = ceil( np.percentile( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data.pseudo_interface_energy, 80 ) + 5 )
x = 100.0
y = low_E_native_fa_intra_rep_pseudo_interface_energy # number 365
if floor(y) < ymin:
    ymin = floor(y - 5)
sc = plt.scatter( x, y, marker='o', s=24, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "pseudo_interface_energy" ], marker='v', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "pseudo_interface_energy" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 324 )
ymin = floor( min( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "total_score" ] ) - 5 )
ymax = ceil( np.percentile( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data.total_score, 80 ) + 5 )
x = 100.0
y = min( low_E_native_fa_intra_rep_data[ "total_score" ] )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_benchmark_50_SSM_no_LCM_reset_am2_5_mpt_fa_intra_rep_data[ "total_score" ], marker='v', s=20, c="orange", clip_on=False )
if floor(y) < ymin:
    ymin = floor(y - 5)
sc = plt.scatter( x, y, marker='o', s=24, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "total_score" ] )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ ymin, ymax ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "Benchmark fa_intra_rep 3ay4 using SSM-50 on Fc glycan and LCM reset, two glycan to protein atoms tightest cst, with ramp, am2, 5 mpt"
plt.suptitle( plot_title, fontsize = 24 )
plt.subplots_adjust(top=0.90)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data, log10_r_squared_to_metric_dict, "log10" )
print_other_data( using_native_full_glycan_benchmark_50_SSM_am2_5_mpt_fa_intra_rep_data )
print "50 25"
print
