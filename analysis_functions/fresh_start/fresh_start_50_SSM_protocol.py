#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 24 } )
import numpy as np
from scipy.stats import linregress, normaltest
import pandas as pd
from collections import Counter
from math import floor, ceil, log10
import sys



def normalize_data( data ):
    new_data = [ ( x - min( data ) / ( max( data ) - min( data ) ) ) for x in data ]

    return new_data


def get_r_of_line_of_best_fit( data, metric, in_type = 1 ):
    """
    Return m, b, r, and p corresponding from linregress( data[Fc_glycan_rmsd], type( data[metric] ) )
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


def get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( data, metric ):
    """
    """
    data_holder = {}
    stdev_data_holder = {}
    bin_size = 1

    for ii in np.arange( 0, int(ceil(max(data["Fc_glycan_rmsd"]))) + bin_size, bin_size ):
        # get data corresponding to decoys that fall within ii <= Fc_glycan_rmsd <= ii_plus_one
        binned_data = data[ ( data["Fc_glycan_rmsd"] >= ii ) & ( data["Fc_glycan_rmsd"] < (ii + 1) ) ]

        # add the average of the metric value for this rmsd bin
        if len( binned_data ) != 0:
            #data_holder[ '_'.join( [ str(ii), str(ii+1) ] ) ] = np.mean( binned_data[metric] )
            data_holder[ ii ] = np.mean( binned_data[metric] )
        # if there is no data (ie, the bin was empty), add None instead
        else:
            #data_holder[ '_'.join( [ str(ii), str(ii+1) ] ) ] = None
            data_holder[ ii ] = None

    # now for each bin that doesn't include None, get the r of the line for this data
    x_data = []
    y_data = []
    for x in np.arange( 0, int(ceil(max(data["Fc_glycan_rmsd"]))) + bin_size, bin_size ):
        # skip None entries
        if data_holder[x] is not None:
            x_data.append( x )
            y_data.append( data_holder[x] )
    m, b, r, p, std_err = linregress( x_data, y_data )
    #if metric == "hbonds":
    #    fig, ax = plt.subplots(figsize=(30,15))
    #    plt.subplot( 111 )
    #    plt.scatter( x_data, y_data )
    #    plt.plot( x_data, [ ( m*x + b ) for x in x_data ], c="red" )
    #    plot_title = "hbonds"
    #    plt.suptitle( plot_title, fontsize = 24 )
    #    plt.subplots_adjust(top=0.87)
    #    plt.savefig( plot_title, dpi=120, transparent=True )

    return r


def print_r_squared_data( data, r_squared_dict, name ):
    r_squared_keys = r_squared_dict.keys()
    r_squared_keys.sort( reverse = True )
    
    for r_squared in r_squared_keys:
        if r_squared > 0.3:
            print r_squared, r_squared_dict[ r_squared ], name


def print_other_data( data ):
    print "MC min:", min( data[ "MonteCarlo_acceptance_rate" ] ), "max:", max( data[ "MonteCarlo_acceptance_rate" ] ), "mean:", round( np.mean( data[ "MonteCarlo_acceptance_rate" ] ), 1 ), "median:", np.median( data[ "MonteCarlo_acceptance_rate" ] )

    top10_total_score_data = data.sort( "total_score" ).head( 10 )
    rmsd_count = len( top10_total_score_data[ top10_total_score_data[ "Fc_glycan_rmsd" ] <= 1.5 ] )
    print "Top10 total_score count Fc_glycan_rmsd <= 1.5:", rmsd_count
    Fnat_count = len( top10_total_score_data[ top10_total_score_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"] >= 85.0 ] )
    print "Top10 total_score count Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A >= 85%:", Fnat_count

    top10_pseudo_interface_energy_data = data.sort( "pseudo_interface_energy" ).head( 10 )
    rmsd_count = len( top10_pseudo_interface_energy_data[ top10_pseudo_interface_energy_data[ "Fc_glycan_rmsd" ] <= 1.5 ] )
    print "Top10 pseudo_interface_energy count Fc_glycan_rmsd <= 1.5:", rmsd_count
    Fnat_count = len( top10_pseudo_interface_energy_data[ top10_pseudo_interface_energy_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"] >= 85.0 ] )
    print "Top10 pseudo_interface_energy count Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A >= 85%:", Fnat_count



low_E_native_data = pd.read_csv( "/Users/Research/pyrosetta_dir/metric_data/fresh_start/get_base_pack_min_pose.csv" )
decoys_below_2A_dict = {}
decoys_above_F_protein_80_Fnat = {}



#############################################################
#### SSM-50 data using full Fc glycan LCM reset and ramp ####
#############################################################
###########
#### 1 ####
###########
# am2_5_mpt
path_to_using_native_full_glycan_50_SSM_am2_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_am2_5_mpt.csv"
using_native_full_glycan_50_SSM_am2_5_mpt_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_am2_5_mpt )
path_to_using_native_full_glycan_50_SSM_am2_5_mpt_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_am2_5_mpt_no_reset.csv"
using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_am2_5_mpt_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_am2_5_mpt_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ] - using_native_full_glycan_50_SSM_am2_5_mpt_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_am2_5_mpt_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_am2_5_mpt_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_am2_5_mpt_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_am2_5_mpt_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_am2_5_mpt_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_5_mpt_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ], c=using_native_full_glycan_50_SSM_am2_5_mpt_data[ "sugar_bb" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_am2_5_mpt_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_am2_5_mpt_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, ramp, am2, 5 mpt - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_am2_5_mpt_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_am2_5_mpt_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_am2_5_mpt_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_am2_5_mpt_data )
print "50 25 no cst"
print "\n\n"




###########
#### 2 ####
###########
# am2_5_mpt_glycan_to_protein_atoms_tightest, only one atom per glycan is constrained
path_to_using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt.csv"
using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt )
path_to_using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset.csv"
using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ] - using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data, metric )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ], c=using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "atom_pair_constraint" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, glycan to protein atoms tightest cst, ramp, am2, 5 mpt - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_glycan_to_protein_atoms_tightest_am2_5_mpt_data )
print "50 25 glycan_to_protein_atoms_tightest"
print "\n\n"




###########
#### 3 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt )
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric

fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ], c=using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "atom_pair_constraint" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, two glycan to protein atoms tightest cst, ramp, am2, 5 mpt - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_data )
print "50 25 two_glycan_to_protein_atoms_tightest"
print "\n\n"




###########
#### 4 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest_half_fa_sol, with ramp, half_fa_sol.sf. Comparing with/without LCM reset
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol )
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data, metric )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ], c=using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "atom_pair_constraint" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, two glycan to protein atoms tightest cst, ramp, am2, 5 mpt, half fa_sol - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_half_fa_sol_data )
print "50 25 two_glycan_to_protein_atoms_tightest, half fa_sol"
print "\n\n"




###########
#### 5 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest_double_hbond, with ramp, double_hbond.sf. Comparing with/without LCM reset
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond )
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset )
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data = using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data.sort("total_score")

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data, metric )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ], c=using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "atom_pair_constraint" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, two glycan to protein atoms tightest cst, ramp, am2, 5 mpt, double hbond - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_hbond_data )
print "50 25 two_glycan_to_protein_atoms_tightest, double hbond"
print "\n\n"




###########
#### 6 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest_triple_hbond, with ramp, triple_hbond.sf. Comparing with/without LCM reset
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond )
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data, metric )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ], c=using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, two glycan to protein atoms tightest cst, ramp, am2, 5 mpt, triple hbond - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_triple_hbond_data )
print "50 25 two_glycan_to_protein_atoms_tightest, triple hbond"
print "\n\n"



###########
#### 7 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest_double_fa_elec, with ramp, double_fa_elec.sf. Comparing with/without LCM reset
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec )
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data, metric )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ], c=using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "atom_pair_constraint" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, two glycan to protein atoms tightest cst, ramp, am2, 5 mpt, double fa_elec - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_elec_data )
print "50 25 two_glycan_to_protein_atoms_tightest, double fa_elec"
print "\n\n"




###########
#### 8 ####
###########
# am2_5_mpt_two_glycan_to_protein_atoms_tightest_double_fa_atr_half_fa_rep, with ramp, double_fa_atr_half_fa_rep.sf. Comparing with/without LCM reset
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep )
path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset = "/Users/Research/pyrosetta_dir/metric_data/fresh_start/using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset.csv"
using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ] - using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "delta_total_score" ] )

metrics = list( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data, metric )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 321 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 322 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ], c=using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ] )
#plt.colorbar(sc)
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 323 )
x = 100.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "pseudo_interface_energy" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ min(ymins) - 1, ymax + 1 ] )

plt.subplot( 324 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 325 )
xmin = floor( min( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 326 )
plt.hist( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" ], histtype="stepfilled" )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A" )
plt.xlim( [ -1, 100 ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "3ay4 using SSM-50 on Fc glycan with LCM reset, two glycan to protein atoms tightest cst, ramp, am2, 5 mpt, double fa_atr and half fa_rep - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_two_glycan_to_protein_atoms_tightest_am2_5_mpt_double_fa_atr_half_fa_rep_data )
print "50 25 two_glycan_to_protein_atoms_tightest, double fa_atr and half fa_rep"
print "\n\n"
