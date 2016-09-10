#!/usr/bin/python

import matplotlib.pyplot as plt
plt.rcParams.update( { "font.size" : 24 } )
import numpy as np
from scipy.stats import linregress, normaltest
import pandas as pd
from collections import Counter
from math import floor, ceil, log10, sqrt
import sys



def stdev( lst ):
    """
    returns the standard deviation of <lst>
    """
    mn = sum( lst ) / len( lst )
    variance = sum( [ ( e - mn ) ** 2 for e in lst ] )
    return sqrt( variance )


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

    top10_total_score_data = data.sort_values( "total_score" ).head( 10 )
    rmsd_count = len( top10_total_score_data[ top10_total_score_data[ "Fc_glycan_rmsd" ] <= 1.5 ] )
    print "Top10 total_score count Fc_glycan_rmsd <= 1.5:", rmsd_count
    #Fnat_count = len( top10_total_score_data[ top10_total_score_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"] >= 85.0 ] )
    #print "Top10 total_score count Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A >= 85%:", Fnat_count

    #top10_pseudo_interface_energy_data = data.sort_values( "pseudo_interface_energy" ).head( 10 )
    #rmsd_count = len( top10_pseudo_interface_energy_data[ top10_pseudo_interface_energy_data[ "Fc_glycan_rmsd" ] <= 1.5 ] )
    #print "Top10 pseudo_interface_energy count Fc_glycan_rmsd <= 1.5:", rmsd_count
    #Fnat_count = len( top10_pseudo_interface_energy_data[ top10_pseudo_interface_energy_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A"] >= 85.0 ] )
    #print "Top10 pseudo_interface_energy count Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A >= 85%:", Fnat_count





##############################################################
#### SSM-50 data using stepwise glycan LCM reset and ramp ####
##############################################################
###########
#### 1 ####
###########
# stepwise_by_removal/glycan_sampling_just_50_sugar_small_moves_with_ramp/am2_2_mpt_up_to_res_2
path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native = "/Users/mlnance/pyrosetta_dir/metric_data/stepwise_by_removal/using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native.csv"
using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native )
#path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/stepwise_by_removal/using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_no_reset.csv"
#using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_no_reset )

#low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ] - using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "delta_total_score" ] )

low_E_APhi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "APhi" ]
APhi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "APhi" ] )
low_E_APsi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "APsi" ]
APsi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "APsi" ] )
print "APhi:", low_E_APhi, "APhi stdev", APhi_stdev, "up_to_res_2"
print "APsi:", low_E_APsi, "APsi stdev", APsi_stdev, "up_to_res_2"
low_E_BPhi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "BPhi" ]
BPhi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "BPhi" ] )
low_E_BPsi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "BPsi" ]
BPsi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "BPsi" ] )
print "BPhi:", low_E_BPhi, "BPhi stdev", BPhi_stdev, "up_to_res_2"
print "BPsi:", low_E_BPsi, "BPsi stdev", BPsi_stdev, "up_to_res_2"


metrics = list( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ], c=using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
xmin = floor( min( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 223 )
xmin = -180
xmax = 180
ymin = -180
ymax = 180
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "APhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "APsi" ] )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "APhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "APsi" ], c=using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] )
plt.colorbar(sc)
plt.xlabel( "APhi" )
plt.xlim( [ xmin, xmax ] )
plt.xticks( np.arange( xmin, xmax + 1, 90.0 ) )
plt.ylabel( "APsi" )
plt.ylim( [ ymin, ymax ] )
plt.yticks( np.arange( ymin, ymax + 1, 90.0 ) )

plt.subplot( 224 )
xmin = -180
xmax = 180
ymin = -180
ymax = 180
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "BPhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "BPsi" ] )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "BPhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "BPsi" ], c=using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] )
plt.colorbar(sc)
plt.xlabel( "BPhi" )
plt.xlim( [ xmin, xmax ] )
plt.xticks( np.arange( xmin, xmax + 1, 90.0 ) )
plt.ylabel( "BPsi" )
plt.ylim( [ ymin, ymax ] )
plt.yticks( np.arange( ymin, ymax + 1, 90.0 ) )


# save the plot
plt.tight_layout()
plot_title = "fa_intra_rep 3ay4 using SSM-50 on Fc glycan with LCM reset, stepwise by removal, up_to_res_2, ramp, am2, 2 mpt - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_2_using_fa_intra_rep_native_data )
print "50 22 stepwise by removal, up_to_res_2"
print "\n\n"




###########
#### 2 ####
###########
# stepwise_by_removal/glycan_sampling_just_50_sugar_small_moves_with_ramp/am2_2_mpt_up_to_res_3
path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native = "/Users/mlnance/pyrosetta_dir/metric_data/stepwise_by_removal/using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native.csv"
using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native )
#path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/stepwise_by_removal/using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_no_reset.csv"
#using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_no_reset )

#low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "pseudo_interface_energy" ] - using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ] - using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "delta_total_score" ] )

low_E_APhi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "APhi" ]
APhi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "APhi" ] )
low_E_APsi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "APsi" ]
APsi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "APsi" ] )
print "APhi:", low_E_APhi, "APhi stdev", APhi_stdev, "up_to_res_3"
print "APsi:", low_E_APsi, "APsi stdev", APsi_stdev, "up_to_res_3"
low_E_BPhi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "BPhi" ]
BPhi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "BPhi" ] )
low_E_BPsi = using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data.loc[ using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ].argmin() ][ "BPsi" ]
BPsi_stdev = np.std( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "BPsi" ] )
print "BPhi:", low_E_BPhi, "BPhi stdev", BPhi_stdev, "up_to_res_3"
print "BPsi:", low_E_BPsi, "BPsi stdev", BPsi_stdev, "up_to_res_3"


metrics = list( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        r = get_r_of_line_of_best_fit( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data, metric, log10 )
        log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ], c=using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
xmin = floor( min( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 223 )
xmin = -180
xmax = 180
ymin = -180
ymax = 180
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "APhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "APsi" ] )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "APhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "APsi" ], c=using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] )
plt.colorbar(sc)
plt.xlabel( "APhi" )
plt.xlim( [ xmin, xmax ] )
plt.xticks( np.arange( xmin, xmax + 1, 90.0 ) )
plt.ylabel( "APsi" )
plt.ylim( [ ymin, ymax ] )
plt.yticks( np.arange( ymin, ymax + 1, 90.0 ) )

plt.subplot( 224 )
xmin = -180
xmax = 180
ymin = -180
ymax = 180
#sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "BPhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "BPsi" ] )
sc = plt.scatter( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "BPhi" ], using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "BPsi" ], c=using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data[ "Fc_glycan_rmsd" ] )
plt.colorbar(sc)
plt.xlabel( "BPhi" )
plt.xlim( [ xmin, xmax ] )
plt.xticks( np.arange( xmin, xmax + 1, 90.0 ) )
plt.ylabel( "BPsi" )
plt.ylim( [ ymin, ymax ] )
plt.yticks( np.arange( ymin, ymax + 1, 90.0 ) )


# save the plot
plt.tight_layout()
plot_title = "fa_intra_rep 3ay4 using SSM-50 on Fc glycan with LCM reset, stepwise by removal, up_to_res_3, ramp, am2, 2 mpt - compared against protocol without reset"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
print_r_squared_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data, r_squared_to_metric_dict, "linear" )
print_r_squared_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data, log10_r_squared_to_metric_dict, "log10" )
print_r_squared_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data, binned_r_squared_to_metric_dict, "binned rmsd" )
print_other_data( using_native_full_glycan_50_SSM_am2_2_mpt_up_to_res_3_using_fa_intra_rep_native_data )
print "50 22 stepwise by removal, up_to_res_3"
print "\n\n"
