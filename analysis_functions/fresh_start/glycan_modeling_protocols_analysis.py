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
    do_plot = False
    if metric == "total_score":
        plot_title = "total_score"
        #do_plot = True
    if do_plot:
        fig, ax = plt.subplots(figsize=(30,15))
        plt.subplot( 111 )
        plt.scatter( x_data, y_data )
        plt.plot( x_data, [ ( m*x + b ) for x in x_data ], c="red" )
        plt.suptitle( plot_title, fontsize = 24 )
        plt.subplots_adjust(top=0.87)
        plt.savefig( plot_title, dpi=120, transparent=True )

    return r



def print_r_squared_data( data, r_squared_dict, name ):
    r_squared_keys = r_squared_dict.keys()
    r_squared_keys.sort( reverse = True )
    
    for r_squared in r_squared_keys:
        if r_squared > 0.5:
            print r_squared, r_squared_dict[ r_squared ], name


def get_top10_data( data ):
    print "MC min:", min( data[ "MonteCarlo_acceptance_rate" ] ), "max:", max( data[ "MonteCarlo_acceptance_rate" ] ), "mean:", round( np.mean( data[ "MonteCarlo_acceptance_rate" ] ), 1 ), "median:", np.median( data[ "MonteCarlo_acceptance_rate" ] )

    top10_total_score_data = data.sort_values( "total_score" ).head( 10 )
    top10_Fc_glycan_rmsd_data = data.sort_values( "Fc_glycan_rmsd" ).head( 10 )
    total_score_rmsd_count_2 = len( top10_total_score_data[ top10_total_score_data[ "Fc_glycan_rmsd" ] <= 2.0 ] )
    print "Top10 total_score count Fc_glycan_rmsd <= 2.0:", total_score_rmsd_count_2
    total_score_rmsd_count_1_point_5 = len( top10_total_score_data[ top10_total_score_data[ "Fc_glycan_rmsd" ] <= 1.5 ] )
    print "Top10 total_score count Fc_glycan_rmsd <= 1.5:", total_score_rmsd_count_1_point_5
    total_score_rmsd_count_1 = len( top10_total_score_data[ top10_total_score_data[ "Fc_glycan_rmsd" ] <= 1.0 ] )
    print "Top10 total_score count Fc_glycan_rmsd <= 1.0:", total_score_rmsd_count_1

    top10_total_score_filenames = top10_total_score_data[ "filename" ]
    top10_Fc_glycan_rmsd_filenames = top10_Fc_glycan_rmsd_data[ "filename" ]
    total_score_num_hits = len( list( set( top10_total_score_filenames ) & set( top10_Fc_glycan_rmsd_filenames ) ) )
    print "The top10 total_score filter got %s of the top 10 decoys with the lowest Fc_glcyan_rmsd" %total_score_num_hits
    
    try:
        # sorted on Fnat, then the top10 of the highest Fnat with the lowest rmsd are kept
        top_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_data = data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )
        top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_data = top_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_data.sort_values( "glycan_rmsd", ascending=True ).head(10)
        #Fnat_count = len( top10_total_score_data[ top10_total_score_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A"] >= 85.0 ] )
        #print "Top10 total_score count Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A >= 85%:", Fnat_count

        #top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_filenames = top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_data[ "filename" ]
        #total_score_num_hits_Fnat = len( list( set( top10_total_score_filenames ) & set( top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_filenames ) ) )
        #print "The top10 total_score filter got %s of the top 10 decoys with the highest Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" %total_score_num_hits_Fnat
    except:
        pass

    try:
        top10_pseudo_interface_energy_data = data.sort_values( "pseudo_interface_energy" ).head( 10 )
        top10_Fc_glycan_rmsd_data = data.sort_values( "Fc_glycan_rmsd" ).head( 10 )
        rmsd_count = len( top10_pseudo_interface_energy_data[ top10_pseudo_interface_energy_data[ "Fc_glycan_rmsd" ] <= 2.0 ] )
        print "Top10 pseudo_interface_energy count Fc_glycan_rmsd <= 2.0:", rmsd_count
        rmsd_count = len( top10_pseudo_interface_energy_data[ top10_pseudo_interface_energy_data[ "Fc_glycan_rmsd" ] <= 1.5 ] )
        print "Top10 pseudo_interface_energy count Fc_glycan_rmsd <= 1.5:", rmsd_count
        rmsd_count = len( top10_pseudo_interface_energy_data[ top10_pseudo_interface_energy_data[ "Fc_glycan_rmsd" ] <= 1.0 ] )
        print "Top10 pseudo_interface_energy count Fc_glycan_rmsd <= 1.0:", rmsd_count

        top10_pseudo_interface_energy_filenames = top10_pseudo_interface_energy_data[ "filename" ]
        top10_Fc_glycan_rmsd_filenames = top10_Fc_glycan_rmsd_data[ "filename" ]
        pseudo_interface_energy_num_hits = len( list( set( top10_pseudo_interface_energy_filenames ) & set( top10_Fc_glycan_rmsd_filenames ) ) )
        print "The top10 pseudo_interface_energy filter got %s of the top 10 decoys with the lowest Fc_glcyan_rmsd" %pseudo_interface_energy_num_hits
    except:
        pass

    return total_score_num_hits, pseudo_interface_energy_num_hits, total_score_rmsd_count_2, total_score_rmsd_count_1_point_5, total_score_rmsd_count_1



def get_top10_Fnat_data( data ):
    top10_total_score_data = data.sort_values( "total_score" ).head( 10 )
    top10_total_score_filenames = top10_total_score_data[ "filename" ]

    top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_data = data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False ).head( 10 )
    top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_filenames = top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_data[ "filename" ]

    Fnat_count = len( top10_total_score_data[ top10_total_score_data["Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A"] >= 85.0 ] )
    print "Top10 total_score count Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A >= 85%:", Fnat_count
    
    top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_filenames = top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_data[ "filename" ]
    total_score_num_hits_Fnat = len( list( set( top10_total_score_filenames ) & set( top10_Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A_filenames ) ) )
    print "The top10 total_score filter got %s of the top 10 decoys with the highest Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" %total_score_num_hits_Fnat





##############################################
#### SSM-N data using different protocols ####
##############################################
'''
###########
#### 0 ####
###########
# glycan_modeling_protocols/protocol_0
path_to_using_native_full_glycan_protocol_0 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_0.csv"
using_native_full_glycan_protocol_0_data = pd.read_csv( path_to_using_native_full_glycan_protocol_0 )
#path_to_using_native_full_glycan_protocol_0_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_0_no_reset.csv"
#using_native_full_glycan_protocol_0_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_0_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_0_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_0_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_0_data[ "total_score" ] - using_native_full_glycan_protocol_0_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_0"
print using_native_full_glycan_protocol_0_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_0_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_0"
print using_native_full_glycan_protocol_0_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_0"
print using_native_full_glycan_protocol_0_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_0_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_0_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_0_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_0_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_0_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_0_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_0_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_0_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_0_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_0_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_0_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_0_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_0_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_0_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_0_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_0_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_0_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_0_data[ "total_score" ], c=using_native_full_glycan_protocol_0_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_0_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_0_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_0_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_0_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_0_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_0"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_0_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_0_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_0_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_0_total_score_num_hits, protocol_0_pseudo_interface_energy_num_hits, protocol_0_total_score_glycan_rmsd_count_2, protocol_0_total_score_glycan_rmsd_count_1_point_5, protocol_0_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_0_data )
protocol_0_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_0_data )
print "protocol_0\n\n\n"



###########
#### 1 ####
###########
# glycan_modeling_protocols/protocol_1
path_to_using_native_full_glycan_protocol_1 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_1.csv"
using_native_full_glycan_protocol_1_data = pd.read_csv( path_to_using_native_full_glycan_protocol_1 )
#path_to_using_native_full_glycan_protocol_1_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_1_no_reset.csv"
#using_native_full_glycan_protocol_1_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_1_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_1_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_1_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_1_data[ "total_score" ] - using_native_full_glycan_protocol_1_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_1"
print using_native_full_glycan_protocol_1_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_1_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_1"
print using_native_full_glycan_protocol_1_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_1"
print using_native_full_glycan_protocol_1_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_1_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_1_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_1_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_1_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_1_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_1_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_1_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_1_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_1_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_1_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_1_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_1_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_1_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_1_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_1_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_1_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_1_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_1_data[ "total_score" ], c=using_native_full_glycan_protocol_1_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_1_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_1_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_1_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_1_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_1_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_1"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_1_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_1_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_1_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_1_total_score_num_hits, protocol_1_pseudo_interface_energy_num_hits, protocol_1_total_score_glycan_rmsd_count_2, protocol_1_total_score_glycan_rmsd_count_1_point_5, protocol_1_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_1_data )
protocol_1_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_1_data )
print "protocol_1\n\n\n"



###########
#### 2 ####
###########
# glycan_modeling_protocols/protocol_2
path_to_using_native_full_glycan_protocol_2 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_2.csv"
using_native_full_glycan_protocol_2_data = pd.read_csv( path_to_using_native_full_glycan_protocol_2 )
#path_to_using_native_full_glycan_protocol_2_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_2_no_reset.csv"
#using_native_full_glycan_protocol_2_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_2_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_2_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_2_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_2_data[ "total_score" ] - using_native_full_glycan_protocol_2_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_2"
print using_native_full_glycan_protocol_2_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_2_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_2"
print using_native_full_glycan_protocol_2_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_2"
print using_native_full_glycan_protocol_2_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_2_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_2_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_2_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_2_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_2_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_2_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_2_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_2_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_2_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_2_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_2_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_2_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_2_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_2_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_2_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_2_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_2_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_2_data[ "total_score" ], c=using_native_full_glycan_protocol_2_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_2_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_2_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_2_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_2_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_2_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_2"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_2_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_2_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_2_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_2_total_score_num_hits, protocol_2_pseudo_interface_energy_num_hits, protocol_2_total_score_glycan_rmsd_count_2, protocol_2_total_score_glycan_rmsd_count_1_point_5, protocol_2_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_2_data )
protocol_2_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_2_data )
print "protocol_2\n\n\n"



###########
#### 3 ####
###########
# glycan_modeling_protocols/protocol_3
path_to_using_native_full_glycan_protocol_3 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_3.csv"
using_native_full_glycan_protocol_3_data = pd.read_csv( path_to_using_native_full_glycan_protocol_3 )
#path_to_using_native_full_glycan_protocol_3_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_3_no_reset.csv"
#using_native_full_glycan_protocol_3_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_3_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_3_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_3_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_3_data[ "total_score" ] - using_native_full_glycan_protocol_3_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_3"
print using_native_full_glycan_protocol_3_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_3_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_3"
print using_native_full_glycan_protocol_3_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_3"
print using_native_full_glycan_protocol_3_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_3_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_3_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_3_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_3_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_3_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_3_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_3_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_3_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_3_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_3_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_3_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_3_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_3_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_3_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_3_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_3_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_3_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_3_data[ "total_score" ], c=using_native_full_glycan_protocol_3_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_3_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_3_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_3_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_3_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_3_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_3"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_3_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_3_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_3_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_3_total_score_num_hits, protocol_3_pseudo_interface_energy_num_hits, protocol_3_total_score_glycan_rmsd_count_2, protocol_3_total_score_glycan_rmsd_count_1_point_5, protocol_3_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_3_data )
protocol_3_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_3_data )
print "protocol_3\n\n\n"



###########
#### 4 ####
###########
# glycan_modeling_protocols/protocol_4
path_to_using_native_full_glycan_protocol_4 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_4.csv"
using_native_full_glycan_protocol_4_data = pd.read_csv( path_to_using_native_full_glycan_protocol_4 )
#path_to_using_native_full_glycan_protocol_4_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_4_no_reset.csv"
#using_native_full_glycan_protocol_4_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_4_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_4_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_4_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_4_data[ "total_score" ] - using_native_full_glycan_protocol_4_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_4"
print using_native_full_glycan_protocol_4_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_4_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_4"
print using_native_full_glycan_protocol_4_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_4"
print using_native_full_glycan_protocol_4_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_4_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_4_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_4_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_4_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_4_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_4_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_4_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_4_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_4_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_4_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_4_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_4_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_4_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_4_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_4_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_4_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_4_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_4_data[ "total_score" ], c=using_native_full_glycan_protocol_4_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_4_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_4_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_4_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_4_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_4_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_4"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_4_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_4_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_4_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_4_total_score_num_hits, protocol_4_pseudo_interface_energy_num_hits, protocol_4_total_score_glycan_rmsd_count_2, protocol_4_total_score_glycan_rmsd_count_1_point_5, protocol_4_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_4_data )
protocol_4_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_4_data )
print "protocol_4\n\n\n"



###########
#### 5 ####
###########
# glycan_modeling_protocols/protocol_5
path_to_using_native_full_glycan_protocol_5 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_5.csv"
using_native_full_glycan_protocol_5_data = pd.read_csv( path_to_using_native_full_glycan_protocol_5 )
#path_to_using_native_full_glycan_protocol_5_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_5_no_reset.csv"
#using_native_full_glycan_protocol_5_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_5_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_5_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_5_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_5_data[ "total_score" ] - using_native_full_glycan_protocol_5_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_5"
print using_native_full_glycan_protocol_5_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_5_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_5"
print using_native_full_glycan_protocol_5_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_5"
print using_native_full_glycan_protocol_5_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_5_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_5_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_5_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_5_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_5_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_5_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_5_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_5_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_5_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_5_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_5_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_5_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_5_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_5_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_5_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_5_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_5_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_5_data[ "total_score" ], c=using_native_full_glycan_protocol_5_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_5_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_5_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_5_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_5_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_5_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_5"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_5_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_5_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_5_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_5_total_score_num_hits, protocol_5_pseudo_interface_energy_num_hits, protocol_5_total_score_glycan_rmsd_count_2, protocol_5_total_score_glycan_rmsd_count_1_point_5, protocol_5_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_5_data )
protocol_5_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_5_data )
print "protocol_5\n\n\n"



###########
#### 6 ####
###########
# glycan_modeling_protocols/protocol_6
path_to_using_native_full_glycan_protocol_6 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_6.csv"
using_native_full_glycan_protocol_6_data = pd.read_csv( path_to_using_native_full_glycan_protocol_6 )
#path_to_using_native_full_glycan_protocol_6_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_6_no_reset.csv"
#using_native_full_glycan_protocol_6_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_6_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_6_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_6_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_6_data[ "total_score" ] - using_native_full_glycan_protocol_6_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_6"
print using_native_full_glycan_protocol_6_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_6_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_6"
print using_native_full_glycan_protocol_6_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_6"
print using_native_full_glycan_protocol_6_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

#print "\n\n\n\n", using_native_full_glycan_protocol_6_data[ ( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ] >= 7 ) & ( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ] <= 8 ) & ( using_native_full_glycan_protocol_6_data[ "total_score" ] >= -535 ) & ( using_native_full_glycan_protocol_6_data[ "total_score" ] <= -520 ) ][ ["filename", "total_score", "Fc_glycan_rmsd" ] ].sort_values( "total_score" )
#sys.exit()

metrics = list( using_native_full_glycan_protocol_6_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_6_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_6_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_6_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_6_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_6_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_6_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_6_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_6_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_6_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_6_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_6_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_6_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_6_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_6_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_6_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_6_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_6_data[ "total_score" ], c=using_native_full_glycan_protocol_6_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_6_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_6_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_6_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_6_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_6_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_6"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_6_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_6_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_6_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_6_total_score_num_hits, protocol_6_pseudo_interface_energy_num_hits, protocol_6_total_score_glycan_rmsd_count_2, protocol_6_total_score_glycan_rmsd_count_1_point_5, protocol_6_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_6_data )
protocol_6_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_6_data )
print "protocol_6\n\n\n"



###########
#### 7 ####
###########
# glycan_modeling_protocols/protocol_7
path_to_using_native_full_glycan_protocol_7 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_7.csv"
using_native_full_glycan_protocol_7_data = pd.read_csv( path_to_using_native_full_glycan_protocol_7 )
#path_to_using_native_full_glycan_protocol_7_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_7_no_reset.csv"
#using_native_full_glycan_protocol_7_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_7_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_7_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_7_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_7_data[ "total_score" ] - using_native_full_glycan_protocol_7_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_7"
print using_native_full_glycan_protocol_7_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_7_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_7"
print using_native_full_glycan_protocol_7_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_7"
print using_native_full_glycan_protocol_7_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_7_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_7_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_7_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_7_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_7_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_7_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_7_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_7_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_7_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_7_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_7_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_7_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_7_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_7_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_7_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_7_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_7_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_7_data[ "total_score" ], c=using_native_full_glycan_protocol_7_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_7_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_7_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_7_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_7_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_7_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_7"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_7_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_7_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_7_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_7_total_score_num_hits, protocol_7_pseudo_interface_energy_num_hits, protocol_7_total_score_glycan_rmsd_count_2, protocol_7_total_score_glycan_rmsd_count_1_point_5, protocol_7_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_7_data )
protocol_7_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_7_data )
print "protocol_7\n\n\n"



###########
#### 8 ####
###########
# glycan_modeling_protocols/protocol_8
path_to_using_native_full_glycan_protocol_8 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_8.csv"
using_native_full_glycan_protocol_8_data = pd.read_csv( path_to_using_native_full_glycan_protocol_8 )
#path_to_using_native_full_glycan_protocol_8_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_8_no_reset.csv"
#using_native_full_glycan_protocol_8_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_8_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_8_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_8_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_8_data[ "total_score" ] - using_native_full_glycan_protocol_8_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_8"
print using_native_full_glycan_protocol_8_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_8_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_8"
print using_native_full_glycan_protocol_8_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_8"
print using_native_full_glycan_protocol_8_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_8_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_8_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_8_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_8_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_8_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_8_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_8_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_8_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_8_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_8_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_8_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_8_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_8_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_8_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_8_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_8_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_8_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_8_data[ "total_score" ], c=using_native_full_glycan_protocol_8_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_8_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_8_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_8_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_8_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_8_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_8"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_8_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_8_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_8_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_8_total_score_num_hits, protocol_8_pseudo_interface_energy_num_hits, protocol_8_total_score_glycan_rmsd_count_2, protocol_8_total_score_glycan_rmsd_count_1_point_5, protocol_8_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_8_data )
protocol_8_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_8_data )
print "protocol_8\n\n\n"



###########
#### 9 ####
###########
# glycan_modeling_protocols/protocol_9
path_to_using_native_full_glycan_protocol_9 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_9.csv"
using_native_full_glycan_protocol_9_data = pd.read_csv( path_to_using_native_full_glycan_protocol_9 )
#path_to_using_native_full_glycan_protocol_9_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_9_no_reset.csv"
#using_native_full_glycan_protocol_9_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_9_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_9_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_9_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_9_data[ "total_score" ] - using_native_full_glycan_protocol_9_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_9"
print using_native_full_glycan_protocol_9_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_9_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_9"
print using_native_full_glycan_protocol_9_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_9"
print using_native_full_glycan_protocol_9_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
for x in range( 1, 11 ):
    print "\n%s : %s\n" %( x-1, x ), using_native_full_glycan_protocol_9_data[ ( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ] >= x-1 ) & ( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ] <= x ) ].sort_values( "total_score", ascending=True )[ ["filename", "total_score", "Fc_glycan_rmsd"] ][ :3 ]


metrics = list( using_native_full_glycan_protocol_9_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_9_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_9_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_9_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_9_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_9_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_9_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_9_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_9_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_9_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_data[ "total_score" ], c=using_native_full_glycan_protocol_9_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_9_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_9_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_9_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_9"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_9_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_9_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_9_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_9_total_score_num_hits, protocol_9_pseudo_interface_energy_num_hits, protocol_9_total_score_glycan_rmsd_count_2, protocol_9_total_score_glycan_rmsd_count_1_point_5, protocol_9_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_9_data )
protocol_9_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_9_data )
print "protocol_9\n\n\n"



############
#### 10 ####
############
# glycan_modeling_protocols/protocol_10
path_to_using_native_full_glycan_protocol_10 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_10.csv"
using_native_full_glycan_protocol_10_data = pd.read_csv( path_to_using_native_full_glycan_protocol_10 )
#path_to_using_native_full_glycan_protocol_10_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_10_no_reset.csv"
#using_native_full_glycan_protocol_10_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_10_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_10_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_10_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_10_data[ "total_score" ] - using_native_full_glycan_protocol_10_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_10"
print using_native_full_glycan_protocol_10_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_10_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_10"
print using_native_full_glycan_protocol_10_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_10"
print using_native_full_glycan_protocol_10_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_10_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_10_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_10_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_10_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_10_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_10_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_10_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_10_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_10_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_10_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_data[ "total_score" ], c=using_native_full_glycan_protocol_10_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_10_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_10_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_10_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_10"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_10_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_10_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_10_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_10_total_score_num_hits, protocol_10_pseudo_interface_energy_num_hits, protocol_10_total_score_glycan_rmsd_count_2, protocol_10_total_score_glycan_rmsd_count_1_point_5, protocol_10_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_10_data )
protocol_10_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_10_data )
print "protocol_10\n\n\n"
'''


############
#### 11 ####
############
# glycan_modeling_protocols/protocol_11
path_to_using_native_full_glycan_protocol_11 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_11.csv"
using_native_full_glycan_protocol_11_data = pd.read_csv( path_to_using_native_full_glycan_protocol_11 )
#path_to_using_native_full_glycan_protocol_11_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_11_no_reset.csv"
#using_native_full_glycan_protocol_11_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_11_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_11_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_11_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_11_data[ "total_score" ] - using_native_full_glycan_protocol_11_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_11"
print using_native_full_glycan_protocol_11_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_11_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_11"
print using_native_full_glycan_protocol_11_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_11"
print using_native_full_glycan_protocol_11_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print "\n", using_native_full_glycan_protocol_11_data[ ( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ] >= 4 ) & ( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ] <= 6 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :5 ]

metrics = list( using_native_full_glycan_protocol_11_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_11_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_11_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_11_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_11_data, metric)
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_11_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_11_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_11_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_11_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_11_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_11_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_11_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_11_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_11_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_11_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_11_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_11_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_11_data[ "total_score" ], c=using_native_full_glycan_protocol_11_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_11_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_11_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_11_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_11_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
#xmin = floor( min( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ] ) )
xmax = ceil( max( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ], histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_11_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_11"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_11_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_11_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_11_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_11_total_score_num_hits, protocol_11_pseudo_interface_energy_num_hits, protocol_11_total_score_glycan_rmsd_count_2, protocol_11_total_score_glycan_rmsd_count_1_point_5, protocol_11_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_11_data )
protocol_11_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_11_data )
print "protocol_11\n\n\n"




'''
###########################
#### COMPARE PROTOCOLS ####
###########################
# protocol_0_total_score_num_hits, protocol_0_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_0_data
# protocol_1_total_score_num_hits, protocol_1_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_1_data
# protocol_2_total_score_num_hits, protocol_2_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_2_data
# protocol_3_total_score_num_hits, protocol_3_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_3_data
# protocol_4_total_score_num_hits, protocol_4_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_4_data
# protocol_5_total_score_num_hits, protocol_5_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_5_data
# protocol_6_total_score_num_hits, protocol_6_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_6_data
# protocol_7_total_score_num_hits, protocol_7_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_7_data
# protocol_8_total_score_num_hits, protocol_8_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_8_data
# protocol_9_total_score_num_hits, protocol_9_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_9_data
# protocol_10_total_score_num_hits, protocol_10_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_10_data
# protocol_11_total_score_num_hits, protocol_11_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_11_data

names = ( "Protocol_0", "Protocol_1", "Protocol_2", "Protocol_3", "Protocol_4", "Protocol_5", "Protocol_6", "Protocol_7", "Protocol_8", "Protocol_9", "Protocol_10", "Protocol_11", )
total_score_protocol_data = ( protocol_0_total_score_num_hits, protocol_1_total_score_num_hits, protocol_2_total_score_num_hits, protocol_3_total_score_num_hits, protocol_4_total_score_num_hits, protocol_5_total_score_num_hits, protocol_6_total_score_num_hits, protocol_7_total_score_num_hits, protocol_8_total_score_num_hits, protocol_9_total_score_num_hits, protocol_10_total_score_num_hits, protocol_11_total_score_num_hits, )
pseudo_interface_energy_protocol_data = ( protocol_0_pseudo_interface_energy_num_hits, protocol_1_pseudo_interface_energy_num_hits, protocol_2_pseudo_interface_energy_num_hits, protocol_3_pseudo_interface_energy_num_hits, protocol_4_pseudo_interface_energy_num_hits, protocol_5_pseudo_interface_energy_num_hits, protocol_6_pseudo_interface_energy_num_hits, protocol_7_pseudo_interface_energy_num_hits, protocol_8_pseudo_interface_energy_num_hits, protocol_9_pseudo_interface_energy_num_hits, protocol_10_pseudo_interface_energy_num_hits, protocol_11_pseudo_interface_energy_num_hits, )
total_score_glycan_rmsd_count_2 = ( protocol_0_total_score_glycan_rmsd_count_2, protocol_1_total_score_glycan_rmsd_count_2, protocol_2_total_score_glycan_rmsd_count_2, protocol_3_total_score_glycan_rmsd_count_2, protocol_4_total_score_glycan_rmsd_count_2, protocol_5_total_score_glycan_rmsd_count_2, protocol_6_total_score_glycan_rmsd_count_2, protocol_7_total_score_glycan_rmsd_count_2, protocol_8_total_score_glycan_rmsd_count_2, protocol_9_total_score_glycan_rmsd_count_2, protocol_10_total_score_glycan_rmsd_count_2, protocol_11_total_score_glycan_rmsd_count_2, )
total_score_glycan_rmsd_count_1_point_5 = ( protocol_0_total_score_glycan_rmsd_count_1_point_5, protocol_1_total_score_glycan_rmsd_count_1_point_5, protocol_2_total_score_glycan_rmsd_count_1_point_5, protocol_3_total_score_glycan_rmsd_count_1_point_5, protocol_4_total_score_glycan_rmsd_count_1_point_5, protocol_5_total_score_glycan_rmsd_count_1_point_5, protocol_6_total_score_glycan_rmsd_count_1_point_5, protocol_7_total_score_glycan_rmsd_count_1_point_5, protocol_8_total_score_glycan_rmsd_count_1_point_5, protocol_9_total_score_glycan_rmsd_count_1_point_5, protocol_10_total_score_glycan_rmsd_count_1_point_5, protocol_11_total_score_glycan_rmsd_count_1_point_5, )
total_score_glycan_rmsd_count_1 = ( protocol_0_total_score_glycan_rmsd_count_1, protocol_1_total_score_glycan_rmsd_count_1, protocol_2_total_score_glycan_rmsd_count_1, protocol_3_total_score_glycan_rmsd_count_1, protocol_4_total_score_glycan_rmsd_count_1, protocol_5_total_score_glycan_rmsd_count_1, protocol_6_total_score_glycan_rmsd_count_1, protocol_7_total_score_glycan_rmsd_count_1, protocol_8_total_score_glycan_rmsd_count_1, protocol_9_total_score_glycan_rmsd_count_1, protocol_10_total_score_glycan_rmsd_count_1, protocol_11_total_score_glycan_rmsd_count_1, )

N = np.arange( len( names ) )
width = 0.3

# for top10 total_score that got top10 glycan_rmsd
fig, ax = plt.subplots(figsize=(28,15))
bins = map(lambda x: x - width/2,range( 1, len( names ) + 1 ) )
tick_bins = map(lambda x: x - width/2,range( 1, len( names ) + 1 ) )
bigger_bins = map(lambda x: x - width,range( 1, len( names ) + 1 ) )
total_score_data = ax.bar( bins, total_score_protocol_data, width, color='b' )
pseudo_interface_energy_data = ax.bar( bigger_bins, pseudo_interface_energy_protocol_data, width, color='r' )
ax.set_xticks( tick_bins )
ax.set_xticklabels( names, rotation_mode="anchor", ha="center" )
ax.set_xlabel( "Top10 total_score capturing Top10 Fc_glycan_rmsd" )
ax.set_ylabel( "Count" )
ax.set_ylim( [ 0, 10 ] )
ax.set_yticks( range( 0, 11 ) )
ax.legend( ( total_score_data[0], pseudo_interface_energy_data[0] ), ( "total_score", "pseudo_interface_energy" ) )

# save the plot
fig.tight_layout()
plot_title = "Comparing protocol variations of 3ay4-glycan modeling"
fig.suptitle( plot_title, fontsize = 20 )
fig.subplots_adjust(top=0.93)
fig.savefig( plot_title, dpi=120, transparent=True )
plt.close()


# for rmsd count below 2, 1.5, and 1
width = 0.25
fig, ax = plt.subplots(figsize=(26,14))
tick_bins = map(lambda x: x - width/2,range( 1, len( names ) + 1 ) )
bins = map(lambda x: x - width/4,range( 1, len( names ) + 1 ) )
bigger_bins = map(lambda x: x - width,range( 1, len( names ) + 1 ) )
evenbigger_bins = map(lambda x: x - width*2,range( 1, len( names ) + 1 ) )
total_score_glycan_rmsd_count_2_data = ax.bar( evenbigger_bins, total_score_glycan_rmsd_count_2, width, color='b', label="glycan_rmsd < 2" )
total_score_glycan_rmsd_count_1_point_5_data = ax.bar( bigger_bins, total_score_glycan_rmsd_count_1_point_5, width, color='r', label="glycan_rmsd < 1.5" )
total_score_glycan_rmsd_count_1_data = ax.bar( bins, total_score_glycan_rmsd_count_1, width, color='g', label="glycan_rmsd < 1" )
ax.set_xticks( tick_bins )
ax.set_xticklabels( names, rotation_mode="anchor", ha="center" )
ax.set_xlabel( "Protocol Number" )
ax.set_ylabel( "Count" )
ax.set_ylim( [ 0, 10 ] )
ax.set_yticks( range( 0, 11 ) )
leg = ax.legend(loc=2, bbox_to_anchor=(1,1))

# save the plot
fig.tight_layout()
plot_title = "Comparing protocol glycan rmsd below threshold variations of 3ay4-glycan modeling"
fig.suptitle( plot_title, fontsize = 20 )
fig.subplots_adjust(top=0.93)
#fig.savefig( plot_title, dpi=120, transparent=True )
fig.savefig( plot_title, bbox_extra_artists=(leg,), bbox_inches="tight", dpi=120, transparent=True )
plt.close()
'''
