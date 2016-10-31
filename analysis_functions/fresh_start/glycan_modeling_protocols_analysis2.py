#!/usr/bin/python

import matplotlib as mpl
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


def get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( data, metric, protocol_num ):
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
    if metric == "hbonds":
        plot_title = "hbonds_%s" %protocol_num
        do_plot = True
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
###########
#### 9 ####
###########
# glycan_modeling_protocols/protocol_9_again
path_to_using_native_full_glycan_protocol_9_again = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_9_again.csv"
using_native_full_glycan_protocol_9_again_data = pd.read_csv( path_to_using_native_full_glycan_protocol_9_again )
#path_to_using_native_full_glycan_protocol_9_again_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_9_again_no_reset.csv"
#using_native_full_glycan_protocol_9_again_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_9_again_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_9_again_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_9_again_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_9_again_data[ "total_score" ] - using_native_full_glycan_protocol_9_again_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_9_again"
print using_native_full_glycan_protocol_9_again_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_9_again_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_9_again"
print using_native_full_glycan_protocol_9_again_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_9_again"
print using_native_full_glycan_protocol_9_again_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
grab_these_structs = []
#for jj in range( 1, 11 ):
#    low_E_struct = using_native_full_glycan_protocol_9_again_data[ ( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ] >= jj - 1 ) & ( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ] <= jj ) ].sort_values( "total_score", ascending=True )
#    print jj-1, jj, "mean:", np.mean( low_E_struct["total_score"] ), "median:", np.median( low_E_struct["total_score"] ), "min:", min( low_E_struct["total_score"] ), "max:", max( low_E_struct["total_score"] ), "range:", max( low_E_struct["total_score"] ) - min( low_E_struct["total_score"] )
#    print "  /home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_9_again/base_structs/%s" %low_E_struct.iloc[0]["filename"]
#    grab_these_structs.append( "/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_9_again/base_structs/protocol_%s.gz\n" %low_E_struct.iloc[0]["filename"] )
#with open( "protocol_9_again_lowest_E_of_1A_binned_rmsd", "wb" ) as fh:
#    fh.writelines( grab_these_structs )

metrics = list( using_native_full_glycan_protocol_9_again_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_9_again_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_9_again_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_9_again_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_9_again_data, metric, "protocol_9_again" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_9_again_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_again_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_again_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_again_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_9_again_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_again_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_9_again_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_again_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_9_again_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_9_again_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_again_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_again_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_9_again_data[ "total_score" ], c=using_native_full_glycan_protocol_9_again_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_again_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_9_again_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_9_again_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_9_again_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_9_again_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

# save the plot
plt.tight_layout()
plot_title = "protocol_9_again"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_9_again_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_9_again_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_9_again_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_9_again_total_score_num_hits, protocol_9_again_pseudo_interface_energy_num_hits, protocol_9_again_total_score_glycan_rmsd_count_2, protocol_9_again_total_score_glycan_rmsd_count_1_point_5, protocol_9_again_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_9_again_data )
protocol_9_again_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_9_again_data )
print "protocol_9_again\n\n\n"


############
#### 10 ####
############
# glycan_modeling_protocols/protocol_10_again
path_to_using_native_full_glycan_protocol_10_again = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_10_again.csv"
using_native_full_glycan_protocol_10_again_data = pd.read_csv( path_to_using_native_full_glycan_protocol_10_again )
#path_to_using_native_full_glycan_protocol_10_again_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_10_again_no_reset.csv"
#using_native_full_glycan_protocol_10_again_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_10_again_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_10_again_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_10_again_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_10_again_data[ "total_score" ] - using_native_full_glycan_protocol_10_again_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_10_again"
print using_native_full_glycan_protocol_10_again_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_10_again_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_10_again"
print using_native_full_glycan_protocol_10_again_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_10_again"
print using_native_full_glycan_protocol_10_again_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_10_again_data[ ( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ] >= 4 ) & ( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ] <= 6 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :5 ]

metrics = list( using_native_full_glycan_protocol_10_again_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_10_again_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_10_again_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_10_again_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_10_again_data, metric, "protocol_10_again" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_10_again_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_again_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_again_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_again_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_10_again_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_again_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_10_again_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_again_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_10_again_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_10_again_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_again_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_again_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_10_again_data[ "total_score" ], c=using_native_full_glycan_protocol_10_again_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_again_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_10_again_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_again_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_10_again_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_10_again_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_10_again_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_10_again_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_10_again_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_10_again_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_10_again_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_10_again"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_10_again_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_10_again_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_10_again_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_10_again_total_score_num_hits, protocol_10_again_pseudo_interface_energy_num_hits, protocol_10_again_total_score_glycan_rmsd_count_2, protocol_10_again_total_score_glycan_rmsd_count_1_point_5, protocol_10_again_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_10_again_data )
protocol_10_again_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_10_again_data )
print "protocol_10_again\n\n\n"



############
#### 12 ####
############
# glycan_modeling_protocols/protocol_12
path_to_using_native_full_glycan_protocol_12 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_12.csv"
using_native_full_glycan_protocol_12_data = pd.read_csv( path_to_using_native_full_glycan_protocol_12 )
#path_to_using_native_full_glycan_protocol_12_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_12_no_reset.csv"
#using_native_full_glycan_protocol_12_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_12_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_12_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_12_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_12_data[ "total_score" ] - using_native_full_glycan_protocol_12_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_12"
print using_native_full_glycan_protocol_12_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_12_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_12"
print using_native_full_glycan_protocol_12_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_12"
print using_native_full_glycan_protocol_12_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_12_data[ using_native_full_glycan_protocol_12_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] >= 100 ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :5 ]

metrics = list( using_native_full_glycan_protocol_12_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_12_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_12_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_12_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_12_data, metric, "protocol_12" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_12_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_12_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_12_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_12_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_12_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_12_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_12_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_12_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_12_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_12_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_12_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_12_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_12_data[ "total_score" ], c=using_native_full_glycan_protocol_12_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_12_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_12_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_12_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_12_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_12_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_12_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_12_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_12_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_12_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_12_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_12_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_12_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_12"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_12_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_12_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_12_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_12_total_score_num_hits, protocol_12_pseudo_interface_energy_num_hits, protocol_12_total_score_glycan_rmsd_count_2, protocol_12_total_score_glycan_rmsd_count_1_point_5, protocol_12_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_12_data )
protocol_12_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_12_data )
print "protocol_12\n\n\n"



############
#### 13 ####
############
# glycan_modeling_protocols/protocol_13
path_to_using_native_full_glycan_protocol_13 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_13.csv"
using_native_full_glycan_protocol_13_data = pd.read_csv( path_to_using_native_full_glycan_protocol_13 )
#path_to_using_native_full_glycan_protocol_13_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_13_no_reset.csv"
#using_native_full_glycan_protocol_13_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_13_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_13_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_13_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_13_data[ "total_score" ] - using_native_full_glycan_protocol_13_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_13"
print using_native_full_glycan_protocol_13_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_13_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_13"
print using_native_full_glycan_protocol_13_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_13"
print using_native_full_glycan_protocol_13_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_13_data[ ( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ] >= 4 ) & ( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ] <= 6 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :5 ]

metrics = list( using_native_full_glycan_protocol_13_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_13_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_13_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_13_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_13_data, metric, "protocol_13" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_13_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_13_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_13_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_13_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_13_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_13_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_13_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_13_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_13_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_13_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_13_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_13_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_13_data[ "total_score" ], c=using_native_full_glycan_protocol_13_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_13_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_13_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_13_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_13_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_13_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_13_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_13_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_13_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_13_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_13_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_13_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_13_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_13"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_13_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_13_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_13_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_13_total_score_num_hits, protocol_13_pseudo_interface_energy_num_hits, protocol_13_total_score_glycan_rmsd_count_2, protocol_13_total_score_glycan_rmsd_count_1_point_5, protocol_13_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_13_data )
protocol_13_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_13_data )
print "protocol_13\n\n\n"



############
#### 14 ####
############
# glycan_modeling_protocols/protocol_14
path_to_using_native_full_glycan_protocol_14 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_14.csv"
using_native_full_glycan_protocol_14_data = pd.read_csv( path_to_using_native_full_glycan_protocol_14 )
#path_to_using_native_full_glycan_protocol_14_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_14_no_reset.csv"
#using_native_full_glycan_protocol_14_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_14_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_14_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_14_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_14_data[ "total_score" ] - using_native_full_glycan_protocol_14_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_14"
print using_native_full_glycan_protocol_14_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_14_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_14"
print using_native_full_glycan_protocol_14_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_14"
print using_native_full_glycan_protocol_14_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_14_data[ ( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ] >= 4 ) & ( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ] <= 6 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :5 ]

metrics = list( using_native_full_glycan_protocol_14_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_14_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_14_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_14_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_14_data, metric, "protocol_14" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_14_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_14_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_14_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_14_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_14_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_14_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_14_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_14_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_14_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_14_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_14_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_14_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_14_data[ "total_score" ], c=using_native_full_glycan_protocol_14_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_14_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_14_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_14_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_14_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_14_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_14_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_14_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_14_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_14_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_14_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_14_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_14_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_14"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_14_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_14_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_14_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_14_total_score_num_hits, protocol_14_pseudo_interface_energy_num_hits, protocol_14_total_score_glycan_rmsd_count_2, protocol_14_total_score_glycan_rmsd_count_1_point_5, protocol_14_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_14_data )
protocol_14_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_14_data )
print "protocol_14\n\n\n"



############
#### 15 ####
############
# glycan_modeling_protocols/protocol_15
path_to_using_native_full_glycan_protocol_15 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_15.csv"
using_native_full_glycan_protocol_15_data = pd.read_csv( path_to_using_native_full_glycan_protocol_15 )
#path_to_using_native_full_glycan_protocol_15_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_15_no_reset.csv"
#using_native_full_glycan_protocol_15_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_15_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_15_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_15_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_15_data[ "total_score" ] - using_native_full_glycan_protocol_15_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_15"
print using_native_full_glycan_protocol_15_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_15_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_15"
print using_native_full_glycan_protocol_15_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_15"
print using_native_full_glycan_protocol_15_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_15_data[ ( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ] >= 7 ) & ( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ] <= 7.5 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]

metrics = list( using_native_full_glycan_protocol_15_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_15_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_15_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_15_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_15_data, metric, "protocol_15" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_15_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_15_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_15_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_15_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_15_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_15_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_15_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_15_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_15_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_15_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_15_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_15_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_15_data[ "total_score" ], c=using_native_full_glycan_protocol_15_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_15_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_15_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_15_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_15_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_15_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_15_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_15_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_15_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_15_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_15_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_15_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_15_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_15"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_15_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_15_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_15_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_15_total_score_num_hits, protocol_15_pseudo_interface_energy_num_hits, protocol_15_total_score_glycan_rmsd_count_2, protocol_15_total_score_glycan_rmsd_count_1_point_5, protocol_15_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_15_data )
protocol_15_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_15_data )
print "protocol_15\n\n\n"



############
#### 16 ####
############
# glycan_modeling_protocols/protocol_16
path_to_using_native_full_glycan_protocol_16 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_16.csv"
using_native_full_glycan_protocol_16_data = pd.read_csv( path_to_using_native_full_glycan_protocol_16 )
#path_to_using_native_full_glycan_protocol_16_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_16_no_reset.csv"
#using_native_full_glycan_protocol_16_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_16_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_16_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_16_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_16_data[ "total_score" ] - using_native_full_glycan_protocol_16_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_16"
print using_native_full_glycan_protocol_16_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_16_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_16"
print using_native_full_glycan_protocol_16_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_16"
print using_native_full_glycan_protocol_16_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_16_data[ ( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ] >= 6 ) & ( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ] <= 8 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
grab_these_structs = []
#for jj in range( 1, 11 ):
#    low_E_struct = using_native_full_glycan_protocol_16_data[ ( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ] >= jj - 1 ) & ( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ] <= jj ) ].sort_values( "total_score", ascending=True )
#    print jj-1, jj, "mean:", np.mean( low_E_struct["total_score"] ), "median:", np.median( low_E_struct["total_score"] ), "min:", min( low_E_struct["total_score"] ), "max:", max( low_E_struct["total_score"] ), "range:", max( low_E_struct["total_score"] ) - min( low_E_struct["total_score"] )
#    print "  /home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_16/base_structs/%s" %low_E_struct.iloc[0]["filename"]
#    grab_these_structs.append( "/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_16/base_structs/protocol_%s.gz\n" %low_E_struct.iloc[0]["filename"] )
#with open( "protocol_16_lowest_E_of_1A_binned_rmsd", "wb" ) as fh:
#    fh.writelines( grab_these_structs )


metrics = list( using_native_full_glycan_protocol_16_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_16_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_16_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_16_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_16_data, metric, "protocol_16" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_16_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_16_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_16_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_16_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_16_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_16_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_16_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_16_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_16_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_16_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_16_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_16_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_16_data[ "total_score" ], c=using_native_full_glycan_protocol_16_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_16_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_16_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_16_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_16_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ -1, 10 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )
#plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_16_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_16_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_16_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_16_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_16_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_16_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_16_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_16_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_16"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_16_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_16_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_16_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_16_total_score_num_hits, protocol_16_pseudo_interface_energy_num_hits, protocol_16_total_score_glycan_rmsd_count_2, protocol_16_total_score_glycan_rmsd_count_1_point_5, protocol_16_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_16_data )
protocol_16_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_16_data )
print "protocol_16\n\n\n"



############
#### 17 ####
############
# glycan_modeling_protocols/protocol_17
path_to_using_native_full_glycan_protocol_17 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_17.csv"
using_native_full_glycan_protocol_17_data = pd.read_csv( path_to_using_native_full_glycan_protocol_17 )
#path_to_using_native_full_glycan_protocol_17_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_17_no_reset.csv"
#using_native_full_glycan_protocol_17_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_17_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_17_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_17_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_17_data[ "total_score" ] - using_native_full_glycan_protocol_17_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_17"
print using_native_full_glycan_protocol_17_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_17_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_17"
print using_native_full_glycan_protocol_17_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_17"
print using_native_full_glycan_protocol_17_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_17_data[ ( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ] >= 6 ) & ( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ] <= 8 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
grab_these_structs = []
#for jj in range( 1, 11 ):
#    low_E_struct = using_native_full_glycan_protocol_17_data[ ( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ] >= jj - 1 ) & ( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ] <= jj ) ].sort_values( "total_score", ascending=True )
#    print jj-1, jj, "mean:", np.mean( low_E_struct["total_score"] ), "median:", np.median( low_E_struct["total_score"] ), "min:", min( low_E_struct["total_score"] ), "max:", max( low_E_struct["total_score"] ), "range:", max( low_E_struct["total_score"] ) - min( low_E_struct["total_score"] )
#    print "  /home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_17/base_structs/%s" %low_E_struct.iloc[0]["filename"]
#    grab_these_structs.append( "/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_17/base_structs/protocol_%s.gz\n" %low_E_struct.iloc[0]["filename"] )
#with open( "protocol_17_lowest_E_of_1A_binned_rmsd", "wb" ) as fh:
#    fh.writelines( grab_these_structs )


metrics = list( using_native_full_glycan_protocol_17_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_17_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_17_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_17_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_17_data, metric, "protocol_17" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_17_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_17_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_17_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_17_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_17_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_17_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_17_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_17_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_17_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_17_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_17_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_17_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_17_data[ "total_score" ], c=using_native_full_glycan_protocol_17_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_17_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_17_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_17_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_17_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "total_score" )
#plt.ylim( [ -560, -520 ] )
plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_17_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_17_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_17_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_17_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_17_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_17_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_17_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_17_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_17"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_17_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_17_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_17_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_17_total_score_num_hits, protocol_17_pseudo_interface_energy_num_hits, protocol_17_total_score_glycan_rmsd_count_2, protocol_17_total_score_glycan_rmsd_count_1_point_5, protocol_17_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_17_data )
protocol_17_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_17_data )
print "protocol_17\n\n\n"


############
#### 18 ####
############
# glycan_modeling_protocols/protocol_18
path_to_using_native_full_glycan_protocol_18 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_18.csv"
using_native_full_glycan_protocol_18_data = pd.read_csv( path_to_using_native_full_glycan_protocol_18 )
#path_to_using_native_full_glycan_protocol_18_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_18_no_reset.csv"
#using_native_full_glycan_protocol_18_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_18_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_18_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_18_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_18_data[ "total_score" ] - using_native_full_glycan_protocol_18_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_18"
print using_native_full_glycan_protocol_18_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_18_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_18"
print using_native_full_glycan_protocol_18_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_18"
print using_native_full_glycan_protocol_18_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_18_data[ ( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ] >= 6 ) & ( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ] <= 8 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
grab_these_structs = []
#for jj in range( 1, 11 ):
#    low_E_struct = using_native_full_glycan_protocol_18_data[ ( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ] >= jj - 1 ) & ( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ] <= jj ) ].sort_values( "total_score", ascending=True )
#    print jj-1, jj, "mean:", np.mean( low_E_struct["total_score"] ), "median:", np.median( low_E_struct["total_score"] ), "min:", min( low_E_struct["total_score"] ), "max:", max( low_E_struct["total_score"] ), "range:", max( low_E_struct["total_score"] ) - min( low_E_struct["total_score"] )
#    print "  /home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_18/base_structs/%s" %low_E_struct.iloc[0]["filename"]
#    grab_these_structs.append( "/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_18/base_structs/protocol_%s.gz\n" %low_E_struct.iloc[0]["filename"] )
#with open( "protocol_18_lowest_E_of_1A_binned_rmsd", "wb" ) as fh:
#    fh.writelines( grab_these_structs )


metrics = list( using_native_full_glycan_protocol_18_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_18_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_18_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_18_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_18_data, metric, "protocol_18" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_18_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_18_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_18_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_18_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_18_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_18_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_18_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_18_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_18_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_18_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_18_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_18_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_18_data[ "total_score" ], c=using_native_full_glycan_protocol_18_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_18_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_18_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_18_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_18_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "total_score" )
#plt.ylim( [ -560, -520 ] )
plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_18_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_18_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_18_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_18_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_18_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_18_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_18_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_18_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_18"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_18_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_18_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_18_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_18_total_score_num_hits, protocol_18_pseudo_interface_energy_num_hits, protocol_18_total_score_glycan_rmsd_count_2, protocol_18_total_score_glycan_rmsd_count_1_point_5, protocol_18_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_18_data )
protocol_18_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_18_data )
print "protocol_18\n\n\n"



############
#### 19 ####
############
# glycan_modeling_protocols/protocol_19
path_to_using_native_full_glycan_protocol_19 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_19.csv"
using_native_full_glycan_protocol_19_data = pd.read_csv( path_to_using_native_full_glycan_protocol_19 )
#path_to_using_native_full_glycan_protocol_19_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_19_no_reset.csv"
#using_native_full_glycan_protocol_19_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_19_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_19_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_19_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_19_data[ "total_score" ] - using_native_full_glycan_protocol_19_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_19"
print using_native_full_glycan_protocol_19_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_19_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_19"
print using_native_full_glycan_protocol_19_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_19"
print using_native_full_glycan_protocol_19_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_19_data[ ( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ] >= 6 ) & ( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ] <= 8 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
grab_these_structs = []
#for jj in range( 1, 11 ):
#    low_E_struct = using_native_full_glycan_protocol_19_data[ ( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ] >= jj - 1 ) & ( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ] <= jj ) ].sort_values( "total_score", ascending=True )
#    print jj-1, jj, "mean:", np.mean( low_E_struct["total_score"] ), "median:", np.median( low_E_struct["total_score"] ), "min:", min( low_E_struct["total_score"] ), "max:", max( low_E_struct["total_score"] ), "range:", max( low_E_struct["total_score"] ) - min( low_E_struct["total_score"] )
#    print "  /home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_19/base_structs/%s" %low_E_struct.iloc[0]["filename"]
#    grab_these_structs.append( "/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_19/base_structs/protocol_%s.gz\n" %low_E_struct.iloc[0]["filename"] )
#with open( "protocol_19_lowest_E_of_1A_binned_rmsd", "wb" ) as fh:
#    fh.writelines( grab_these_structs )


metrics = list( using_native_full_glycan_protocol_19_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_19_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_19_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_19_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_19_data, metric, "protocol_19" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_19_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_19_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_19_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_19_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_19_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_19_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_19_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_19_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_19_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_19_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_19_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_19_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_19_data[ "total_score" ], c=using_native_full_glycan_protocol_19_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_19_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_19_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_19_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_19_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "total_score" )
#plt.ylim( [ -560, -520 ] )
plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_19_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_19_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_19_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_19_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_19_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_19_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_19_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_19_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_19"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_19_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_19_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_19_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_19_total_score_num_hits, protocol_19_pseudo_interface_energy_num_hits, protocol_19_total_score_glycan_rmsd_count_2, protocol_19_total_score_glycan_rmsd_count_1_point_5, protocol_19_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_19_data )
protocol_19_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_19_data )
print "protocol_19\n\n\n"



############
#### 20 ####
############
# glycan_modeling_protocols/protocol_20
path_to_using_native_full_glycan_protocol_20 = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_20.csv"
using_native_full_glycan_protocol_20_data = pd.read_csv( path_to_using_native_full_glycan_protocol_20 )
#path_to_using_native_full_glycan_protocol_20_no_reset = "/Users/mlnance/pyrosetta_dir/metric_data/glycan_modeling_protocols/using_native_full_glycan_protocol_20_no_reset.csv"
#using_native_full_glycan_protocol_20_no_reset_data = pd.read_csv( path_to_using_native_full_glycan_protocol_20_no_reset )

low_E_native_pseudo_interface_energy = np.mean( using_native_full_glycan_protocol_20_data[ "pseudo_interface_energy" ] - using_native_full_glycan_protocol_20_data[ "delta_pseudo_interface_energy" ] )
low_E_native_total_score = np.mean( using_native_full_glycan_protocol_20_data[ "total_score" ] - using_native_full_glycan_protocol_20_data[ "delta_total_score" ] )

print "***Top10 total_score protocol_20"
print using_native_full_glycan_protocol_20_data.sort_values( "total_score", ascending=True )[ :10 ].sort_values( "total_score" )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ]
#print using_native_full_glycan_protocol_20_data.sort_values( "total_score", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
print "***Top 10 glycan_rmsd protocol_20"
print using_native_full_glycan_protocol_20_data.sort_values( "Fc_glycan_rmsd", ascending=True )[ [ "total_score", "Fc_glycan_rmsd", "filename", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
print "***Top 10 Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A protocol_20"
print using_native_full_glycan_protocol_20_data.sort_values( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", ascending=False )[ [ "total_score", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A", "Fc_glycan_rmsd", "filename" ] ][ :10 ]
#print using_native_full_glycan_protocol_20_data[ ( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ] >= 6 ) & ( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ] <= 8 ) ].sort_values( "total_score", ascending=True )[ [ "filename", "total_score", "Fc_glycan_rmsd", "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ] ][ :10 ]
grab_these_structs = []
#for jj in range( 1, 11 ):
#    low_E_struct = using_native_full_glycan_protocol_20_data[ ( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ] >= jj - 1 ) & ( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ] <= jj ) ].sort_values( "total_score", ascending=True )
#    print jj-1, jj, "mean:", np.mean( low_E_struct["total_score"] ), "median:", np.median( low_E_struct["total_score"] ), "min:", min( low_E_struct["total_score"] ), "max:", max( low_E_struct["total_score"] ), "range:", max( low_E_struct["total_score"] ) - min( low_E_struct["total_score"] )
#    print "  /home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_20/base_structs/%s" %low_E_struct.iloc[0]["filename"]
#    grab_these_structs.append( "/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/glycan_modeling_protocols/protocol_20/base_structs/protocol_%s.gz\n" %low_E_struct.iloc[0]["filename"] )
#with open( "protocol_20_lowest_E_of_1A_binned_rmsd", "wb" ) as fh:
#    fh.writelines( grab_these_structs )


metrics = list( using_native_full_glycan_protocol_20_data.columns.values )
r_squared_to_metric_dict = {}
log10_r_squared_to_metric_dict = {}
binned_r_squared_to_metric_dict = {}
for metric in metrics:
    if metric == "atom_pair_constraint" or metric == "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_10A":
#    if metric != "filename" and metric != "Fc_glycan_rmsd":
        ## check normality of data
        #z, pval = normaltest( using_native_full_glycan_protocol_20_data[ metric ] )
        #if pval >= 0.05:
        r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_20_data, metric )
        r_squared_to_metric_dict[ r**2 ] = metric

        #r = get_r_of_line_of_best_fit( using_native_full_glycan_protocol_20_data, metric, log10 )
        #log10_r_squared_to_metric_dict[ r**2 ] = metric

for metric in metrics:
    if metric != "filename" and not metric.startswith( "delta" ) and metric != "Fc_glycan_rmsd":
        binned_r = get_r_of_line_of_best_fit_binned_Fc_glycan_rmsd( using_native_full_glycan_protocol_20_data, metric, "protocol_20" )
        binned_r_squared_to_metric_dict[ binned_r**2 ] = metric


fig, ax = plt.subplots(figsize=(40,25))
plt.subplot( 221 )
x = 0.0
y = low_E_native_pseudo_interface_energy
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_20_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_20_no_reset_data[ "pseudo_interface_energy" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_20_data[ "pseudo_interface_energy" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_20_data[ "pseudo_interface_energy" ], c=using_native_full_glycan_protocol_20_data[ "sugar_bb" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_20_no_reset_data[ "pseudo_interface_energy" ])), floor(min(using_native_full_glycan_protocol_20_data[ "pseudo_interface_energy" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_20_data[ "pseudo_interface_energy" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_20_data[ "pseudo_interface_energy" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "pseudo_interface_energy" )
plt.ylim( [ -25, 0 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 222 )
x = 0.0
y = low_E_native_total_score
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
#sc = plt.scatter( using_native_full_glycan_protocol_20_no_reset_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_20_no_reset_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_20_data[ "total_score" ] )
#sc = plt.scatter( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ], using_native_full_glycan_protocol_20_data[ "total_score" ], c=using_native_full_glycan_protocol_20_data[ "delta_hbonds" ] )
#plt.colorbar(sc)
#ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_20_no_reset_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_20_data[ "total_score" ])) ]
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_20_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_20_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_rmsd" )
#plt.xlim( [ -1, 10 ] )
plt.xlim( [ -1, 20 ] )
plt.ylabel( "total_score" )
#plt.ylim( [ -560, -520 ] )
plt.ylim( [ -560, -460 ] )
#plt.ylim( [ min(ymins) - 5, ymax + 5 ] )

plt.subplot( 223 )
xmin = 0
xmax = 12
#xmin = floor( min( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ] ) )
#xmax = ceil( max( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ] ) )
bins = np.arange( xmin, xmax + 0.5, 0.5 )
plt.hist( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ], bins=bins, histtype="stepfilled" )
#plt.hist( using_native_full_glycan_protocol_20_data[ "Fc_glycan_rmsd" ], bins, histtype="stepfilled" )
#plt.xticks( bins )
plt.xlabel( "Fc_glycan_rmsd" )
plt.xlim( [ xmin, xmax ] )
plt.ylabel( "count" )

plt.subplot( 224 )
x = 100.0
y = low_E_native_total_score
#sc = plt.scatter( using_native_full_glycan_protocol_20_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_20_data[ "total_score" ], marker='v', linewidth='0', s=20, c="orange", clip_on=False )
sc = plt.scatter( x, y, marker='D', s=36, c="red", clip_on=False )
sc = plt.scatter( using_native_full_glycan_protocol_20_data[ "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" ], using_native_full_glycan_protocol_20_data[ "total_score" ] )
ymins = [ floor(y), floor(min(using_native_full_glycan_protocol_20_data[ "total_score" ])), floor(min(using_native_full_glycan_protocol_20_data[ "total_score" ])) ]
ymax = ceil( np.percentile(using_native_full_glycan_protocol_20_data[ "total_score" ], 20) )
plt.xlabel( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_8A" )
plt.xlim( [ 101, 0 ] )
plt.ylabel( "total_score" )
plt.ylim( [ -560, -520 ] )


# save the plot
plt.tight_layout()
plot_title = "protocol_20"
#plot_title = "fa_intra_rep 3ay4 using SSM-200 on Fc glycan with LCM reset, using ideal pop data with native omega, no pack, min before mc, ramp, am3, 3 mpt, Gal_5A_1A_tol cst"
plt.suptitle( plot_title, fontsize = 36 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()

# print data
#print_r_squared_data( using_native_full_glycan_protocol_20_data, r_squared_to_metric_dict, "linear" )
#print_r_squared_data( using_native_full_glycan_protocol_20_data, log10_r_squared_to_metric_dict, "log10" )
#print_r_squared_data( using_native_full_glycan_protocol_20_data, binned_r_squared_to_metric_dict, "binned rmsd" )
protocol_20_total_score_num_hits, protocol_20_pseudo_interface_energy_num_hits, protocol_20_total_score_glycan_rmsd_count_2, protocol_20_total_score_glycan_rmsd_count_1_point_5, protocol_20_total_score_glycan_rmsd_count_1 = get_top10_data( using_native_full_glycan_protocol_20_data )
protocol_20_total_score_Fnat_num_hits = get_top10_Fnat_data( using_native_full_glycan_protocol_20_data )
print "protocol_20\n\n\n"







###########################
#### COMPARE PROTOCOLS ####
###########################
# protocol_9_again_total_score_num_hits, protocol_9_again_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_9_again_data
# protocol_10_again_total_score_num_hits, protocol_10_again_pseudo_interface_energy_num_hits = using_native_full_glycan_protocol_10_again_data

names = ( "Protocol_9_again", "Protocol_10_again", "Protocol_12", "Protocol_13", "Protocol_14", "Protocol_15", "Protocol_16", "Protocol_17", "Protocol_18", "Protocol_19", "Protocol_20", )
total_score_protocol_data = ( protocol_9_again_total_score_num_hits, protocol_10_again_total_score_num_hits, protocol_12_total_score_num_hits, protocol_13_total_score_num_hits, protocol_14_total_score_num_hits, protocol_15_total_score_num_hits, protocol_16_total_score_num_hits, protocol_17_total_score_num_hits, protocol_18_total_score_num_hits, protocol_19_total_score_num_hits, protocol_20_total_score_num_hits, )
pseudo_interface_energy_protocol_data = ( protocol_9_again_pseudo_interface_energy_num_hits, protocol_10_again_pseudo_interface_energy_num_hits, protocol_12_pseudo_interface_energy_num_hits, protocol_13_pseudo_interface_energy_num_hits, protocol_14_pseudo_interface_energy_num_hits, protocol_15_pseudo_interface_energy_num_hits, protocol_16_pseudo_interface_energy_num_hits, protocol_17_pseudo_interface_energy_num_hits, protocol_18_pseudo_interface_energy_num_hits, protocol_19_pseudo_interface_energy_num_hits, protocol_20_pseudo_interface_energy_num_hits, )
total_score_glycan_rmsd_count_2 = ( protocol_9_again_total_score_glycan_rmsd_count_2, protocol_10_again_total_score_glycan_rmsd_count_2, protocol_12_total_score_glycan_rmsd_count_2, protocol_13_total_score_glycan_rmsd_count_2, protocol_14_total_score_glycan_rmsd_count_2, protocol_15_total_score_glycan_rmsd_count_2, protocol_16_total_score_glycan_rmsd_count_2, protocol_17_total_score_glycan_rmsd_count_2, protocol_18_total_score_glycan_rmsd_count_2, protocol_19_total_score_glycan_rmsd_count_2, protocol_20_total_score_glycan_rmsd_count_2, )
total_score_glycan_rmsd_count_1_point_5 = ( protocol_9_again_total_score_glycan_rmsd_count_1_point_5, protocol_10_again_total_score_glycan_rmsd_count_1_point_5, protocol_12_total_score_glycan_rmsd_count_1_point_5, protocol_13_total_score_glycan_rmsd_count_1_point_5, protocol_14_total_score_glycan_rmsd_count_1_point_5, protocol_15_total_score_glycan_rmsd_count_1_point_5, protocol_16_total_score_glycan_rmsd_count_1_point_5, protocol_17_total_score_glycan_rmsd_count_1_point_5, protocol_18_total_score_glycan_rmsd_count_1_point_5, protocol_19_total_score_glycan_rmsd_count_1_point_5, protocol_20_total_score_glycan_rmsd_count_1_point_5, )
total_score_glycan_rmsd_count_1 = ( protocol_9_again_total_score_glycan_rmsd_count_1, protocol_10_again_total_score_glycan_rmsd_count_1, protocol_12_total_score_glycan_rmsd_count_1, protocol_13_total_score_glycan_rmsd_count_1, protocol_14_total_score_glycan_rmsd_count_1, protocol_15_total_score_glycan_rmsd_count_1, protocol_16_total_score_glycan_rmsd_count_1, protocol_17_total_score_glycan_rmsd_count_1, protocol_18_total_score_glycan_rmsd_count_1, protocol_19_total_score_glycan_rmsd_count_1, protocol_20_total_score_glycan_rmsd_count_1, )

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
plot_title = "Comparing more protocol variations of 3ay4-glycan modeling"
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
plot_title = "Comparing more protocol glycan rmsd below threshold variations of 3ay4-glycan modeling"
fig.suptitle( plot_title, fontsize = 20 )
fig.subplots_adjust(top=0.93)
#fig.savefig( plot_title, dpi=120, transparent=True )
fig.savefig( plot_title, bbox_extra_artists=(leg,), bbox_inches="tight", dpi=120, transparent=True )
plt.close()
