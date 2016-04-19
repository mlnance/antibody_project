#!/usr/bin/python
__author__ = "morganlnance"


import argparse

parser = argparse.ArgumentParser(description="Use Rosetta to calculate RMSD between a native pose and a directory of structures")
parser.add_argument("native_pdb_filename", help="the filename of the PDB structure to serve as the native structure")
parser.add_argument("structure_dir", help="where do the structures to which I am comparing the native live?")
parser.add_argument("resulting_filename", type=str, help="what do you want the resulting csv file to be called? This program will add the .csv extension for you")
input_args = parser.parse_args()



#################
#### IMPORTS ####
#################

from antibody_functions import initialize_rosetta, load_pose
from rosetta import Pose, get_fa_scorefxn
from rosetta.core.scoring import CA_rmsd

import os, sys, csv
try:
    import pandas as pd
except:
    pass



#######################
#### RMSD PROTOCOL ####
#######################

# check the structure directory
working_dir = os.getcwd() + '/'
if os.path.isdir( input_args.structure_dir ):
    if not input_args.structure_dir.endswith( '/' ):
        structure_dir = input_args.structure_dir + '/'
    else:
        structure_dir = input_args.structure_dir
else:
    print "It seems like you gave me an incorrect path, exiting"
    sys.exit()

# get all the structure names from the structure directory
os.chdir( structure_dir )
structures = []
structure_names = []
for f in os.listdir( os.getcwd() ):
    if f.endswith( ".pdb" ):
        structures.append( os.path.abspath( f ) )
        structure_names.append( f.split( '/' )[-1] )
os.chdir( working_dir )

# inform the user of the structure directory and number of files to be analyzed
num_structs = len( structure_names )
print "Analyzing", num_structs, "structures from", structure_dir
print


# check and load native pose
initialize_rosetta()
try:
    if os.path.isfile( input_args.native_pdb_filename ):
        native = Pose()
        native.assign( load_pose( input_args.native_pdb_filename ) )
except:
    print "It appears", input_args.native_pdb_filename, "is not a valid pdb file. Exiting"
    sys.exit()
    

# make a scorefunction
sf = get_fa_scorefxn()

# collect the data
rmsds = []
scores = []
pdb_names = []

# for if Pandas doesn't work
header = [ "pdb_names", "score", "rmsd" ]
mutant_df_data = []
mutant_df_data.append( header )

# run the data collection while printing out decoy_num
decoy_num = 1
for pdb in structures:
    print "Working on decoy number", decoy_num, "of", num_structs
    
    # load the mutant Pose
    mutant = Pose()
    mutant.assign( load_pose( pdb ) )
    
    # collect the data
    name = pdb.split( '/' )[-1]
    score = sf( mutant )
    rmsd = CA_rmsd( native, mutant )
    
    # for if Pandas does work
    pdb_names.append( name )
    rmsds.append( rmsd )
    scores.append( score )
    
    # for if Pandas doesn't work
    mutant_data = [ name, score, rmsd ]
    mutant_df_data.append( mutant_data )
    
    # up the decoy_num counter
    decoy_num += 1


# check out the passed result filename for .csv extension
filename = input_args.resulting_filename
if not filename.endswith( ".csv" ):
    filename += ".csv"

# dump the data into a DataFrame
try:
    df = pd.DataFrame( index=pdb_names )
    df["score"] = scores
    df["rmsd"] = rmsds
    print df
    df.to_csv( filename, index=True, index_label = "pdb_name" )

# or, dump the three lists into three separate .csv files
except:
    with open( filename, "wb" ) as f:
        writer = csv.writer( f )
        writer.writerows( mutant_df_data )
