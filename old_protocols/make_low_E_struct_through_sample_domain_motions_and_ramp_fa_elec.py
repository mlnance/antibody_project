#!/usr/bin/python

# this takes a pose and runs through only domain motion rounds to make a low-energy structure
# using this to take the native crystal structure of 3ay4 and do nstruct of 1000 to get a low energy structure
# this low energy structure will then be the base structure to make mutations, and then this protocol will be run on all of the mutations
import argparse

# parse and store args
parser = argparse.ArgumentParser(description="Use PyRosetta to sample domain motions")
parser.add_argument("pdb_filename", type=str, help="the filename of the PDB structure to be minimized")
parser.add_argument("structure_directory", type=str, help="where do you want your decoys to be dumped? Each PDB will have its own directory there")
parser.add_argument("--nstruct", type=int, default=1000, help="how many decoy structures do you want to create? Default = 1000")
parser.add_argument("--trials", type=int, default=100, help="how many trials do you want to complete during each call of the protocol? Default is 100.")
input_args = parser.parse_args()


import os
import sys

# check to see that all of the files exist (PDB and loops file)
if not os.path.isfile( input_args.pdb_filename ):
    print input_args.pdb_filename, "does not exist, exiting"
    sys.exit()

## check the structure directory
# if the path is not a directory, exit
if not os.path.isdir( input_args.structure_directory ):
    print input_args.structure_directory, "is not a directory, exiting"
    sys.exit()
# otherwise, if it is a directory, add a '/' to the end of the path if needed
else:
    if input_args.structure_directory.endswith( '/' ):
        main_structure_dir = input_args.structure_directory
    else:
        main_structure_dir = input_args.structure_directory + '/'


# good to go, import needed functions
from antibody_protocols import *
from rosetta import PyJobDistributor

# get the full path to the PDB and the filename
orig_pdb_filename_full_path = input_args.pdb_filename
orig_pdb_filename = orig_pdb_filename_full_path.split( '/' )[-1]
orig_pdb_name = orig_pdb_filename.split( ".pdb" )[0]

# create structure directory where decoys will be dumped
structure_dir = main_structure_dir + orig_pdb_name
if not os.path.isdir( structure_dir ):
    os.mkdir( structure_dir )
decoy_pdb_name = structure_dir + '/' + orig_pdb_name

# sets up the input pdb as being the base pose
native_pose = load_pose( orig_pdb_filename_full_path )
working_pose = load_pose( orig_pdb_filename_full_path )

# make the necessary score function
sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, working_pose )


# create and use the PyJobDistributor
jd = PyJobDistributor( decoy_pdb_name, input_args.nstruct, sf )
jd.native_pose = native_pose

# pymol watcher
pmm.keep_history( True )

while not jd.job_complete:
    # sample domain motions
    working_pose = make_rigid_body_moves( sf, working_pose, trials = input_args.trials, verbose = True )
    
    jd.output_decoy( working_pose )
    
