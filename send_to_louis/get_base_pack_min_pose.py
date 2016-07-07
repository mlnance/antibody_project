#!/usr/bin/python
__author__ = "morganlnance"

'''
Plans for this code:
1) take in a native PDB structure
2) run single pack and minimization with base_nstruct=1000, dumping the structures into base_structs
3) take lowest E from (2) and turn it into lowest_E_single_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb, and dump into into lowest_E_structs dir
'''



import argparse

# parse and store args
parser = argparse.ArgumentParser(description="Use PyRosetta to pack and minimize a structure into a low-energy conformation")
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure")
parser.add_argument("structure_directory", type=str, help="where do you want your decoys to be dumped? Each PDB will have its own directory there")
parser.add_argument("utility_directory", type=str, help="where do the utility files live? Give me the directory.")
parser.add_argument("base_nstruct", type=int, help="how many decoy structures do you want to create to get a base native structure?")
input_args = parser.parse_args()



###################################
#### CHECK ALL INPUT ARGUMENTS ####
###################################

import os
import sys

# check to see that all of the files exist
if not os.path.isfile( input_args.native_pdb_file ):
    print "Your argument", input_args.native_pdb_file, "for native_pdb_file does not exist, exiting"
    sys.exit()

# check the structure directory
if not os.path.isdir( input_args.structure_directory ):
    print "Your argument", input_args.structure_directory, "for structure_directory is not a directory, exiting"
    sys.exit()
else:
    if input_args.structure_directory.endswith( '/' ):
        main_structure_dir = input_args.structure_directory
    else:
        main_structure_dir = input_args.structure_directory + '/'

# check the utility directory
if not os.path.isdir( input_args.utility_directory ):
    print "Your argument", input_args.utility_directory, "for utility_directory is not a directory, exiting"
    sys.exit()

# add the utility directory to the system path for loading of modules
sys.path.append( input_args.utility_directory )

# make the needed directories if needed
# base_structs and lowest_E_structs
# input_args.structure_directory as the base directory
base_structs_dir = main_structure_dir + "base_structs/"
lowest_E_structs_dir = main_structure_dir + "lowest_E_structs/"

if not os.path.isdir( base_structs_dir ):
    os.mkdir( base_structs_dir )
if not os.path.isdir( lowest_E_structs_dir ):
    os.mkdir( lowest_E_structs_dir )

# relay information to user
print
print "Native PDB filename:\t\t", input_args.native_pdb_file.split( '/' )[-1]
print "Main structure directory:\t", main_structure_dir
print "Base structures directory:\t", base_structs_dir
print "Lowest E structures directory:\t", lowest_E_structs_dir
print



################################
#### INITIAL PROTOCOL SETUP ####
################################

# good to go, import needed functions
from antibody_functions import load_pose, make_all_mutations, \
    initialize_rosetta, apply_sugar_constraints_to_sf
from antibody_protocols import make_base_pack_min_pose
from file_mover_based_on_fasc import main as get_lowest_E_from_fasc
from rosetta import Pose, get_fa_scorefxn, PyJobDistributor, PyMOL_Mover

# initialize Rosetta ( comes from antibody_functions )
initialize_rosetta()

# get the full path to the original native PDB filename
orig_pdb_filename_full_path = input_args.native_pdb_file
orig_pdb_filename = orig_pdb_filename_full_path.split( '/' )[-1]
orig_pdb_name = orig_pdb_filename.split( ".pdb" )[0]

# make the directory for the native PDB in the base_structs_dir
# this is where packed and minimized versions of the native will lie
structure_dir = base_structs_dir + orig_pdb_name
if not os.path.isdir( structure_dir ):
    os.mkdir( structure_dir )
decoy_pdb_name = structure_dir + '/' + orig_pdb_name

# sets up the input native PDB as being the base pose
native_pose = Pose()
native_pose.assign( load_pose( orig_pdb_filename_full_path ) )
native_pose.pdb_info().name( "base" )

# make the necessary score function
sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, native_pose )

# pymol stuff
pmm = PyMOL_Mover()
pmm.keep_history( True )
pmm.apply( native_pose )



#######################################
#### BASE NATIVE POSE CONSTRUCTION ####
#######################################

# create and use the PyJobDistributor
jd = PyJobDistributor( decoy_pdb_name, input_args.base_nstruct, sf )
jd.native_pose = native_pose

# make base_nstruct of the native doing one pack and minimization to get a standard low E structure
print "Making a low E base pose by packing and minimizing the passed native pose..."
decoy_num = 1
while not jd.job_complete:
    # grab a fresh copy of the native pose
    fresh_native = Pose()
    fresh_native.assign( native_pose )

    # make a working pose
    working_pose = Pose()
    
    # pack and minimize
    working_pose.assign( make_base_pack_min_pose( sf, fresh_native, trials = 2, verbose = True, pmm = pmm ) )
    
    # inform user of decoy number
    print "\tFinished with decoy %s" %str( decoy_num )
    decoy_num += 1
    
    # dump the decoy
    jd.output_decoy( working_pose )
    
# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = decoy_pdb_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 5 )
