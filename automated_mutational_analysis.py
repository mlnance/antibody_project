#!/usr/bin/python
__author__ = "morganlnance"

'''
Plans for this code:
1) take in a native PDB structure
2) run single pack and minimization with base_nstruct=1000, dumping the structures into base_structs
3) take lowest E from (2) and turn it into lowest_E_single_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb, and dump into into lowest_E_structs dir
4) for each mutant in the given mutant_list_file...
5) make the mutation and save the structure into mut_structs
6) run two pack and minimizations on the mutant with mut_nstruct=1000, dumping the structures into base_structs
'''

# this takes a pose and runs through only pack and minimization rounds to make a low-energy structure
# using this to take the native crystal structure of 3ay4 and do nstruct of 1000 to get a low energy structure
# this low energy structure will then be the base structure to make mutations, and then this protocol will be run on all of the mutations

import argparse

# parse and store args
parser = argparse.ArgumentParser(description="Use PyRosetta to pack and minimize a structure into a low-energy conformation")
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure")
parser.add_argument("structure_directory", type=str, help="where do you want your decoys to be dumped? Each PDB will have its own directory there")
parser.add_argument("mutant_list_file", help="give me the list of mutations you want to have made")
parser.add_argument("--base_native_file", default=None, help="if you have it, give the the path to the low E base native structure you want to use. Default = None = make a low E base native structure.")
parser.add_argument("--base_nstruct", type=int, default=1000, help="how many decoy structures do you want to create to get a base native structure? Default = 1000")
parser.add_argument("--mut_nstruct", type=int, default=1000, help="how many decoy structures do you want to create to get a low E mutant structure? Default = 1000")
input_args = parser.parse_args()



###################################
#### CHECK ALL INPUT ARGUMENTS ####
###################################

import os
import sys
sys.path.append( "/Users/Research/antibody_project/utility_functions" )

# check to see that all of the files exist
if not os.path.isfile( input_args.native_pdb_file ):
    print "Your argument", input_args.native_pdb_file, "for native_pdb_file does not exist, exiting"
    sys.exit()
if input_args.base_native_file is not None:
    if not os.path.isfile( input_args.base_native_file ):
        print "Your argument", input_args.base_native_file, "for base_native_file does not exist, exiting"
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

# make the needed directories if needed
# base_structs, lowest_E_structs, and mut_structs with
# input_args.structure_directory as the base directory
base_structs_dir = main_structure_dir + "base_structs/"
lowest_E_structs_dir = main_structure_dir + "lowest_E_structs/"
mut_structs_dir = main_structure_dir + "mut_structs/"

if not os.path.isdir( base_structs_dir ):
    os.mkdir( base_structs_dir )
if not os.path.isdir( lowest_E_structs_dir ):
    os.mkdir( lowest_E_structs_dir )
if not os.path.isdir( mut_structs_dir ):
    os.mkdir( mut_structs_dir )

# relay information to user
print
print "Native PDB filename:\t\t", input_args.native_pdb_file.split( '/' )[-1]
if input_args.base_native_file is not None:
    print "Low E Base native PDB:\t", input_args.base_native_file
print "Main structure directory:\t", main_structure_dir
print "Base structures directory:\t", base_structs_dir
print "Lowest E structures directory:\t", lowest_E_structs_dir
print "Mutant structures directory:\t", mut_structs_dir
print



################################
#### INITIAL PROTOCOL SETUP ####
################################

# good to go, import needed functions
from antibody_functions import load_pose, make_all_mutations, \
    initialize_rosetta, apply_sugar_constraints_to_sf, \
    make_my_new_symmetric_antibody, make_my_new_asymmetric_antibody
from antibody_protocols import make_base_pack_min_pose
from file_mover_based_on_fasc import main as get_lowest_E_from_fasc
from rosetta import Pose, get_fa_scorefxn, PyJobDistributor

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

# make the necessary score function
sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, native_pose )



#######################################
#### BASE NATIVE POSE CONSTRUCTION ####
#######################################

# if the user did not pass a base native PDB structure to use
# say for instance, from an earlier run of this program
if input_args.base_native_file is None:
    # create and use the PyJobDistributor
    jd = PyJobDistributor( decoy_pdb_name, input_args.base_nstruct, sf )
    jd.native_pose = native_pose
    
    # make base_nstruct of the native doing one pack and minimization to get a standard low E structure
    print "Making a low E base pose by packing and minimizing the passed native pose..."
    while not jd.job_complete:
        # grab a fresh copy of the native pose
        working_pose = Pose()
        working_pose.assign( native_pose )
        
        # pack and minimize
        working_pose.assign( make_base_pack_min_pose( sf, working_pose, trials = 2, verbose = True ) )
        
        # dump the decoy
        jd.output_decoy( working_pose )
        
    # move the lowest E pack and minimized native structure into the lowest_E_structs dir
    fasc_filename = decoy_pdb_name + ".fasc"
    lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 5 )
    
# otherwise, they already gave a lowest_E_native_filename
else:
    lowest_E_native_filename = input_args.base_native_file



##################################
#### MUTANT POSE CONSTRUCTION ####
##################################

# load the lowest_E_native_filename as the new native pose
low_E_native_pose = Pose()
low_E_native_pose.assign( load_pose( lowest_E_native_filename ) )

# ensure the mutation list filename is accurate and can be opened
try:
    mutant_lines = []
    f = open( input_args.mutant_list_file, 'rb' )
    lines = f.readlines()

    # for each line specifying a mutation
    for line in lines:
        # take off the carriage return
        line = line.rstrip()
        # skip over comments
        if line != '' and line[0] != '#':
            # mutation specifications should be split by whitespace
            symmetry = line.split( ' ' )[1]
            
            # check the symmetry - if it's not 'sym' or 'asym,' exit
            # '' is fine for a symmetry designation - it defaults to 'sym'
            if symmetry == '' or symmetry == "sym" or symmetry == "asym":
                mutant_lines.append( line )
            else:
                print "ERROR in the following line:", line
                print "'%s'" %symmetry, "isn't a valid a symmetrical designation. Please put 'sym' or 'asym'"
                print "Exiting"
                sys.exit()
except:
    print
    print
    raise


# for each filename, make mut_nstruct decoys by pack and minimization
for mut_line in mutant_lines:
    # pull apart the needed information from the mutation lines
    mut_line = mut_line.split( ' ' )
    mutation = mut_line[ 0 ]
    symmetry = mut_line[ 1 ]
    
    # make the directory for the mutant PDB in the base_structs_dir
    structure_dir = base_structs_dir + mutation
    if not os.path.isdir( structure_dir ):
        os.mkdir( structure_dir )
    mut_decoy_pdb_name = structure_dir + '/' + mutation
    
    # make the necessary score function
    sf = get_fa_scorefxn()
    sf = apply_sugar_constraints_to_sf( sf, low_E_native_pose )
    
    # create and use the PyJobDistributor
    jd = PyJobDistributor( mut_decoy_pdb_name, input_args.mut_nstruct, sf )
    jd.native_pose = low_E_native_pose
    
    # make mut_nstruct of the native doing one pack and minimization to get a standard low E structure
    print "Making %s decoys of a %s %s mutant..." %( input_args.mut_nstruct, symmetry, mutation )
    while not jd.job_complete:
        # grab a fresh copy of the mutant native pose
        mut_working_pose = Pose()
        mut_working_pose.assign( low_E_native_pose )
        
        if symmetry == '' or symmetry == "sym":
            mut_working_pose.assign( make_my_new_symmetric_antibody( mutation, sf, mut_working_pose,
                                                                     pack_around_mut = True,
                                                                     dump_pose = False,
                                                                     dump_dir = None ) )
        elif symmetry == "asym":
            mut_working_pose.assign( make_my_new_asymmetric_antibody( mutation, sf, mut_working_pose,
                                                                      pack_around_mut = True,
                                                                      dump_pose = False,
                                                                      dump_dir = None ) )
            
        # dump the decoy
        jd.output_decoy( mut_working_pose )

    # move the lowest E pack and minimized mutant structure into the lowest_E_structs dir
    fasc_filename = mut_decoy_pdb_name + ".fasc"
    lowest_E_mutant_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 5 )
