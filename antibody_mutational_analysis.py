#!/usr/bin/python
__author__ = "morganlnance"

'''
Plans for this code:
1) take in a low-energy pose
2) for each mutant in the given mutant_list_file...
3) make the mutation and save the structure into mut_structs
4) run one pack and minimization for residues within 20 Angstroms of the mutation with mut_nstruct=1000, dumping the structures into base_structs
'''



import argparse

# parse and store args
parser = argparse.ArgumentParser(description="Use PyRosetta to make point mutations in an already pack and minimized pose.")
parser.add_argument("low_E_pdb_file", type=str, help="the filename of the lowest-energy native PDB file.")
parser.add_argument("structure_directory", type=str, help="where do you want your decoys to be dumped? Each PDB will have its own directory there")
parser.add_argument("mutant_list_file", help="give me the list of mutations you want to have made")
parser.add_argument("utility_directory", type=str, help="where do the utility files live? Give me the directory.")
parser.add_argument("mut_nstruct", type=int, help="how many decoy structures do you want to create to get a low E mutant structure?")
input_args = parser.parse_args()



###################################
#### CHECK ALL INPUT ARGUMENTS ####
###################################

import os
import sys

# check to see that all of the files exist
if not os.path.isfile( input_args.low_E_pdb_file ):
    print "Your argument", input_args.low_E_pdb_file, "for low_E_pdb_file does not exist, exiting"
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
print "Native PDB filename:\t\t", input_args.low_E_pdb_file.split( '/' )[-1]
print "Main structure directory:\t", main_structure_dir
print "Base structures directory:\t", base_structs_dir
print "Lowest E structures directory:\t", lowest_E_structs_dir
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
low_E_pdb_filename_full_path = input_args.low_E_pdb_file
low_E_pdb_filename = low_E_pdb_filename_full_path.split( '/' )[-1]
low_E_pdb_name = low_E_pdb_filename.split( ".pdb" )[0]

# load the lowest_E_native_filename as the new native pose
low_E_native_pose = Pose()
low_E_native_pose.assign( load_pose( low_E_pdb_filename_full_path ) )

# make the necessary score function
sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, low_E_native_pose )



##################################
#### MUTANT POSE CONSTRUCTION ####
##################################

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
