#!/usr/bin/python
__author__ = "morganlnance"

'''
Plans for this code:
1) take in a low-energy native PDB structure that has been packed and minimized
2) make the mutation(s) specified
3) run a pack and minimization around the mutation site within 20 Angstroms
5) dump the mutant pose in mut_structs
'''




import argparse

# parse and store args
parser = argparse.ArgumentParser(description="Use PyRosetta mutate the given Pose with mutations specified by a file")
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure")
parser.add_argument("structure_directory", type=str, help="where do you want your decoys to be dumped? Each PDB will have its own directory there")
parser.add_argument("mutations_file", help="give me a file containing all mutations you want to make.")
parser.add_argument("utility_directory", type=str, help="where do the utility files live? Give me the directory.")
input_args = parser.parse_args()



###################################
#### CHECK ALL INPUT ARGUMENTS ####
###################################

import os, sys

# check to see that all of the files exist
# check the native file
if not os.path.isfile( input_args.native_pdb_file ):
    print "Your argument", input_args.native_pdb_file, "for native_pdb_file does not exist, exiting"
    sys.exit()

# check the mutations file
if not os.path.isfile( input_args.mutations_file ):
    print "Your argument", input_args.mutations_file, "for mutations_file does not exist. exiting"
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

# get the full path to the original native PDB filename
orig_pdb_filename_full_path = input_args.native_pdb_file
orig_pdb_filename = orig_pdb_filename_full_path.split( '/' )[-1]
orig_pdb_name = orig_pdb_filename.split( ".pdb" )[0]

# make the needed mutant structure directory
mut_structs_dir = main_structure_dir + "mutations_of_" + orig_pdb_name
if not os.path.isdir( mut_structs_dir ):
    os.mkdir( mut_structs_dir )


# relay information to user
info_file_details = []
info_file_details.append( "Native PDB filename:\t\t%s\n" %input_args.native_pdb_file.split( '/' )[-1] )
info_file_details.append( "Main structure directory:\t%s\n" %main_structure_dir )
info_file_details.append( "Mut structure directory:\t%s\n" %mut_structs_dir )
info_file = ''.join( info_file_details )
print "\n", info_file, "\n"

# write out the info file with the collected info from above
info_filename = main_structure_dir + "protocol_run.info"
with open( info_filename, "wb" ) as fh:
    fh.write( "Info for this run of %s\n\n" %__file__ )
    fh.write( info_file )




################################
#### INITIAL PROTOCOL SETUP ####
################################

# good to go, import needed functions
from antibody_functions import load_pose, initialize_rosetta, \
    read_3ay4_mutation_file, mutate_residue
from rosetta import Pose, get_fa_scorefxn


# initialize Rosetta ( comes from antibody_functions )
initialize_rosetta()

# sets up the input native PDB as being the base pose
native_pose = Pose()
native_pose.assign( load_pose( orig_pdb_filename_full_path ) )

# make the necessary score function
sf = get_fa_scorefxn()



##################################
#### MUTANT POSE CONSTRUCTION ####
##################################

# read in the mutations to be made from the mutations_file
all_mutations = read_3ay4_mutation_file( input_args.mutations_file )

# for each set of mutations, make and dump a mutant pose
for mutations in all_mutations:
    # get a fresh copy of the native Pose
    mutant = Pose()
    mutant.assign( native_pose )

    # split on '+' to ensure you are making each single point mutation on this Pose
    mutations = mutations.split( '+' )

    # make each single point mutation
    for mutation_str in mutations:
        # make the mutation depending on how the mutation was designated (with a chain id or not)
        # here a chain id was designated
        if len( mutation_str.split( '_' ) ) == 2:
            # get the mutation and the chain id based on how the string should be structured
            # mutation_chainid such as A123T_A
            mutation = mutation_str.split( '_' )[ 0 ]
            chain_id = mutation_str.split( '_' )[ 1 ]

            # get the original amino acid, the PDB number, the new amino acid, and the chain id from the mutation
            # the original amino acid is always the first letter
            orig_AA = mutation[ 0 ].upper()
            # the new amino acid is always at the end
            new_AA = mutation[ -1 ].upper()
            # the pdb res num is always in the middle
            pdb_num = int( mutation[ 1 : -1 ] )
            pose_num = native_pose.pdb_info().pdb2pose( chain_id, pdb_num )

            # ensure the orginal amino acid specified is the same as the one already in the Pose
            # if not, skip the mutation and print to screen
            if orig_AA != native_pose.residue( pose_num ).name1():
                print "Hold up! What you said was the original amino acid is actually incorrect!!"
                print "You told me there was originally a", orig_AA, "at PDB position", pdb_num, "chain", chain_id, "but there actually was a", native_pose.residue( pose_num ).name1(), ". Check your input. Exiting."
                pass
            # if they're the same, continue
            else:
                mutant.assign( mutate_residue( pdb_num, new_AA, mutant, sf, pdb_num = True, pdb_chain = chain_id, pack_radius = 20 ) )

        # else, no chain id. This is a symmetrical mutant on chain A and B
        else:
            # get the mutation id based on how the string should be structured
            # mutation such as A123T
            mutation = mutation_str

            # get the original amino acid, the PDB number, and the new amino acid id from the mutation
            # the original amino acid is always the first letter
            orig_AA = mutation[ 0 ].upper()
            # the new amino acid is always at the end
            new_AA = mutation[ -1 ].upper()
            # the pdb res num is always in the middle
            pdb_num = int( mutation[ 1 : -1 ] )

            # the pose_num is collected from chain A and chain B
            pose_num1 = native_pose.pdb_info().pdb2pose( 'A', pdb_num )
            pose_num2 = native_pose.pdb_info().pdb2pose( 'B', pdb_num )

            # ensure the orginal amino acid specified is the same as the two already in the Pose
            # if not, skip the mutation and print to screen
            if orig_AA != native_pose.residue( pose_num1 ).name1():
                print "Hold up! What you said was the original amino acid is actually incorrect!!"
                print "You told me there was originally a", orig_AA, "at PDB position", pdb_num, "chain A, but there actually was a", native_pose.residue( pose_num1 ).name1(), ". Check your input. Exiting."
                pass
            elif orig_AA != native_pose.residue( pose_num2 ).name1():
                print "Hold up! What you said was the original amino acid is actually incorrect!!"
                print "You told me there was originally a", orig_AA, "at PDB position", pdb_num, "chain B, but there actually was a", native_pose.residue( pose_num2 ).name1(), ". Check your input. Exiting."
                pass
            # if they're the same, continue
            else:
                mutant.assign( mutate_residue( pdb_num, new_AA, mutant, sf, pdb_num = True, pdb_chain = 'A', pack_radius = 20 ) )
                mutant.assign( mutate_residue( pdb_num, new_AA, mutant, sf, pdb_num = True, pdb_chain = 'B', pack_radius = 20 ) )

    # update the mutant's name
    mutant_name = '+'.join( mutations )
    mutant.pdb_info().name( mutant_name )

    # dump the mutant Pose
    dump_name = mut_structs_dir + '/' + mutant_name + ".pdb"
    mutant.dump_file( dump_name )
