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
parser.add_argument("--scorefxn_file", default=None, type=str, help="/path/to/the .sf scorefxn space-delimited file that tells me which scoring weights beyond the norm you want to use")
parser.add_argument("--fa_intra_rep", action="store_true", help="do you want to set fa_intra_rep to 0.440 in the sf?")
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

# collect and create necessary directories for use in metric calculations
working_dir = os.getcwd() + '/'
metrics_dump_dir = working_dir + "base_metrics_dir"
try:
    os.mkdir( metrics_dump_dir )
except:
    pass



################################
#### INITIAL PROTOCOL SETUP ####
################################

# good to go, import needed functions
from rosetta import Pose, get_fa_scorefxn, PyJobDistributor, \
    PyMOL_Mover, MoveMap, MinMover
from rosetta.core.scoring import score_type_from_name

from antibody_functions import load_pose, \
    initialize_rosetta, apply_sugar_constraints_to_sf, \
    make_pack_rotamers_mover, make_fa_scorefxn_from_file, \
    hold_chain_and_res_designations_3ay4

from file_mover_based_on_fasc import main as get_lowest_E_from_fasc
from get_pose_metrics_on_native import main as get_pose_metrics_on_native


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
native_pose.pdb_info().name( "native" )

# automatically populates chain and residue information into holder for native 3ay4
native_pose_info = hold_chain_and_res_designations_3ay4()
native_pose_info.native()


# use the scorefxn_file to set up additional weights
if input_args.scorefxn_file is not None:
    sf = make_fa_scorefxn_from_file( input_args.scorefxn_file )
# else create a fa_scorefxn
else:
    sf = get_fa_scorefxn()

# fa_intra_rep should always 0.440 since that's what I've been using
if input_args.fa_intra_rep:
    sf.set_weight( score_type_from_name( "fa_intra_rep" ), 0.440 )


# pymol stuff
pmm = PyMOL_Mover()
pmm.keep_history( True )
pmm.apply( native_pose )


# relay information to user
info_file_details = []
info_file_details.append( "Native PDB filename:\t\t\t%s\n" %input_args.native_pdb_file.split( '/' )[-1] )
info_file_details.append( "Creating this many decoys:\t\t%s\n" %str( input_args.base_nstruct ) )
info_file_details.append( "ScoreFunction file used?:\t\t%s\n" %str( input_args.scorefxn_file ).split( '/' )[-1] )
info_file_details.append( "Main structure directory:\t\t%s\n" %main_structure_dir )
info_file_details.append( "Base structure directory:\t\t%s\n" %base_structs_dir )
info_file_details.append( "Lowest E structure directory:\t\t%s\n" %lowest_E_structs_dir )
info_file_details.append( "\nScore weights used in sf:\n%s\n" %( "\n".join( [ "%s: %s" %( str( name ), sf.get_weight( name ) ) for name in sf.get_nonzero_weighted_scoretypes() ] ) ) )
info_file = ''.join( info_file_details )
print "\n", info_file, "\n"

# write out the info file with the collected info from above
info_filename = main_structure_dir + "protocol_run.info"
with open( info_filename, "wb" ) as fh:
    fh.write( "Info for this run of %s\n\n" %__file__ )
    fh.write( info_file )




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
    # make a working pose
    working_pose = Pose()
    working_pose.assign( native_pose )
    working_pose.pdb_info().name( "base_%s" %str( decoy_num ) )
    pmm.apply( working_pose )

    # instantiate the 3ay4 information holder class object
    working_pose_info = hold_chain_and_res_designations_3ay4()

    # see antibody_functions for more information on this hard-coded function
    working_pose_info.native()


    for ii in range( 2 ):
        # pack
        pack_rotamers_mover = make_pack_rotamers_mover( sf, working_pose,
                                                        apply_sf_sugar_constraints = False,
                                                        pack_branch_points = True )
        pack_rotamers_mover.apply( working_pose )

        # minimize
        mm = MoveMap()
        mm.set_bb( True )
        mm.set_chi( True )
        mm.set_branches( True )

        min_mover = MinMover( movemap_in = mm,
                              scorefxn_in = sf,
                              min_type_in = "dfpmin_strong_wolfe",
                              tolerance_in = 0.01,
                              use_nb_list_in = True )
        min_mover.apply( working_pose )

    pmm.apply( working_pose )

    # inform user of decoy number
    print "\tFinished with decoy %s" %str( decoy_num )
    decoy_num += 1

    # collect additional metric data
    try:
        metrics = get_pose_metrics_on_native( working_pose,
                                              working_pose_info,
                                              native_pose,
                                              native_pose_info,
                                              sf,
                                              2, # Fc-FcR interface JUMP_NUM
                                              jd.current_num,
                                              metrics_dump_dir,
                                              input_args.utility_directory )
    except:
        metrics = ''
        pass

    # add the metric data to the .fasc file
    jd.additional_decoy_info = metrics
    
    # dump the decoy
    jd.output_decoy( working_pose )
    
# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = decoy_pdb_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 5 )
