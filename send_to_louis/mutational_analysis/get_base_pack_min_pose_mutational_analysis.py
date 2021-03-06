#!/usr/bin/python
__author__ = "morganlnance"



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
from rosetta.core.scoring import fa_rep, fa_intra_rep

from antibody_functions import load_pose, \
    initialize_rosetta, make_pack_rotamers_mover, \
    hold_chain_and_res_designations_3ay4, make_RotamerTrialsMover

from file_mover_based_on_fasc import main as get_lowest_E_from_fasc
from get_pose_metrics_on_native import main as get_pose_metrics_on_native


# initialize Rosetta ( comes from antibody_functions )
initialize_rosetta()

# get the full path to the original native PDB filename
orig_pdb_filename_full_path = input_args.native_pdb_file
orig_pdb_filename = orig_pdb_filename_full_path.split( '/' )[-1]
orig_pdb_name = orig_pdb_filename.split( ".pdb" )[0]

# this is where the double pack and minimized versions of the native will be dumped
structure_dir = base_structs_dir
decoy_pdb_name = structure_dir + orig_pdb_name

# this is where the single pack and minimized decoy will be dumped (the predecessor to double)
# see long comment in jd code for clarification
round1_structure_dir = base_structs_dir + "/round1_decoys/"
if not os.path.isdir( round1_structure_dir ):
    os.mkdir( round1_structure_dir )
round1_decoy_pdb_name = round1_structure_dir + orig_pdb_name

# sets up the input native PDB as being the base pose
native_pose = Pose()
native_pose.assign( load_pose( orig_pdb_filename_full_path ) )
native_pose.pdb_info().name( "native" )

# automatically populates chain and residue information into holder for native 3ay4
native_pose_info = hold_chain_and_res_designations_3ay4()
native_pose_info.native()


# fa_intra_rep should always 0.440 since that's what I've been using for sugars
sf = get_fa_scorefxn()
sf.set_weight( fa_intra_rep, 0.440 )
orig_fa_rep = sf.get_weight( fa_rep )


# pymol stuff
pmm = PyMOL_Mover()
pmm.keep_history( True )
pmm.apply( native_pose )


# relay information to user
info_file_details = []
info_file_details.append( "Native PDB filename:\t\t\t%s\n" %input_args.native_pdb_file.split( '/' )[-1] )
info_file_details.append( "Creating this many decoys:\t\t%s\n" %str( input_args.base_nstruct ) )
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
while not jd.job_complete:
    # make a working pose by cloning a fresh native_pose
    working_pose = native_pose.clone()
    working_pose.pdb_info().name( "base_%s" %str( jd.current_num ) )
    pmm.apply( working_pose )

    # instantiate the 3ay4 information holder class object
    # this is for calculating metrics
    # see antibody_functions for more information on this hard-coded function
    working_pose_info = hold_chain_and_res_designations_3ay4()
    working_pose_info.native()

    # collect all non-branch point residue numbers
    moveable_residues = [ res_num for res_num in range( 1, working_pose.n_residue() + 1 ) if working_pose.residue( res_num ).is_branch_point() == False ]


    # running this twice because the mutational analysis SHOULDN'T have to adjust the native input pose each time it gets run.
    # this means I need to dump two decoys per turn. The first pack/min round dumps the decoy that would be used as the native
    # for each call of the mutation script. This is because the lowE "native" would get mutated then pack/min in this way again
    # meaning, at the end of the day, it got pack/min twice. So as a far comparison, a base lowE native decoy with no mutations
    # should have been pack/min twice. For good measure, I don't want things that were sent on different paths/seeds, so the
    # starting native will be the first round dump of the double pack/min decoy that got the lowest score. Therefore, dump
    # round 1 decoy into one directory and dump round 2 decoy into the main structure dir. Only the final round 2 decoy will
    # have a fasc file, which is fine

    
    for ii in range( 2 ):
        # pack all residues except for branch point residues
        pack_rotamers_mover = make_RotamerTrialsMover( moveable_residues = moveable_residues, 
                                                       sf = sf,
                                                       input_pose = working_pose,
                                                       pack_radius = None )
        pack_rotamers_mover.apply( working_pose )
        pmm.apply( working_pose )

        # we packed the side chains, so minimize them (only residues marked by moveable_residues)
        # keep the backbone the same as the crystal because 1) we can, 2) it's easier, 3) we want to see what mutations do to packing more so
        mm = MoveMap()
        mm.set_bb( False )
        for res_num in moveable_residues:
            mm.set_chi( res_num, True )
        # gradient min of native input_pose
        for jj in range( 3 ):
            if jj == 0:
                sf.set_weight( fa_rep, orig_fa_rep * 0.1 )
            elif jj == 1:
                sf.set_weight( fa_rep, orig_fa_rep * 0.33 )
            elif jj == 2:
                sf.set_weight( fa_rep, orig_fa_rep )
            min_mover = MinMover( movemap_in = mm,
                                  scorefxn_in = sf,
                                  min_type_in = "lbfgs_armijo_nonmonotone",
                                  tolerance_in = 0.001,
                                  use_nb_list_in = True )
            min_mover.max_iter( 2500 )
            min_mover.apply( working_pose )
            print "working_pose", sf( working_pose ), jj
            pmm.apply( working_pose )

        # dump the round1 decoy into the appropriate directory
        if ii == 0:
            dump_name = round1_structure_dir + jd.current_name.split( '/' )[-1]
            working_pose.dump_file( dump_name )

    # inform user of decoy number
    print "\tFinished with decoy %s" %str( jd.current_num )

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
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 10 )
