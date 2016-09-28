#!/usr/bin/python
__author__ = "morganlnance"



#########################
#### PARSE ARGUMENTS ####
#########################

import argparse

# parse and store input arguments
parser = argparse.ArgumentParser(description="Use PyRosetta to glycosylate a pose and find a low E structure")
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure.")
#parser.add_argument("glyco_file", type=str, help="/path/to/the .iupac glycan file to be used.")
parser.add_argument("utility_dir", type=str, help="where do your utility files live? Give me the directory.")
parser.add_argument("structure_dir", type=str, help="where do you want to dump the decoys made during this protocol?")
parser.add_argument("nstruct", type=int, help="how many decoys do you want to make using this protocol?")
parser.add_argument("scorefxn_file", type=str, help="which scorefxn weights do you want to use on top of a standard full atom scorefunction?")
parser.add_argument("protocol_num", type=int, help="which protocol number from native_3ay4_glycan_modeling_protocols do you want to run?")
parser.add_argument("--verbose", "-v", action="store_true", default=False, help="do you want the program to print out pose scores during the protocol?")
parser.add_argument("--watch_for_convergence", "-c", action="store_true", default=False, help="do you want the program to print out pose scores using a static sf during the protocol? So you can watch if you are converging")
input_args = parser.parse_args()



##########################
#### CHECK ALL INPUTS ####
##########################

# check for validity of file paths
import os, sys

# check the utility directory
if not os.path.isdir( input_args.utility_dir ):
    print "\nYour utility_dir argument( %s ) does not exist. Please check your input. Exiting." %input_args.utility_dir
    sys.exit()

# add the utility directory to the system path for loading of modules
sys.path.append( input_args.utility_dir )


## check the validity of the passed arguments
# make sure the structure_dir passed is valid
if os.path.isdir( input_args.structure_dir ):
    if not input_args.structure_dir.endswith( '/' ):
        main_structure_dir = input_args.structure_dir + '/'
    else:
        main_structure_dir = input_args.structure_dir
else:
    print
    print "It seems that the structure_dir argument you gave me ( %s ) does not exist. Please check your input or create this directory before running this protocol." %input_args.structure_dir
    sys.exit()

# make sure the files passed actually exist
if input_args.native_pdb_file is not None:
    if not os.path.isfile( input_args.native_pdb_file ):
        print "\nYour native_pdb_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.native_pdb_file
        sys.exit()
#if input_args.glyco_file is not None:
#    if not os.path.isfile( input_args.glyco_file ):
#        print "\nYour glyco_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.glyco_file
#        sys.exit()



##################################
#### CREATE NECESSARY OBJECTS ####
##################################

# imports 
from antibody_functions import initialize_rosetta, load_pose, \
    make_fa_scorefxn_from_file, native_Fc_glycan_nums_except_core_GlcNAc
from rosetta import MoveMap, PyMOL_Mover

# initialize Rosetta
initialize_rosetta()

# PyMOL_Mover
pmm = PyMOL_Mover()
pmm.keep_history( True )

# load the input pose
native_pose = load_pose( input_args.native_pdb_file )
#native_pose.pdb_info().name( "native_pose" )
pmm.apply( native_pose )

# create the desired scorefxn
main_sf = make_fa_scorefxn_from_file( input_args.scorefxn_file )





#########################
#### JOB DISTRIBUTOR ####
#########################

# instantiate the proper Protocol_X object
if input_args.protocol_num == 0:
    from native_3ay4_glycan_modeling_protocols import Protocol_0

    # create the necessary minimization and overall movement MoveMap for Protocol_0
    min_mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        min_mm.set_bb( res_num, True )
        min_mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            min_mm.set_branches( res_num, False )

    # Protocol_0
    protocol = Protocol_0( min_mm, main_sf, input_args.structure_dir, input_args.utility_dir )
    protocol.write_protocol_info_file( native_pose )
    protocol.verbose = True

# create an appropriate decoy_name using the protocol.dump_dir and input_args.protocol_num
decoy_name = protocol.base_structs_dir + "protocol_%s_decoy" %input_args.protocol_num



##########################
#### PYJOBDISTRIBUTOR ####
##########################

# imports
from rosetta import PyJobDistributor

# create and use the PyJobDistributor object
jd = PyJobDistributor( decoy_name, input_args.nstruct, main_sf )
jd.native_pose = native_pose
cur_decoy_num = 1

print "Running Protocol %s in a PyJobDistributor..." %input_args.protocol_num
while not jd.job_complete:
    # run the appropriate protocol
    testing_pose = native_pose.clone()
    testing_pose.pdb_info().name( "p%s_decoy%s" %( input_args.protocol_num, cur_decoy_num ) )
    testing_pose.assign( protocol.apply( testing_pose ) )

    # collect additional metric data
    try:
        metrics = get_pose_metrics_on_native( testing_pose, 
                                              testing_pose_info, 
                                              native_pose, 
                                              native_pose_info, 
                                              main_sf, 
                                              2, # Fc-FcR interface JUMP_NUM
                                              jd.current_num, 
                                              metrics_dump_dir, 
                                              input_args.utility_dir, 
                                              MC_acceptance_rate = mc_acceptance, 
                                              native_constraint_file = input_args.native_constraint_file )
    except:
        metrics = ''
        pass

    # add the metric data to the .fasc file
    jd.additional_decoy_info = metrics

    # dump the decoy
    jd.output_decoy( testing_pose )
    cur_decoy_num += 1
