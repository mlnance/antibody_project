#!/usr/bin/python
__author__ = "morganlnance"



#########################
#### PARSE ARGUMENTS ####
#########################

import argparse

# parse and store input arguments
parser = argparse.ArgumentParser(description="Use PyRosetta to glycosylate a pose and find a low E structure")
#parser.add_argument("argument_file", type=str, help="/path/to/the protocol argument file")
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure.")
#parser.add_argument("glyco_file", type=str, help="/path/to/the .iupac glycan file to be used.")
parser.add_argument("utility_dir", type=str, help="where do your utility files live? Give me the directory.")
parser.add_argument("structure_dir", type=str, help="where do you want to dump the decoys made during this protocol?")
parser.add_argument("nstruct", type=int, help="how many decoys do you want to make using this protocol?")
#parser.add_argument("scorefxn_file", type=str, help="which scorefxn weights do you want to use on top of a standard full atom scorefunction?")
parser.add_argument("protocol_num", type=int, help="which protocol number from native_3ay4_glycan_modeling_protocols do you want to run?")
parser.add_argument("--verbose", "-v", action="store_true", default=False, help="do you want the program to print out pose scores during the protocol?")
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
from native_3ay4_glycan_modeling_protocol_functions import get_fa_scorefxn_with_given_weights, \
    load_pose, initialize_rosetta, native_Fc_glycan_nums_except_core_GlcNAc
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

# collect and create necessary directories for use in metric calculations
working_dir = os.getcwd() + '/'
metrics_dump_dir = working_dir + "metrics_dump_dir"
try:
    os.mkdir( metrics_dump_dir )
except:
    pass




#########################
#### JOB DISTRIBUTOR ####
#########################

# instantiate the proper Protocol_X object
from native_3ay4_glycan_modeling_protocols import Model3ay4Glycan
if input_args.protocol_num == 0:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_0 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_0 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_5A_1A_tol.cst"
    GlycanModelProtocol.verbose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 1:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_1 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_1 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_5A_1A_tol.cst"
    GlycanModelProtocol.verbose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 2:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_2 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_2 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 3:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_3 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_3 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 100
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 4:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_4 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_4 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 100
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )


elif input_args.protocol_num == 5:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_5 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_5 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

# else I haven't made this protocol number yet
else:
    print "\nI haven't created the protocol number you gave me yet.\n"
    sys.exit()

# create an appropriate decoy_name using the Protocol_X.dump_dir and input_args.protocol_num
decoy_name = GlycanModelProtocol.base_structs_dir + "protocol_%s_decoy" %input_args.protocol_num



##########################
#### PYJOBDISTRIBUTOR ####
##########################

# imports
from rosetta import PyJobDistributor

# create and use the PyJobDistributor object
jd = PyJobDistributor( decoy_name, input_args.nstruct, sf )
jd.native_pose = native_pose
cur_decoy_num = 1

print "Running Protocol %s in a PyJobDistributor..." %input_args.protocol_num
while not jd.job_complete:
    # run the appropriate protocol
    working_pose = native_pose.clone()
    working_pose.pdb_info().name( "p%s_decoy%s" %( input_args.protocol_num, cur_decoy_num ) )
    working_pose.assign( GlycanModelProtocol.apply( working_pose ) )

    # collect additional metric data
    try:
        # all this is here until I update get_pose_metrics(_on_native)
        from antibody_functions import hold_chain_and_res_designations_3ay4
        from get_pose_metrics_on_native import main as get_pose_metrics_on_native
        native_pose_info = hold_chain_and_res_designations_3ay4()
        native_pose_info.native()
        working_pose_info = hold_chain_and_res_designations_3ay4()
        working_pose_info.native()
        # metric calculations
        metrics = get_pose_metrics_on_native( working_pose, 
                                              working_pose_info, 
                                              native_pose, 
                                              native_pose_info, 
                                              sf, 
                                              2, # Fc-FcR interface JUMP_NUM
                                              jd.current_num, 
                                              metrics_dump_dir, 
                                              input_args.utility_dir, 
                                              MC_acceptance_rate = GlycanModelProtocol.mc_acceptance, 
                                              native_constraint_file = GlycanModelProtocol.constraint_file )
    except:
        metrics = ''
        pass

    # add the metric data to the .fasc file
    jd.additional_decoy_info = metrics

    # dump the decoy
    jd.output_decoy( working_pose )
    cur_decoy_num += 1
