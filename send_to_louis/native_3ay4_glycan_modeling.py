#!/usr/bin/python
__author__="morganlnance"



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
parser.add_argument("protocol_num", type=int, help="which protocol number do you want to run?")
parser.add_argument("--verbose", "-v", action="store_true", default=False, help="do you want the program to print out pose scores during the protocol?")
parser.add_argument("--zip_decoy", "-z", action="store_true", default=False, help="do you want to zip up the dumped decoy pdb?")
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
    load_pose, initialize_rosetta, native_Fc_glycan_nums_except_core_GlcNAc, \
    native_Fc_glycan_nums
from rosetta import MoveMap, PyMOL_Mover

# utility directory function
from file_mover_based_on_fasc import main as get_lowest_E_from_fasc


# initialize Rosetta
initialize_rosetta()

# PyMOL_Mover
pmm = PyMOL_Mover()
pmm.keep_history( True )

# load the input pose
native_pose = load_pose( input_args.native_pdb_file )
#native_pose.pdb_info().name( "native_pose" )
pmm.apply( native_pose )





#########################
#### JOB DISTRIBUTOR ####
#########################

# instantiate the proper Protocol_X object
from native_3ay4_glycan_modeling_protocol import Model3ay4Glycan
if input_args.protocol_num == 9:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_9 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_9 creation and argument setting
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
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 10:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_10 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_10 creation and argument setting
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
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

# Skipping 11 here because I had done 11 before (it was 10 but without the constraint)
elif input_args.protocol_num == 12:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_12 version
    #######################################################################
    #### CORE GlcNAc CAN MOVE IN THIS PROTOCOL BUT IT STARTS AT NATIVE ####
    #######################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_12 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 13:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_13 version
    ########################################################################################################
    #### CORE GlcNAc CAN MOVE IN THIS PROTOCOL BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    ########################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_13 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 14:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_14 version
    ###############################################
    #### CORE GlcNAc CAN MOVE IN THIS PROTOCOL ####
    ###############################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_14 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 15:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_15 version
    ####################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN ####
    ####################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_15 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 16:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_16 version
    ##################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN AND CORE GlcNAc CAN MOVE BUT STARTS AT NATIVE ####
    ##################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_16 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 17:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_17 version
    #############################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    #############################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_17 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 18:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_18 version
    #############################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN AND CORE GlcNAc CAN MOVE ####
    #############################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_18 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 19:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_19 version
    #############################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    #############################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_19 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.spin_using_ideal_omegas = True  # meaning use only 180, 60, -60. Don't sample within +/- something of these values
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 20:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_20 version
    #############################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    #############################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_20 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.spin_using_ideal_omegas = False  # meaning use 180, 60, -60 but sample within +/- 0-15 of these values
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = False
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 21:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_21 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_21 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = False
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 22:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_22 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_22 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 23:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_23 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_23 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    # testing if the constraint is necessary
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 24:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_24 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    # fa_sol is 10% of normal value
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0, "fa_sol" : 0.09375 } )

    # Protocol_24 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

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

# run the appropriate protocol
print "Running Protocol %s in a PyJobDistributor..." %input_args.protocol_num
while not jd.job_complete:
    # get a fresh pose object
    working_pose = native_pose.clone()
    # name to use specifically in PyMOL
    GlycanModelProtocol.pmm_name = "p%s_decoy%s" %( input_args.protocol_num, cur_decoy_num ) 
    # name to use when dumping the decoy. Should include full path
    working_pose.pdb_info().name( jd.current_name )
    # apply the protocol
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
                                              GlycanModelProtocol.metrics_dump_dir, 
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
    # zip up the decoy pose, if desired
    if input_args.zip_decoy:
        os.popen( "gzip %s" %working_pose.pdb_info().name() )

    # increment the decoy number counter
    cur_decoy_num += 1

    #for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
    #    print "Res num", res_num, "Reset phi:", GlycanModelProtocol.reset_pose.phi( res_num ), "end phi:", working_pose.phi( res_num )
    #    print "Res num", res_num, "Reset psi:", GlycanModelProtocol.reset_pose.psi( res_num ), "end psi:", working_pose.psi( res_num )
    #    print "Res num", res_num, "Reset omega:", GlycanModelProtocol.reset_pose.omega( res_num ), "end omega:", working_pose.omega( res_num )
    #    print

# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = decoy_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, GlycanModelProtocol.lowest_E_structs_dir, 10 )
