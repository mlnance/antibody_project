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
### PURELY FOR TESTING HOW MANY MOVES NEED TO HAPPEN BEFORE THE TOTAL E FLATTENS ####
if input_args.protocol_num == 100:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_100 version
    ###########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG1 Fc DATA ####
    ###########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_100 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.outer_trials = 1
    GlycanModelProtocol.inner_trials = 20
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    # default.table now includes IgG1 Fc stats, so this shouldn't need to be used
    GlycanModelProtocol.set_native_core_omegas_to_stats = False
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
    GlycanModelProtocol.make_move = True
    GlycanModelProtocol.watch_E_vs_trial = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

### PURELY FOR TESTING HOW MANY MOVES NEED TO HAPPEN BEFORE THE TOTAL E FLATTENS ####
elif input_args.protocol_num == 101:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_101 version
    ###########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG1 Fc DATA ####
    ###########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_101 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.outer_trials = 1
    GlycanModelProtocol.inner_trials = 40
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    # default.table now includes IgG1 Fc stats, so this shouldn't need to be used
    GlycanModelProtocol.set_native_core_omegas_to_stats = False
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
    GlycanModelProtocol.make_move = True
    GlycanModelProtocol.watch_E_vs_trial = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

### PURELY FOR TESTING HOW MANY MOVES NEED TO HAPPEN BEFORE THE TOTAL E FLATTENS ####
elif input_args.protocol_num == 102:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_102 version
    ###########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG1 Fc DATA ####
    ###########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_102 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.outer_trials = 1
    GlycanModelProtocol.inner_trials = 80
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    # default.table now includes IgG1 Fc stats, so this shouldn't need to be used
    GlycanModelProtocol.set_native_core_omegas_to_stats = False
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
    GlycanModelProtocol.make_move = True
    GlycanModelProtocol.watch_E_vs_trial = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

### PURELY FOR TESTING HOW MANY MOVES NEED TO HAPPEN BEFORE THE TOTAL E FLATTENS ####
elif input_args.protocol_num == 103:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_103 version
    ###########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG1 Fc DATA ####
    ###########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_103 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.outer_trials = 1
    GlycanModelProtocol.inner_trials = 160
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    # default.table now includes IgG1 Fc stats, so this shouldn't need to be used
    GlycanModelProtocol.set_native_core_omegas_to_stats = False
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
    GlycanModelProtocol.make_move = True
    GlycanModelProtocol.watch_E_vs_trial = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

### PURELY FOR TESTING HOW MANY MOVES NEED TO HAPPEN BEFORE THE TOTAL E FLATTENS ####
elif input_args.protocol_num == 104:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_104 version
    ###########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG1 Fc DATA ####
    ###########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_104 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.outer_trials = 1
    GlycanModelProtocol.inner_trials = 320
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    # default.table now includes IgG1 Fc stats, so this shouldn't need to be used
    GlycanModelProtocol.set_native_core_omegas_to_stats = False
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
    GlycanModelProtocol.make_move = True
    GlycanModelProtocol.watch_E_vs_trial = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

### PURELY FOR TESTING HOW MANY MOVES NEED TO HAPPEN BEFORE THE TOTAL E FLATTENS ####
elif input_args.protocol_num == 105:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_105 version
    ###########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG1 Fc DATA ####
    ###########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_105 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.outer_trials = 1
    GlycanModelProtocol.inner_trials = 640
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    # default.table now includes IgG1 Fc stats, so this shouldn't need to be used
    GlycanModelProtocol.set_native_core_omegas_to_stats = False
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
    GlycanModelProtocol.make_move = True
    GlycanModelProtocol.watch_E_vs_trial = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

### PURELY FOR TESTING HOW MANY MOVES NEED TO HAPPEN BEFORE THE TOTAL E FLATTENS ####
elif input_args.protocol_num == 106:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_106 version
    ###########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG1 Fc DATA ####
    ###########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_106 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.outer_trials = 1
    GlycanModelProtocol.inner_trials = 1280
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    # default.table now includes IgG1 Fc stats, so this shouldn't need to be used
    GlycanModelProtocol.set_native_core_omegas_to_stats = False
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
    GlycanModelProtocol.make_move = True
    GlycanModelProtocol.watch_E_vs_trial = True

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
Fc_glycan_rmsd_lowest_seen = None
protocol_movie_poses = []
protocol_movie_decoy_num = 0

# run the appropriate protocol
print "Running Protocol %s in a PyJobDistributor..." %input_args.protocol_num
while not jd.job_complete:
    # get a fresh pose object
    working_pose = native_pose.clone()
    # name to use specifically in PyMOL
    GlycanModelProtocol.pmm_name = "p%s_decoy%s" %( input_args.protocol_num, jd.current_num ) 
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

    # extra work
    # check out the movie poses and keep track of those corresponding to the structure with the lowest rmsd
    if GlycanModelProtocol.make_movie:
        # hacky way of pulling out Fc_glycan_rmsd from the extra metrics I calculate
        # the number following Fc_glycan_rmsd: in the metrics string is the actual Fc_glycan_rmsd
        metrics_split = metrics.split( ' ' )
        Fc_glycan_rmsd = float( metrics_split[ metrics_split.index( "Fc_glycan_rmsd:" ) + 1 ] )
        # if the Fc_glycan_rmsd for this decoy is the lowest seen, store the Fc_glycan_rmsd and keep its movie_poses
        # update the Fc_glycan_lowest_seen if this is the first decoy done
        if Fc_glycan_rmsd_lowest_seen is None:
            Fc_glycan_rmsd_lowest_seen = Fc_glycan_rmsd
            protocol_movie_decoy_num = jd.current_num
            # grab the poses associated with this decoy's movie
            protocol_movie_poses = GlycanModelProtocol.movie_poses
        # if this decoy is lower in Fc_glycan_rmsd than the lowest seen, update it
        elif Fc_glycan_rmsd < Fc_glycan_rmsd_lowest_seen:
            Fc_glycan_rmsd_lowest_seen = Fc_glycan_rmsd
            protocol_movie_decoy_num = jd.current_num
            # grab the poses associated with this decoy's movie
            protocol_movie_poses = GlycanModelProtocol.movie_poses    

    # dump the decoy
    jd.output_decoy( working_pose )
    # zip up the decoy pose, if desired
    if input_args.zip_decoy:
        os.popen( "gzip %s" %working_pose.pdb_info().name() )

    #for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
    #    print "Res num", res_num, "Reset phi:", GlycanModelProtocol.reset_pose.phi( res_num ), "end phi:", working_pose.phi( res_num )
    #    print "Res num", res_num, "Reset psi:", GlycanModelProtocol.reset_pose.psi( res_num ), "end psi:", working_pose.psi( res_num )
    #    print "Res num", res_num, "Reset omega:", GlycanModelProtocol.reset_pose.omega( res_num ), "end omega:", working_pose.omega( res_num )
    #    print


# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = decoy_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, GlycanModelProtocol.lowest_E_structs_dir, 10 )

# dump the movie poses, if desired
if GlycanModelProtocol.make_movie:
    # make a directory in the movie_poses_dir using the decoy number
    movie_dir = GlycanModelProtocol.movie_poses_dir + "decoy_%s_movie/" %protocol_movie_decoy_num
    os.mkdir( movie_dir )
    for pose in protocol_movie_poses:
        pose.dump_file( movie_dir + pose.pdb_info().name() )

sys.exit()
