#!/usr/bin/python
__author__ = "morganlnance"

'''
NOTES
STARTING POSE (3ay4 without Fc) is a base structure that was acquired from making 1000 pack/min decoys of native 3ay4 (with Fc sugar on) and then manually removing the Fc sugars of the lowest energy decoy
'''

'''
SAMPLE INPUT
run glycan_sampling_just_sugar_small_moves_on_native.py project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII_removed_Fc_sugar.pdb ~/antibody_project/carbohydrate_files/3ay4_Fc_Glycan.iupac project_utility_files/ ~/pyrosetta_dir/test_pdb_dir/ 1 10 5
'''


# global values
kT = 0.8


#########################
#### PARSE ARGUMENTS ####
#########################

import argparse

# parse and store input arguments
parser = argparse.ArgumentParser(description="Use PyRosetta to glycosylate a pose and find a low E structure")
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure.")
parser.add_argument("glyco_file", type=str, help="/path/to/the .iupac glycan file to be used.")
parser.add_argument("utility_dir", type=str, help="where do your utility files live? Give me the directory.")
parser.add_argument("structure_dir", type=str, help="where do you want to dump the decoys made during this protocol?")
parser.add_argument("nstruct", type=int, help="how many decoys do you want to make using this protocol?")
parser.add_argument("num_sugar_small_move_trials", type=int, help="how many SugarSmallMoves do you want to make within the Fc glycan?")
parser.add_argument("num_moves_per_trial", type=int, help="how many SugarSmallMoves do you want to make within one trial?")
parser.add_argument("--random_reset", action="store_true", help="do you want to randomly reset the phi, psi, and omega values of the Fc glycan? (Excluding core GlcNAc)")
parser.add_argument("--ramp_sf", action="store_true", help="do you want to ramp up the fa_atr term and ramp down the fa_rep term?")
parser.add_argument("--native_constraint_file", default=None, type=str, help="/path/to/the .cst constraint file you want to use for the protocol")
parser.add_argument("--scorefxn_file", default=None, type=str, help="/path/to/the .sf scorefxn space-delimited file that tells me which scoring weights beyond the norm you want to use")
parser.add_argument("--angle_multiplier", type=float, default=1.0, help="by how much do you want to increase the angle_max for the SugarSmallMover? Default is 1.0 (standard value)" )
parser.add_argument("--verbose", "-v", action="store_true", default=False, help="do you want the program to print out pose scores during the protocol?")
input_args = parser.parse_args()



##########################
#### CHECK ALL INPUTS ####
##########################

# check for validity of file paths
import os, sys, random

# check the utility directory
if not os.path.isdir( input_args.utility_dir ):
    print "Your argument", input_args.utility_dir, "for utility_directory is not a directory, exiting"
    sys.exit()

# add the utility directory to the system path for loading of modules
sys.path.append( input_args.utility_dir )


## collect and create necessary directories for use in metric calculations
# use a unique id as suffix to metric dir as often many instances are run at the same time
from util import id_generator

working_dir = os.getcwd() + '/'
metrics_dir_id = str( id_generator() )
metrics_dump_dir = working_dir + "sugar_small_%s_dir_%s" %( str( input_args.num_sugar_small_move_trials ), metrics_dir_id )
try:
    os.mkdir( metrics_dump_dir )
except:
    pass

## check the validity of the passed arguments
# make sure the structure_dir passed is valid
if os.path.isdir( input_args.structure_dir ):
    if not input_args.structure_dir.endswith( '/' ):
        main_structure_dir = input_args.structure_dir + '/'
    else:
        main_structure_dir = input_args.structure_dir
else:
    print
    print "It seems that the directory you gave me ( %s ) does not exist. Please check your input or create this directory before running this protocol." %input_args.structure_dir
    sys.exit()

# make the needed directories if needed
# base_structs and lowest_E_structs
# input_args.structure_directory as the base directory
base_structs_dir = main_structure_dir + "base_structs/"
lowest_E_structs_dir = main_structure_dir + "lowest_E_structs/"

if not os.path.isdir( base_structs_dir ):
    try:
        os.mkdir( base_structs_dir )
    except:
        pass
if not os.path.isdir( lowest_E_structs_dir ):
    try:
        os.mkdir( lowest_E_structs_dir )
    except:
        pass



#################
#### IMPORTS ####
#################

# Rosetta worker functions
from rosetta import Pose, pose_from_file, get_fa_scorefxn, \
    PyMOL_Mover, MonteCarlo, PyJobDistributor, MoveMap
from rosetta.numeric.random import random_range, uniform
from rosetta.core.scoring import score_type_from_name
from rosetta.protocols.simple_moves import ConstraintSetMover
from rosetta import SmallMover, MinMover
from toolbox import get_hbonds

# Rosetta functions I wrote out
from antibody_functions import initialize_rosetta, \
    get_fa_scorefxn_with_given_weights, make_pack_rotamers_mover, \
    load_pose, ramp_score_weight, make_fa_scorefxn_from_file, \
    SugarSmallMover, hold_chain_and_res_designations_3ay4

# utility functions
from file_mover_based_on_fasc import main as get_lowest_E_from_fasc
from get_pose_metrics_on_native import main as get_pose_metrics_on_native



################################
#### INITIAL PROTOCOL SETUP ####
################################

# initialize Rosetta
initialize_rosetta()

# load up the poses given from the arguments passed
native_pose = Pose()
native_pose.assign( load_pose( input_args.native_pdb_file ) )

# automatically populates chain and residue information into holder for native 3ay4
native_pose_info = hold_chain_and_res_designations_3ay4()
native_pose_info.native()

# get the full path of the native PDB name
native_pdb_filename_full_path = input_args.native_pdb_file
native_pdb_filename = native_pdb_filename_full_path.split( '/' )[-1]
native_pdb_name = native_pdb_filename.split( ".pdb" )[0]

# change the name of the native PDB name
native_pose_name = "native_pose"
native_pose.pdb_info().name( native_pose_name )


# make the directory for the working PDBs in the base_structs_dir
structure_dir = base_structs_dir + native_pdb_name
if not os.path.isdir( structure_dir ):
    try:
        os.mkdir( structure_dir )
    except:
        pass
working_pose_decoy_name = structure_dir + '/' + native_pdb_name + "_just_%s_sugar_small_moves" %input_args.num_sugar_small_move_trials


# use the scorefxn_file to set up additional weights
if input_args.scorefxn_file is not None:
    main_sf = make_fa_scorefxn_from_file( input_args.scorefxn_file )
# else create a fa_scorefxn
else:
    main_sf = get_fa_scorefxn()

# fa_intra_rep should always be 0.440 since that's what I've been using
main_sf.set_weight( score_type_from_name( "fa_intra_rep" ), 0.440 )


# set up constraints from the passed constraint file
if input_args.native_constraint_file is not None:
    # add an appropriate weight to the main_sf if atom_pair_constraint is 0
    if main_sf.get_weight( score_type_from_name( "atom_pair_constraint" ) ) == 0:
        main_sf.set_weight( score_type_from_name( "atom_pair_constraint" ), 1.0 )

    if input_args.verbose:
        print "Setting up a ConstraintSetMover"

    # create ConstraintSetMover
    constraint_file_used = True
    constraint_setter = ConstraintSetMover()
    constraint_setter.constraint_file( input_args.native_constraint_file )
else:
    constraint_file_used = False


# pymol stuff
pmm = PyMOL_Mover()
pmm.keep_history( True )
pmm.apply( native_pose )



# relay information to user
info_file_details = []
info_file_details.append( "Native PDB filename:\t\t\t%s\n" %input_args.native_pdb_file.split( '/' )[-1] )
info_file_details.append( "Sugar filename:\t\t\t\t%s\n" %input_args.glyco_file.split( '/' )[-1] )
info_file_details.append( "Creating this many decoys:\t\t%s\n" %str( input_args.nstruct ) )
info_file_details.append( "Number of SugarSmallMove trials:\t%s\n" %str( input_args.num_sugar_small_move_trials ) )
info_file_details.append( "Number of SugarSmallMoves per trial:\t%s\n" %str( input_args.num_moves_per_trial ) )
info_file_details.append( "Random reset of Fc glycan?:\t\t%s\n" %str( input_args.random_reset ) )
info_file_details.append( "Using score ramping?:\t\t\t%s\n" %str( input_args.ramp_sf ) )
info_file_details.append( "Native constraint file used?:\t\t%s\n" %str( input_args.native_constraint_file ).split( '/' )[-1] )
info_file_details.append( "ScoreFunction file used?:\t\t%s\n" %str( input_args.scorefxn_file ).split( '/' )[-1] )
info_file_details.append( "Angle multiplier used:\t\t\t%s\n" %str( input_args.angle_multiplier ) )
info_file_details.append( "Main structure directory:\t\t%s\n" %main_structure_dir )
info_file_details.append( "Base structure directory:\t\t%s\n" %base_structs_dir )
info_file_details.append( "Lowest E structure directory:\t\t%s\n" %lowest_E_structs_dir )
if input_args.ramp_sf:
    score_types = main_sf.get_nonzero_weighted_scoretypes()
    info_file_details.append( "\nScore weights used in main_sf:\n" )
    for score_type in score_types:
        if str( score_type ) == "fa_atr": 
            info_file_details.append( "%s: %s * 2 in ramp\n" %( str( score_type ), str( main_sf.get_weight( score_type ) ) ) )
        elif str( score_type ) == "fa_rep": 
            info_file_details.append( "%s: %s * 0.5 in ramp\n" %( str( score_type ), str( main_sf.get_weight( score_type ) ) ) )
        else:
            info_file_details.append( "%s: %s\n" %( str( score_type ), str( main_sf.get_weight( score_type ) ) ) )    
else:
    info_file_details.append( "\nScore weights used in main_sf:\n%s\n" %( "\n".join( [ "%s: %s" %( str( name ), main_sf.get_weight( name ) ) for name in main_sf.get_nonzero_weighted_scoretypes() ] ) ) )
info_file = ''.join( info_file_details )
print "\n", info_file, "\n"

# write out the info file with the collected info from above
info_filename = main_structure_dir + "protocol_run.info"
with open( info_filename, "wb" ) as fh:
    fh.write( "Info for this run of %s\n\n" %__file__ )
    fh.write( info_file )



#########################
#### JOB DISTRIBUTOR ####
#########################

# create and use the PyJobDistributor object
jd = PyJobDistributor( working_pose_decoy_name, input_args.nstruct, main_sf )
jd.native_pose = native_pose
cur_decoy_num = 1

print "Running SugarSmallMover PyJobDistributor..."

while not jd.job_complete:
    # get a fresh copy of the working pose to be used in this protocol
    testing_pose = Pose()
    testing_pose.assign( native_pose )
    testing_pose.pdb_info().name( "decoy_num_" + str( cur_decoy_num ) )
    pmm.apply( testing_pose )


    #####################################
    #### ATTACH INFO TO TESTING POSE ####
    #####################################

    # instantiate the 3ay4 information holder class object
    testing_pose_info = hold_chain_and_res_designations_3ay4()

    # see antibody_functions for more information on this hard-coded function
    testing_pose_info.native()



    #########################
    #### ADD CONSTRAINTS ####
    #########################

    # apply the constraints to the pose if a .cst file was passed
    if constraint_file_used:
        try:
            constraint_setter.apply( testing_pose )
        except:
            print "It appears there is something wrong with your constraint file. Please check your input."
            sys.exit()



    ################################
    #### Fc GLYCAN RANDOM RESET ####
    ################################

    if input_args.random_reset:
        # this is for each residue of the Fc glycan except the core GlcNAc
        # changing omega on residues that don't have an omega torsion doesn't affect them
        for res_num in testing_pose_info.native_Fc_glycan_nums_except_core_GlcNAc:
            # pick three random integers between 0 and 360 - ie. -180 to 180
            reset_phi_num = float( random_range( 0, 360 ) )
            reset_psi_num = float( random_range( 0, 360 ) )
            reset_omega_num = float( random_range( 0, 360 ) )

            # reset the phi, psi, and omega values for the residue
            testing_pose.set_phi( res_num, reset_phi_num )
            testing_pose.set_psi( res_num, reset_psi_num )
            testing_pose.set_omega( res_num, reset_omega_num )

        pmm.apply( testing_pose )
        if input_args.verbose:
            print "score of random reset:", main_sf( testing_pose )



    #################################
    #### Fc GLYCAN AREA PACK/MIN ####
    #################################

    # make backbone and chi MoveMap for the Fc sugars
    min_mm = MoveMap()
    for res_num in testing_pose_info.native_Fc_glycan_nums: 
        min_mm.set_bb( res_num, True )
    for res_num in testing_pose_info.native_Fc_glycan_nums: 
        min_mm.set_chi( res_num, True )

    # add in the Fc glycan branch points discluding the ASN connection
    for branch_point in testing_pose_info.native_Fc_glycan_branch_point_nums:
        min_mm.set_branches( branch_point, True )

    # make and apply MinMover
    Fc_glycan_min_mover = MinMover( movemap_in = min_mm, 
                                    scorefxn_in = main_sf, 
                                    min_type_in = "dfpmin_strong_wolfe", 
                                    tolerance_in = 0.01, 
                                    use_nb_list_in = True )

    # pack the Fc sugars and around them within 20 Angstroms
    pack_rotamers_mover = make_pack_rotamers_mover( main_sf, testing_pose, 
                                                    apply_sf_sugar_constraints = False,
                                                    pack_branch_points = True, 
                                                    residue_range = testing_pose_info.native_Fc_glycan_nums, 
                                                    use_pack_radius = True, 
                                                    pack_radius = 20 )

    # do 2 pack/mins to try to get to a low-energy structure of and around the Fc glycan
    # Poses that don't get to a negative score here will likely be outliers
    for ii in range( 1, 2 + 1 ):
        # pack
        pack_rotamers_mover.apply( testing_pose )
        if input_args.verbose:
            print "score of pre-pack:", ii, main_sf( testing_pose )

        # minimize
        Fc_glycan_min_mover.apply( testing_pose )
        if input_args.verbose:
            print "score of pre-min:", ii, main_sf( testing_pose )
        
        pmm.apply( testing_pose )



    #########################
    #### SugarSmallMover ####
    #########################

    ## use the SugarSmallMover find a local sugar minima
    # make an appropriate MonteCarlo object
    mc = MonteCarlo( testing_pose, main_sf, kT )

    # make a SmallMover as to get the angle_max value
    sm = SmallMover()
    angle_max = sm.get_angle_max( 'L' ) * input_args.angle_multiplier

    # raise the fa_atr term and lower the fa_rep term in the ScoreFunction for ramping, if desired
    if input_args.ramp_sf:
        # store the original fa_atr, fa_rep, and fa_elec weights
        FA_ATR_ORIG = main_sf.get_weight( score_type_from_name( "fa_atr" ) )
        FA_REP_ORIG = main_sf.get_weight( score_type_from_name( "fa_rep" ) )

        # double the fa_atr weight
        FA_ATR_NEW = FA_ATR_ORIG * 2
        main_sf.set_weight( score_type_from_name( "fa_atr" ), FA_ATR_NEW )

        # half the fa_rep weight
        FA_REP_NEW = FA_REP_ORIG * 0.5
        main_sf.set_weight( score_type_from_name( "fa_rep" ), FA_REP_NEW )

    # run the SugarSmallMover a range of times using a MonteCarlo object to accept or reject the move
    num_ssm_accept = 0
    num_mc_checks = 0
    for ii in range( 1, input_args.num_sugar_small_move_trials + 1 ):
        # if score ramping is desired
        if input_args.ramp_sf:
            # ramp up or down the appropriate scoring terms and get it back to the MonteCarlo object
            main_sf = ramp_score_weight( main_sf, 
                                         "fa_atr", 
                                         FA_ATR_ORIG, 
                                         ii - 1, 
                                         input_args.num_sugar_small_move_trials )
            main_sf = ramp_score_weight( main_sf, 
                                         "fa_rep", 
                                         FA_REP_ORIG, 
                                         ii - 1, 
                                         input_args.num_sugar_small_move_trials )
            mc.score_function( main_sf )

        # print current score
        if input_args.verbose:
            print "\nstarting score:", main_sf( testing_pose )

        # for as many moves per trial as desired
        for jj in range( input_args.num_moves_per_trial ):
            # pick a random Fc glycan residue except the core GlcNAc
            res_num = random.choice( testing_pose_info.native_Fc_glycan_nums_except_core_GlcNAc )

            # apply the SugarSmallMover and change phi, psi, and omega
            testing_pose.assign( SugarSmallMover( res_num, testing_pose, angle_max, 
                                                  set_phi = True, 
                                                  set_psi = True, 
                                                  set_omega = True ) )
        if input_args.verbose:
            print "score after SugarSmallMover:", main_sf( testing_pose )

        # pack the Fc sugars and around them within 20 Angstroms every other trial
        # trials run 1 through num_trials + 1, so pack on every odd trial
        if ii % 2 != 0:
            # use previously-made pack_rotamers_mover
            pack_rotamers_mover.apply( testing_pose )
            if input_args.verbose:
                print "score after pack:", main_sf( testing_pose )

        # minizmize the sugars using the previously-made Fc_glycan_min_mover
        Fc_glycan_min_mover.apply( testing_pose )
        if input_args.verbose:
            print "score after min:", main_sf( testing_pose )

        # accept or reject the total move using the MonteCarlo object
        if mc.boltzmann( testing_pose ):
            num_ssm_accept += 1
            pmm.apply( testing_pose )
        num_mc_checks += 1

        # print out the MC acceptance rate every 3 trials
        mc_acceptance = round( ( float( num_ssm_accept ) / float( num_mc_checks ) * 100 ), 2 )
        if ii % 3 == 0:
            if input_args.verbose:
                print "Moves made so far:", num_mc_checks, 
                print "  Moves accepted:", num_ssm_accept, 
                print "  Acceptance rate:", mc_acceptance

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


# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = working_pose_decoy_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 5 )
