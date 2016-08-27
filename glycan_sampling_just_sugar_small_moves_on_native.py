#!/usr/bin/python
__author__ = "morganlnance"

'''
NOTES
STARTING POSE (3ay4 without Fc) is a base structure that was acquired from making 1000 pack/min decoys of native 3ay4 (with Fc sugar on) and then manually removing the Fc sugars of the lowest energy decoy
'''

'''
SAMPLE INPUT
run glycan_sampling_just_sugar_small_moves_on_native.py project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb project_glyco_files/3ay4_Fc_Glycan.iupac project_utility_files/ ~/pyrosetta_dir/test_pdb_dir/ 1 3 3
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
#parser.add_argument("--just_chain_A", action="store_true", help="do you want the SugarSmallMover (and any reset) to act on the chain A glycan only?")
parser.add_argument("--LCM_reset", action="store_true", help="do you want a LinkageConformerMover to reset the phi, psi, and omega values of the Fc glycan? (Excluding core GlcNAc)")
parser.add_argument("--use_population_ideal_LCM_reset", action="store_true", help="do you want the LinkageConformerMover to reset the phi, psi, and omega values of the Fc glycan to ideal values (stdev of 0) using population weights? (Excluding core GlcNAc)")
parser.add_argument("--use_ideal_LCM_reset", action="store_true", help="do you want the LinkageConformerMover to reset the phi, psi, and omega values of the Fc glycan to ideal values (stdev of 0) of only the highest population weight? (Excluding core GlcNAc)")
parser.add_argument("--light_reset", action="store_true", help="do you want the random reset to be light? +/- 10-15 degrees discluding 0")
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
if input_args.glyco_file is not None:
    if not os.path.isfile( input_args.glyco_file ):
        print "\nYour glyco_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.glyco_file
        sys.exit()
if input_args.scorefxn_file is not None:
    if not os.path.isfile( input_args.scorefxn_file ):
        print "\nYour scorefxn_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.scorefxn_file
        sys.exit()
if input_args.native_constraint_file is not None:
    if not os.path.isfile( input_args.native_constraint_file ):
        print "\nYour native_constraint_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.native_constraint_file
        sys.exit()


## double check the logic of the passed arguments
# make sure that only one kind of reset argument was set to True
if input_args.light_reset is True and input_args.LCM_reset is True:
    print "\nYou asked for both a light reset and an LCM reset. Please only specify one type of reset. Exiting."
    sys.exit()

# if a different LCM reset option was selected, make sure only one was and that LCM_reset was set to True
if input_args.use_ideal_LCM_reset is True and input_args.use_population_ideal_LCM_reset is True:
    print "\nYou asked for both an ideal LCM reset and an ideal LCM reset using population weights. Please only specify one type of LCM reset option. Exiting"
    sys.exit()
if input_args.use_ideal_LCM_reset is True and input_args.LCM_reset is not True:
    print "You asked for a specific type of LCM reset ( use_ideal_LCM_reset ), but did not tell me you actually want an LCM reset. Please specify the --LCM_reset option. Exiting for safety."
    sys.exit()
if input_args.use_population_ideal_LCM_reset is True and input_args.LCM_reset is not True:
    print "You asked for a specific type of LCM reset ( use_population_ideal_LCM_reset ), but did not tell me you actually want an LCM reset. Please specify the --LCM_reset option. Exiting for safety."
    sys.exit()

# check that num_sugar_small_move_trials is 10 or more if ramp_sf is set to True
if input_args.ramp_sf == True and input_args.num_sugar_small_move_trials < 10:
    print "\nIf you are ramping the ScoreFunction then you need to do 10 or more sugar_small_move_trials for the math to work. Exiting."
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

# collect and create necessary directories for use in metric calculations
working_dir = os.getcwd() + '/'
metrics_dump_dir = working_dir + "sugar_small_dir_%s" %str( input_args.num_sugar_small_move_trials )
try:
    os.mkdir( metrics_dump_dir )
except:
    pass



#################
#### IMPORTS ####
#################

# Rosetta worker functions
from rosetta import Pose, pose_from_file, get_fa_scorefxn, \
    PyMOL_Mover, MonteCarlo, PyJobDistributor, MoveMap
from rosetta.numeric.random import random_range, uniform
from rosetta.core.id import phi_dihedral, psi_dihedral, omega_dihedral
from rosetta.core.scoring import score_type_from_name
from rosetta.core.pose.carbohydrates import set_glycosidic_torsion
from rosetta.protocols.simple_moves import ConstraintSetMover
from rosetta.protocols.carbohydrates import LinkageConformerMover
from rosetta import SmallMover, MinMover
from toolbox import get_hbonds

# Rosetta functions I wrote out
from antibody_functions import initialize_rosetta, \
    get_fa_scorefxn_with_given_weights, make_pack_rotamers_mover, \
    load_pose, ramp_score_weight, make_fa_scorefxn_from_file, \
    SugarSmallMover, hold_chain_and_res_designations_3ay4, \
    set_3ay4_Fc_glycan_except_core_GlcNAc_to_ideal_LCM_phi_psi_omega, \
    get_res_nums_within_radius_of_residue_list

# utility functions
from file_mover_based_on_fasc import main as get_lowest_E_from_fasc
from get_pose_metrics_on_native import main as get_pose_metrics_on_native
from pose_metrics_util import Fc_glycan_rmsd



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
#info_file_details.append( "Just move (and reset) chain A?:\t%s\n" %str( input_args.just_chain_A ) )
info_file_details.append( "LCM reset of Fc glycan?:\t\t%s\n" %str( input_args.LCM_reset ) )
info_file_details.append( "Use main ideal in LCM reset?:\t\t%s\n" %str( input_args.use_ideal_LCM_reset ) )
info_file_details.append( "Use population ideals in LCM reset?:\t%s\n" %str( input_args.use_population_ideal_LCM_reset ) )
info_file_details.append( "Light reset of Fc glycan?:\t\t%s\n" %str( input_args.light_reset ) )
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
    #### Fc GLYCAN RESET ####
    #########################

    # if user wants an LCM reset
    if input_args.LCM_reset:
        # for each residue except core GlcNAc
        for res_num in testing_pose_info.native_Fc_glycan_nums_except_core_GlcNAc:
            # make a MoveMap for this single residue
            res_mm = MoveMap()

            # set bb to True ( phi, psi, omega )
            res_mm.set_bb( res_num, True )

            # make an appropriate LinkageConformerMover
            lcm = LinkageConformerMover()
            lcm.set_movemap( res_mm )

            # if the user only wants to use ideals, but within the different population clusters
            if input_args.use_population_ideal_LCM_reset:
                # set_idealize_torsions uses ideal values instead of sampling from stdev
                lcm.set_idealize_torsions( True )
                # use_conformer_population_stats is True by default, but setting for clarity
                lcm.set_use_conformer_population_stats( True )

                # apply the LCM
                lcm.apply( testing_pose )

            # elif the user wants to use ideals of only the highest population statistics
            elif input_args.use_ideal_LCM_reset:
                # this is hardcoded data at the moment
                # setting residue 3 on chain D and F (the Man with the branch) to the native phi for now to see if that would help get better decoys
                testing_pose = set_3ay4_Fc_glycan_except_core_GlcNAc_to_ideal_LCM_phi_psi_omega( testing_pose, set_3_D_and_F_phi_to_native = True )

                # this option of LCM reset doesn't actually use the LCM, just ideal data from it
                # so there is no actual call to lcm.apply()

            # else, standard LCM reset using (by default) population data and a stdev of 1
            else:
                # setting these options for clarity
                lcm.set_use_conformer_population_stats( True )
                lcm.set_x_standard_deviations( 1 )

                # apply the LCM
                lcm.apply( testing_pose )

        pmm.apply( testing_pose )
        if input_args.verbose:
            print "score of LCM reset:", main_sf( testing_pose )

    # if user wants a light reset by perturbing starting glycan phi/psi/omega
    if input_args.light_reset:
        # for each residue except core GlcNAc
        for res_num in testing_pose_info.native_Fc_glycan_nums_except_core_GlcNAc:

            # determine if the phi, psi, and omega change will be positive or negative
            if random_range( 0, 1 ) == 1:
                phi_mult = 1
            else:
                phi_mult = -1
            if random_range( 0, 1 ) == 1:
                psi_mult = 1
            else:
                psi_mult = -1
            if random_range( 0, 1 ) == 1:
                omega_mult = 1
            else:
                omega_mult = -1

            # create the new phi, psi, and omega
            reset_phi_num = testing_pose.phi( res_num ) + ( random_range( 10, 15 ) * phi_mult )
            reset_psi_num = testing_pose.psi( res_num ) + ( random_range( 10, 15 ) * psi_mult )
            reset_omega_num = testing_pose.omega( res_num ) + ( random_range( 10, 15 ) * omega_mult )

            # reset the phi, psi, and omega values for the residue
            # using this instead of set_phi (etc) because this explicitly calls align_virtual_atoms_in_carbohydrate_residue
            set_glycosidic_torsion( phi_dihedral, testing_pose, res_num, reset_phi_num )
            set_glycosidic_torsion( psi_dihedral, testing_pose, res_num, reset_psi_num )
            set_glycosidic_torsion( omega_dihedral, testing_pose, res_num, reset_omega_num )

        pmm.apply( testing_pose )                                                      
        if input_args.verbose:
            print "score of light reset:", main_sf( testing_pose )



    #############################
    #### Fc GLYCAN RMSD CALC ####
    #############################

    # find out how much the reset (if any) moved the glycan
    '''
    pre_SSM_metrics = []
    if input_args.LCM_reset or input_args.light_reset:
        Fc_glycan_rmsd_after_reset = Fc_glycan_rmsd( testing_pose, testing_pose_info.native_Fc_glycan_chains, 
                                                     native_pose, native_pose_info.native_Fc_glycan_chains, 
                                                     cur_decoy_num, metrics_dump_dir )
    
        pre_SSM_metrics.append( "Fc_glycan_rmsd_after_reset:" )
        pre_SSM_metrics.append( str( Fc_glycan_rmsd_after_reset ) )
    
        if input_args.verbose:
            print "Fc glycan RMSD after reset:", Fc_glycan_rmsd_after_reset
    '''



    #################################
    #### Fc GLYCAN AREA PACK/MIN ####
    #################################

    '''
    # make backbone and chi MoveMap for the Fc sugars
    min_mm = MoveMap()
    for res_num in testing_pose_info.native_Fc_glycan_nums: 
        min_mm.set_bb( res_num, True )
        min_mm.set_chi( res_num, True )

    # add in the Fc glycan branch points discluding the ASN connection
    for branch_point in testing_pose_info.native_Fc_glycan_branch_point_nums:
        min_mm.set_branches( branch_point, True )

    # pack the Fc sugars and around them within 20 Angstroms
    pack_rotamers_mover = make_pack_rotamers_mover( main_sf, testing_pose, 
                                                    apply_sf_sugar_constraints = False,
                                                    pack_branch_points = True, 
                                                    residue_range = testing_pose_info.native_Fc_glycan_nums, 
                                                    use_pack_radius = True, 
                                                    pack_radius = 20 )
    '''
    # make a MoveMap
    residue_range = get_res_nums_within_radius_of_residue_list( testing_pose_info.native_Fc_glycan_nums, testing_pose, 10, 
                                                                include_res_nums = False )
    min_mm = MoveMap()
    for res_num in residue_range:
        min_mm.set_bb( res_num, True )
        min_mm.set_chi( res_num, True )
    for branch_point in testing_pose_info.native_Fc_glycan_branch_point_nums:
        min_mm.set_branches( branch_point, False )

    # make the MinMover
    Fc_glycan_min_mover = MinMover( movemap_in = min_mm, 
                                    scorefxn_in = main_sf, 
                                    min_type_in = "dfpmin_strong_wolfe", 
                                    tolerance_in = 0.01, 
                                    use_nb_list_in = True )

    # pack around the Fc sugars within 25 Angstroms
    pack_rotamers_mover = make_pack_rotamers_mover( main_sf, testing_pose, 
                                                    apply_sf_sugar_constraints = False,
                                                    pack_branch_points = True, 
                                                    residue_range = residue_range )

    # do 2 pack/mins to try to get to a low-energy structure of and around the Fc glycan
    # Poses that don't get to a negative score here will likely be outliers
    '''
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
    '''


    #############################
    #### Fc GLYCAN RMSD CALC ####
    #############################

    # find out how much the pack/min moved the glycan
    '''
    Fc_glycan_rmsd_after_packmin = Fc_glycan_rmsd( testing_pose, testing_pose_info.native_Fc_glycan_chains, 
                                                 native_pose, native_pose_info.native_Fc_glycan_chains, 
                                                 cur_decoy_num, metrics_dump_dir )
    
    pre_SSM_metrics.append( "Fc_glycan_rmsd_after_packmin:" )
    pre_SSM_metrics.append( str( Fc_glycan_rmsd_after_packmin ) )
    
    if input_args.verbose:
        print "Fc glycan RMSD after pack/min rounds:", Fc_glycan_rmsd_after_packmin
    '''



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
    mc_acceptance = None
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
        '''
        if ii % 2 != 0:
            # use previously-made pack_rotamers_mover
            pack_rotamers_mover.apply( testing_pose )
            if input_args.verbose:
                print "score after pack:", main_sf( testing_pose )

        # minizmize the sugars using the previously-made Fc_glycan_min_mover
        Fc_glycan_min_mover.apply( testing_pose )
        if input_args.verbose:
            print "score after min:", main_sf( testing_pose )
        '''

        # accept or reject the total move using the MonteCarlo object
        if mc.boltzmann( testing_pose ):
            num_ssm_accept += 1
            pmm.apply( testing_pose )
        num_mc_checks += 1

        # print out the MC acceptance rate every 3 trials and on the last trial
        mc_acceptance = round( ( float( num_ssm_accept ) / float( num_mc_checks ) * 100 ), 2 )
        if ii % 3 == 0 or ii == input_args.num_sugar_small_move_trials:
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

    # add the pre-metric calculations to the metrics data
    #pre_SSM_metrics = ' '.join( pre_SSM_metrics )
    #metrics += " %s" %pre_SSM_metrics

    # add the metric data to the .fasc file
    jd.additional_decoy_info = metrics

    # dump the decoy
    jd.output_decoy( testing_pose )
    cur_decoy_num += 1


# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = working_pose_decoy_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 5 )
