#!/usr/bin/python
__author__ = "morganlnance"

'''
NOTES
STARTING POSE (3ay4 without Fc) is a base structure that was acquired from making 1000 pack/min decoys of native 3ay4 (with Fc sugar on) and then manually removing the Fc sugars of the lowest energy decoy
'''

'''
SAMPLE INPUT
run glycan_sampling_just_GRM.py pdb_copies_dont_touch/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb pdb_copies_dont_touch/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII_removed_Fc_sugar.pdb /Users/Research/antibody_project/send_to_louis/project_glyco_files/3ay4_Fc_Glycan.iupac /Users/Research/antibody_project/send_to_louis/project_utility_files/ /Users/Research/pyrosetta_dir/test_pdb_dir/ 1 10
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
parser.add_argument("working_pdb_file", type=str, help="the filename of the PDB structure to be glycosylated.")
parser.add_argument("glyco_file", type=str, help="/path/to/the .iupac glycan file to be used.")
parser.add_argument("utility_dir", type=str, help="where do your utility files live? Give me the directory.")
parser.add_argument("structure_dir", type=str, help="where do you want to dump the decoys made during this protocol?")
parser.add_argument("nstruct", type=int, help="how many decoys do you want to make using this protocol?")
parser.add_argument("num_GRM_moves", type=int, help="how many moves do you want to make within the Fc glycan during one GlycanRelaxMover trial?")
parser.add_argument("--random_reset", action="store_true", help="do you want to randomly reset the phi and psi values of the Fc glycan? (Excluding core GlcNAc)")
parser.add_argument("--verbose", "-v", action="store_true", default=False, help="do you want the program to print out pose scores during the protocol?")
input_args = parser.parse_args()



##########################
#### CHECK ALL INPUTS ####
##########################

import os, sys

# collect and create necessary directories for use in metric calculations
working_dir = os.getcwd() + '/'
metrics_dump_dir = working_dir + "GRM_%s_dir" %str( input_args.num_GRM_moves )
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

# check the utility directory
if not os.path.isdir( input_args.utility_dir ):
    print "Your argument", input_args.utility_dir, "for utility_directory is not a directory, exiting"
    sys.exit()

# add the utility directory to the system path for loading of modules
sys.path.append( input_args.utility_dir )

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

# relay information to user
info_file_details = []
info_file_details.append( "Native PDB filename:\t\t%s\n" %input_args.native_pdb_file.split( '/' )[-1] )
info_file_details.append( "Working PDB filename:\t\t%s\n" %input_args.working_pdb_file.split( '/' )[-1] )
info_file_details.append( "Sugar filename:\t\t\t%s\n" %input_args.glyco_file.split( '/' )[-1] )
info_file_details.append( "Creating this many decoys:\t%s\n" %str( input_args.nstruct ) )
info_file_details.append( "Number of GRM moves:\t\t%s\n" %str( input_args.num_GRM_moves ) )
info_file_details.append( "Random reset of Fc glycan?:\t%s\n" %str( input_args.random_reset ) )
info_file_details.append( "Main structure directory:\t%s\n" %main_structure_dir )
info_file_details.append( "Base structure directory:\t%s\n" %base_structs_dir )
info_file_details.append( "Lowest E structure directory:\t%s\n" %lowest_E_structs_dir )
info_file = ''.join( info_file_details )
print "\n", info_file, "\n"

# write out the info file with the collected info from above
info_filename = main_structure_dir + "protocol_run.info"
with open( info_filename, "wb" ) as fh:
    fh.write( "Info for this run of %s\n\n" %__file__ )
    fh.write( info_file )



#################
#### IMPORTS ####
#################

# Rosetta worker functions
from rosetta import Pose, pose_from_file, get_fa_scorefxn, \
    PyMOL_Mover, MonteCarlo, PyJobDistributor
from rosetta.numeric.random import random_range, uniform
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file
from rosetta.core.scoring import score_type_from_name
from rosetta.protocols.carbohydrates import GlycanRelaxMover
from rosetta import MinMover

# Rosetta functions I wrote out
from antibody_functions import initialize_rosetta, \
    get_fa_scorefxn_with_given_weights, make_pack_rotamers_mover, \
    make_movemap_for_range, load_pose, get_phi_psi_omega_of_res, \
    ramp_score_weight

# utility functions
from file_mover_based_on_fasc import main as get_lowest_E_from_fasc
from pose_metrics_util import get_pose_metrics



################################
#### INITIAL PROTOCOL SETUP ####
################################

# initialize Rosetta
initialize_rosetta()

## load up the poses given from the arguments passed
# native pose ( for comparison, really )
native_pose = Pose()
native_pose.assign( load_pose( input_args.native_pdb_file ) )

# get the full path of the native PDB name
native_pdb_filename_full_path = input_args.native_pdb_file
native_pdb_filename = native_pdb_filename_full_path.split( '/' )[-1]
native_pdb_name = native_pdb_filename.split( ".pdb" )[0]

# change the name of the native PDB name
native_pose_name = "native_pose"
native_pose.pdb_info().name( native_pose_name )


# load up the working pose
working_pose = Pose()
working_pose.assign( load_pose( input_args.working_pdb_file ) )

# get the full path of the working pose PDB name
working_pdb_filename_full_path = input_args.working_pdb_file
working_pdb_filename = working_pdb_filename_full_path.split( '/' )[-1]
working_pdb_name = working_pdb_filename.split( ".pdb" )[0]

# change the name of the working PDB name
working_pose_name = "working_pose"
working_pose.pdb_info().name( working_pose_name )

# make the directory for the working PDBs in the base_structs_dir
structure_dir = base_structs_dir + working_pdb_name
if not os.path.isdir( structure_dir ):
    try:
        os.mkdir( structure_dir )
    except:
        pass
working_pose_decoy_name = structure_dir + '/' + working_pdb_name + "_glycosylated_then_just_%s_grm" %input_args.num_GRM_moves


# collect the core GlcNAc values from the native pose
A_phi, A_psi, A_omega = get_phi_psi_omega_of_res( native_pose, 216 )
B_phi, B_psi, B_omega = get_phi_psi_omega_of_res( native_pose, 440 )

# collect the chain id's of the Fc glycan from the native pose
native_Fc_glycan_chains = [ 'D', 'E', 'F', 'G' ]

## get some numbers that will be used in pieces of this protocol
# this number is used later for resetting the core glycan
n_res_no_Fc_glycan = working_pose.n_residue()

# these numbers are of just the receptor glycan - they should be ignored sometimes
# and get the residue numbers of branch points - they should also be ignored sometimes
FcR_seqpos_nums = []
FcR_branch_point_nums = []
working_pose_chains = []
for res in working_pose:
    # residue numbers
    if res.is_carbohydrate():
        FcR_seqpos_nums.append( res.seqpos() )
        
    # branch points
    if res.is_branch_point():
        FcR_branch_point_nums.append( res.seqpos() )
        
    # chain id's
    res_chain = working_pose.pdb_info().chain( res.seqpos() )
    if res_chain not in working_pose_chains:
        working_pose_chains.append( res_chain )

# get the number of chains because the Pose renumbers its chains after glycosylation
num_working_pose_chains = len( working_pose_chains )

# get a standard fa_scorefxn for protein stuff
sf = get_fa_scorefxn()

# adjust the standard fa_scorefxn for sugar stuff
sugar_sf = get_fa_scorefxn_with_given_weights( "fa_intra_rep", 0.440 )



# pymol stuff
pmm = PyMOL_Mover()
pmm.keep_history( True )
pmm.apply( native_pose )
pmm.apply( working_pose )



#########################
#### JOB DISTRIBUTOR ####
#########################

# create and use the PyJobDistributor object
jd = PyJobDistributor( working_pose_decoy_name, input_args.nstruct, sugar_sf )
jd.native_pose = native_pose
cur_decoy_num = 0

print "Running GRM PyJobDistributor..."

while not jd.job_complete:
    # get a fresh copy of the working pose to be used in this protocol
    testing_pose = Pose()
    testing_pose.assign( working_pose )



    ##########################
    #### GLYCOSYLATE POSE #### 
    ##########################

    # glycosylate the given testing_pose
    # 69 and 284 are the two ASN297 residues from 3ay4 ( pose numbering system, not PDB )
    glycosylate_these_ASN = [ 69, 284 ]
    for ASN in glycosylate_these_ASN:
        glycosylate_pose_by_file( testing_pose, ASN, "ND2", input_args.glyco_file )

    testing_pose.pdb_info().name( "decoy_num_" + str( cur_decoy_num ) )
    pmm.apply( testing_pose )
    if input_args.verbose:
        print
        print "score of glycosylated pose", sf( testing_pose )



    ###################################
    #### COLLECT Fc GLYCAN NUMBERS ####
    ###################################

    # pull out some core sugar residue numbers
    n_res_Fc_glycan = testing_pose.n_residue()
    num_sugars_added = n_res_Fc_glycan - n_res_no_Fc_glycan
    size_of_one_glycan = num_sugars_added / 2
    A_core_GlcNAc = n_res_no_Fc_glycan + 1
    B_core_GlcNAc = n_res_no_Fc_glycan + size_of_one_glycan + 1    

    # get the res nums and branch points of the Fc sugars added
    Fc_sugar_nums = []
    Fc_sugar_nums_except_core_GlcNAc = []
    Fc_glycan_branch_point_nums = []

    for res in testing_pose:
        # if the residue is a carbohydrate
        if res.is_carbohydrate():
            # if the residue number is not in the FcR
            if res.seqpos() not in FcR_seqpos_nums:
                Fc_sugar_nums.append( res.seqpos() )

            # if the residue is a branch point
            if res.is_branch_point():
                # if it's not a branch point found in the FcR
                if res.seqpos() not in FcR_branch_point_nums:
                    # if it's a branch point that is a sugar ( ie. not the linking ASN )
                    if res.is_carbohydrate():
                        Fc_glycan_branch_point_nums.append( res.seqpos() )

    # get a list of the Fc sugars discluding the core GlcNAc residues
    for res_num in Fc_sugar_nums:
        if res_num != A_core_GlcNAc and res_num != B_core_GlcNAc:
            Fc_sugar_nums_except_core_GlcNAc.append( res_num )

    # get the chain id's of the glycosylated testing_pose
    testing_pose_chains = []
    for res in testing_pose:
        # chain id's
        res_chain = testing_pose.pdb_info().chain( res.seqpos() )
        if res_chain not in testing_pose_chains:
            testing_pose_chains.append( res_chain )
    num_testing_pose_chains = len( testing_pose_chains )

    # pull out the chain id's based on how many chains were added to testing_pose
    # Rosetta seems to just rename the chains, so the last X chains added should be the glycan
    num_chains_added = num_testing_pose_chains - num_working_pose_chains
    Fc_glycan_chains = testing_pose_chains[ ( -1 * num_chains_added ) : ]



    ###########################
    #### CORE GlcNAc RESET ####
    ###########################

    ## reset both of the core GlcNAc residue of the glycosylated testing_pose
    # chain A
    testing_pose.set_phi( A_core_GlcNAc, A_phi )
    testing_pose.set_psi( A_core_GlcNAc, A_psi )
    testing_pose.set_omega( A_core_GlcNAc, A_omega )

    # chain B
    testing_pose.set_phi( B_core_GlcNAc, B_phi )
    testing_pose.set_psi( B_core_GlcNAc, B_psi )
    testing_pose.set_omega( B_core_GlcNAc, B_omega )

    pmm.apply( testing_pose )
    if input_args.verbose:
        print "score of glyco reset", sf( testing_pose )



    ################################
    #### Fc GLYCAN RANDOM RESET ####
    ################################

    if input_args.random_reset:
        # for each residue except the core GlcNAc in Fc glycan
        # this method gives the range for the first Fc glycan without core GlcNAc
        # this will only work if the two glycans are the same on both sides
        for res_num in Fc_sugar_nums_except_core_GlcNAc[ : ( size_of_one_glycan - 1 ) ]:
            # pick three random integers between 0 and 360
            reset_phi_num = float( random_range( 0, 360 ) )
            reset_psi_num = float( random_range( 0, 360 ) )
            reset_omega_num = float( random_range( 0, 360 ) )

            # reset the phi, psi, and omega values for the residue
            # side A
            testing_pose.set_phi( res_num, reset_phi_num )
            testing_pose.set_psi( res_num, reset_psi_num )
            #testing_pose.set_omega( res_num, reset_omega_num )

            # side B
            testing_pose.set_phi( res_num + size_of_one_glycan - 1, reset_phi_num )
            testing_pose.set_psi( res_num + size_of_one_glycan - 1, reset_phi_num )
            #testing_pose.set_omega( res_num + size_of_one_glycan - 1, reset_phi_num )

        pmm.apply( testing_pose )
        if input_args.verbose:
            print "score of random reset", sf( testing_pose )



    #################################
    #### Fc GLYCAN AREA PACK/MIN ####
    #################################

    # make MoveMap for the Fc sugars
    min_mm = make_movemap_for_range( Fc_sugar_nums, 
                                     allow_bb_movement = True, 
                                     allow_chi_movement = True )

    # add in the branch points myself ( does not include the two ASN residues )
    for branch_point in Fc_glycan_branch_point_nums:
        min_mm.set_branches( branch_point, True )

    # make and apply MinMover
    min_mover = MinMover( movemap_in = min_mm, 
                          scorefxn_in = sugar_sf, 
                          min_type_in = "dfpmin_strong_wolfe", 
                          tolerance_in = 0.01, 
                          use_nb_list_in = True )

    # pack the Fc sugars and around them within 20 Angstroms
    pack_rotamers_mover = make_pack_rotamers_mover( sugar_sf, testing_pose, 
                                                    apply_sf_sugar_constraints = False,
                                                    pack_branch_points = True, 
                                                    residue_range = Fc_sugar_nums, 
                                                    use_pack_radius = True, 
                                                    pack_radius = 20 )

    # do 2 pack/mins to try to get to a low-energy structure
    # those that don't get to a negative score will probably just be outliers?
    for ii in range( 2 ):
        # pack
        pack_rotamers_mover.apply( testing_pose )
        if input_args.verbose:
            print "score of pack", sf( testing_pose )

        # minimize
        min_mover.apply( testing_pose )
        if input_args.verbose:
            print "score of min", sf( testing_pose )
        
        pmm.apply( testing_pose )



    ##########################
    #### GlycanRelaxMover ####
    ##########################

    ## use the GlycanRelaxMover find a local sugar minima
    # make a MoveMap for these Fc sugars allowing only bb movement
    grm_mm = make_movemap_for_range( Fc_sugar_nums_except_core_GlcNAc, 
                                     allow_bb_movement = True, 
                                     allow_chi_movement = False )

    # add in the branch points myself ( does not include the two ASN residues )
    for branch_point in Fc_glycan_branch_point_nums:
        grm_mm.set_branches( branch_point, True )

    # make an appropriate GlycanRelaxMover
    grm = GlycanRelaxMover( mm = grm_mm, scorefxn = sugar_sf, rounds = input_args.num_GRM_moves )

    # make an appropriate MonteCarlo object
    mc = MonteCarlo( testing_pose, sugar_sf, kT )

    # apply the GlycanRelaxMover
    grm.apply( testing_pose )

    # pack the Fc sugars and around them within 20 Angstroms
    pack_rotamers_mover = make_pack_rotamers_mover( sugar_sf, testing_pose, 
                                                    apply_sf_sugar_constraints = False,
                                                    pack_branch_points = True, 
                                                    residue_range = Fc_sugar_nums, 
                                                    use_pack_radius = True, 
                                                    pack_radius = 20 )
    pack_rotamers_mover.apply( testing_pose )
    if input_args.verbose:
        print "score after pack", sf( testing_pose )

    # make MoveMap for the Fc sugars
    min_mm = make_movemap_for_range( Fc_sugar_nums, 
                                     allow_bb_movement = True, 
                                     allow_chi_movement = True )

    # add in the branch points myself ( does not include the two ASN residues )
    for branch_point in Fc_glycan_branch_point_nums:
        min_mm.set_branches( branch_point, True )

    # make and apply MinMover
    min_mover = MinMover( movemap_in = min_mm, 
                          scorefxn_in = sugar_sf, 
                          min_type_in = "dfpmin_strong_wolfe", 
                          tolerance_in = 0.01, 
                          use_nb_list_in = True )
    min_mover.apply( testing_pose )
    if input_args.verbose:
        print "score after min", sf( testing_pose )

    # collect additional metric data
    try:
        metrics = get_pose_metrics( testing_pose, 
                                    native_pose, 
                                    sugar_sf, 
                                    2, # interface JUMP_NUM
                                    Fc_glycan_chains, 
                                    native_Fc_glycan_chains, 
                                    jd.current_num, 
                                    metrics_dump_dir )
    except:
        metrics = ''
        pass

    # add the metric data to the .fasc file
    jd.additional_decoy_info = metrics

    # dump the decoy
    jd.output_decoy( testing_pose )
    cur_decoy_num += 1


# remove the metrics_dump_dir
try:
    os.rmdir( metrics_dump_dir )
except:
    pass

# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = working_pose_decoy_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 5 )
