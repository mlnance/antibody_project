#!/usr/bin/python
__author__="morganlnance"



#########################
#### PARSE ARGUMENTS ####
#########################

import argparse

# parse and store input arguments
parser = argparse.ArgumentParser(description="Use PyRosetta to glycosylate a pose and find a low E structure of a single mutation")
parser.add_argument("native_pdb_file", type=str, help="the filename of the double pack/min native PDB structure.")
parser.add_argument("working_pdb_file", type=str, help="the filename of the working single pack/min PDB structure.")
parser.add_argument("mutation_file", type=str, help="/path/to/the file that contains all the mutations you want to make.")
parser.add_argument("utility_dir", type=str, help="where do your utility files live? Give me the directory.")
parser.add_argument("structure_dir", type=str, help="where do you want to dump the decoys made during this protocol?")
parser.add_argument("nstruct", type=int, help="how many decoys of a single mutant do you want to make?")
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
if not os.path.isfile( input_args.native_pdb_file ):
    print "\nYour native_pdb_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.native_pdb_file
    sys.exit()
if not os.path.isfile( input_args.working_pdb_file ):
    print "\nYour working_pdb_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.working_pdb_file
    sys.exit()
if not os.path.isfile( input_args.mutation_file ):
    print "\nYour mutation_file argument ( %s ) does not exist. Please check your input. Exiting" %input_args.mutation_file
    sys.exit()



##################################
#### CREATE NECESSARY OBJECTS ####
##################################

# imports 
from native_3ay4_glycan_modeling_protocol_functions import get_fa_scorefxn_with_given_weights, \
    load_pose, initialize_rosetta, native_Fc_glycan_nums_except_core_GlcNAc, \
    native_Fc_glycan_nums, make_RotamerTrialsMover
from rosetta import MoveMap, PyMOL_Mover, RotamerTrialsMover, \
    score_type_from_name, Pose, MinMover
from rosetta.core.scoring import fa_rep
# utility directory function
from file_mover_based_on_fasc import main as get_lowest_E_from_fasc


# make the appropriate directories
# create the structure_dir and lowest_E_structs dir
base_structs_dir = main_structure_dir + "base_structs/"
lowest_E_structs_dir = main_structure_dir + "lowest_E_structs/"

# make the needed dump directories
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

# create a working directory to be used for metric calculations
metrics_dump_dir = main_structure_dir + "metrics_dump_dir"
try:
    os.mkdir( metrics_dump_dir )
except:
    pass


# initialize Rosetta
initialize_rosetta()

# PyMOL_Mover
pmm = PyMOL_Mover()
pmm.keep_history( True )

# load the native lowE pose
# this should have been packed and minimized twice
low_E_native_pose = load_pose( input_args.native_pdb_file )

# load the working pose
# this should have been the first round of pack/min of the lowE pose
working_pose = load_pose( input_args.working_pdb_file )

# make a ScoreFunction
sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )
orig_fa_rep = sf.get_weight( fa_rep )


# nothing being done to the input native structure as that should have already been
# packed and minimized TWICE using the format of pack/min that is used in this
# mutational protocol. If this script had to pack/min the input native, then each
# mutant would be based on a different native, which is not what I want
# the input working_pose will have been packed and minimized twice after this
# mutational code runs. This way we are comparing apples to apples
# essentially, the low_E_native_pose serves as the lowE control and the
# working_pose serves as the experimental. both will have been packed and
# minimized the same way twice at the end of the day


###############################
#### PREPARE FOR MUTATIONS ####
###############################

# imports
from antibody_functions import read_mutation_file, get_rank_order_of_list, \
    mutate_residue, calc_interface_sasa, get_interface_score

# get the list of mutations to make
all_mutations_holder = read_mutation_file( input_args.mutation_file )
all_mutations_to_ratio_dict = all_mutations_holder.mutation_to_ratio

# create an appropriate decoy_name using the mutation_file name
mutation_filename = input_args.mutation_file.split( '/' )[-1].split( '.' )[0]
decoy_dir = base_structs_dir + mutation_filename + '/'
if not os.path.isdir( decoy_dir ):
    os.mkdir( decoy_dir )
decoy_name = decoy_dir + mutation_filename

# data holders
# again, the low_E_native_pose is the pose that has ALREADY been packed and minimized TWICE
nat_tot_score = sf( low_E_native_pose )
nat_intf_score = get_interface_score( 2, sf, low_E_native_pose )
nat_intf_sasa = calc_interface_sasa( low_E_native_pose, 2 ) 


#########################
#### JOB DISTRIBUTOR ####
#########################
from rosetta import PyJobDistributor

# create and use the PyJobDistributor object
jd = PyJobDistributor( decoy_name, input_args.nstruct, sf )
jd.native_pose = low_E_native_pose

while not jd.job_complete:
    # to hold all of the mutants to look at later
    for mutation_set, ratio in all_mutations_to_ratio_dict.items():
        # get a fresh pose object by cloning the working_pose
        # again, the working_pose is the pose that has ALREADY been packed and minimized ONCE
        mutant_pose = working_pose.clone()
        mutant_pose.pdb_info().name( main_structure_dir + mutation_set + ".pdb" )
        print "%s..." %mutation_set

    	###########################
    	#### MAKE THE MUTATION ####
    	###########################
    	# each mutation_set can include multiple point mutations, so make them one by one
        mutations = mutation_set.split( '-' )
        for mutation in mutations:
            try:
                # if this mutation requires a specific chain
                if '_' in mutation:
                    orig_residue = mutation[0]
                    pdb_num = int( mutation.split( '_' )[0][1 : -1] )
                    new_residue = mutation.split( '_' )[0][-1]
                    pdb_chain = mutation.split( '_' )[-1]
                    # mutate
                    mutant_pose.assign( mutate_residue( pdb_num, new_residue, mutant_pose, sf, 
                                                        pdb_num = True, 
                                                        pdb_chain = pdb_chain, 
                                                        do_pack = False, 
                                                        do_min = False ) )
                # otherwise, this mutation is on both chain A and B
                else:
                    orig_residue = mutation[0]
                    pdb_num = int( mutation[1 : -1] )
                    new_residue = mutation[-1]
                    # mutate chain A
                    mutant_pose.assign( mutate_residue( pdb_num, new_residue, mutant_pose, sf, 
                                                        pdb_num = True, 
                                                        pdb_chain = 'A',
                                                        do_pack = False, 
                                                        do_min = False ) )
                    # mutate chain B
                    mutant_pose.assign( mutate_residue( pdb_num, new_residue, mutant_pose, sf, 
                                                        pdb_num = True, 
                                                        pdb_chain = 'B', 
                                                        do_pack = False, 
                                                        do_min = False ) )
            except:
                # just skip mutations that don't work for whatever reason
                continue

        ###########################
    	#### GRADIENT PACK/MIN ####
    	###########################
        # collect all non-branch point residue numbers
        moveable_residues = [ res_num for res_num in range( 1, mutant_pose.n_residue() + 1 ) if mutant_pose.residue( res_num ).is_branch_point() == False ]

        # pack all residues except for branch point residues
        pack_rotamers_mover = make_RotamerTrialsMover( moveable_residues = moveable_residues,
                                                       sf = sf,
                                                       input_pose = mutant_pose,
                                                       pack_radius = None )
        pack_rotamers_mover.apply( mutant_pose )

        # we packed the side chains, so minimize them (only residues marked by moveable_residues)
        # keep the backbone the same as the crystal because 1) we can, 2) it's easier, 3) we want to see what mutations do to packing more so
        mm = MoveMap()
        mm.set_bb( False )
        for res_num in moveable_residues:
            mm.set_chi( res_num, True )
        # gradient min of mutant_pose
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
            min_mover.apply( mutant_pose )
            print "mutant_pose", sf( mutant_pose ), jj

        # collect additional metric data
        try:
            # all this is here until I update get_pose_metrics(_on_native)
            from antibody_functions import hold_chain_and_res_designations_3ay4
            from get_pose_metrics_on_native import main as get_pose_metrics_on_native
            low_E_native_pose_info = hold_chain_and_res_designations_3ay4()
            low_E_native_pose_info.native()
            working_pose_info = hold_chain_and_res_designations_3ay4()
            working_pose_info.native()
            # metric calculations
            metrics = get_pose_metrics_on_native( working_pose, 
                                                  working_pose_info, 
                                                  low_E_native_pose, 
                                                  low_E_native_pose_info, 
                                                  sf, 
                                                  2, # Fc-FcR interface JUMP_NUM
                                                  jd.current_num, 
                                                  metrics_dump_dir, 
                                                  input_args.utility_dir, 
                                                  MC_acceptance_rate = None, 
                                                  native_constraint_file = None )
        except:
            metrics = ''
            pass

        # add the metric data to the .fasc file
        jd.additional_decoy_info = metrics

        jd.output_decoy( mutant_pose )

# move the lowest E pack and minimized native structure into the lowest_E_structs dir
fasc_filename = decoy_name + ".fasc"
lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, lowest_E_structs_dir, 10 )
