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
parser.add_argument("mutation_file", type=str, help="/path/to/the file that contains all the mutations you want to make.")
parser.add_argument("utility_dir", type=str, help="where do your utility files live? Give me the directory.")
parser.add_argument("structure_dir", type=str, help="where do you want to dump the decoys made during this protocol?")
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
if input_args.mutation_file is not None:
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


# initialize Rosetta
initialize_rosetta()

# PyMOL_Mover
pmm = PyMOL_Mover()
pmm.keep_history( True )

# load the input pose
native_pose = load_pose( input_args.native_pdb_file )
low_E_native_pose = native_pose.clone()
#native_pose.pdb_info().name( "native_pose" )
#pmm.apply( native_pose )

# make a ScoreFunction
sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )
orig_fa_rep = sf.get_weight( fa_rep )

# collect all non-branch point residue numbers
moveable_residues = [ res_num for res_num in range( 1, native_pose.n_residue() + 1 ) if native_pose.residue( res_num ).is_branch_point() == False ]



######################
#### PACK AND MIN ####
######################
# the input native pose should have already gone through this type of pack min ONCE
# this does it again one more time so that it is on the same page as the mutants that get this treatment once after mutation
# pack all residues except for branch point residues
pack_rotamers_mover = make_RotamerTrialsMover( moveable_residues = moveable_residues,
                                               sf = sf,
                                               input_pose = low_E_native_pose,
                                               pack_radius = None )
pack_rotamers_mover.apply( low_E_native_pose )

# we packed the side chains, so minimize them (only residues marked by moveable_residues)
# keep the backbone the same as the crystal because 1) we can, 2) it's easier, 3) we want to see what mutations do to packing more so
mm = MoveMap()
mm.set_bb( False )
for res_num in moveable_residues:
    mm.set_chi( res_num, True )
# gradient min of native_pose
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
    min_mover.apply( low_E_native_pose )
    print "low_E_native_pose", sf( low_E_native_pose ), jj
    pmm.apply( low_E_native_pose )



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
decoy_name = main_structure_dir + mutation_filename

# data holders
mutation_names = []
dtotal_scores = []
nat_tot_score = sf( low_E_native_pose )
dinterface_scores = []
nat_intf_score = get_interface_score( 2, sf, low_E_native_pose )
dintf_sasa = []
nat_intf_sasa = calc_interface_sasa( low_E_native_pose, 2 ) 
binding_ratios = []
# True if dintf_score is negative if ratio > 1 and False if dintf_score is positive if ratio > 1 (or vice versa)
hits = []

# to hold all of the mutants to look at later
for mutation_set, ratio in all_mutations_to_ratio_dict.items():
    # get a fresh pose object by cloning the native_pose
    mutant_pose = native_pose.clone()
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
                                                   input_pose = low_E_native_pose,
                                                   pack_radius = None )
    pack_rotamers_mover.apply( low_E_native_pose )

    # we packed the side chains, so minimize them (only residues marked by moveable_residues)
    # keep the backbone the same as the crystal because 1) we can, 2) it's easier, 3) we want to see what mutations do to packing more so
    mm = MoveMap()
    mm.set_bb( False )
    for res_num in moveable_residues:
        mm.set_chi( res_num, True )
    # gradient min of native_pose
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
        min_mover.apply( low_E_native_pose )
        print "mutant_pose", sf( low_E_native_pose ), jj

    # dump the packed and minimized mutant_pose
    mutant_pose.dump_pdb( mutant_pose.pdb_info().name() )
    mutant_pose.pdb_info().name( mutation_set )

    ##############
    #### DATA ####
    ##############
    mutation_names.append( mutation_set )
    dtotal_scores.append( sf( mutant_pose ) - nat_tot_score )
    dintf_score = get_interface_score( 2, sf, mutant_pose ) - nat_intf_score
    dinterface_scores.append( dintf_score )
    dintf_sasa.append( calc_interface_sasa( mutant_pose, 2 ) - nat_intf_sasa )
    if ratio is not None:
        binding_ratios.append( ratio )
    else:
        binding_ratios.append( None )
    if ratio is not None:
        # if the ratio is higher than one and the dintf_score is negative, we got it
        if ratio > 1 and dintf_score < 0:
            hits.append( True )
        # if the ratio is higher than one and the dintf_score is positive, we missed it
        elif ratio > 1 and dintf_score > 0:
            hits.append( False )
        # if the ratio is less than one and the dintf_score is negative, we missed it
        elif ratio < 1 and dintf_score < 0:
            hits.append( False )
        # if the ratio is less than one and the dintf_score is positive, we got it
        elif ratio < 1 and dintf_score > 0:
            hits.append( True )
        elif ratio == 1:
            hits.append( None )
    else:
        hits.append( None )

# normalize the dinterface_scores and binding_ratios
#normalized_dinterface_scores = [ ( score - min( dinterface_scores ) ) / ( max( dinterface_scores ) - min( dinterface_scores ) ) for score in dinterface_scores ]
# need the opposite order of normalized_dinterface_scores because a value of 0 should correspond to a bad mut (and the reverse) which won't happen with regular normalization since the "good" mutations result in a negative dintf_score, ie a normalized value closer to 0
#normalized_reverse_dinterface_scores = [ 1 - score for score in normalized_dinterface_scores ]
#normalized_binding_ratios = [ ( score - min( binding_ratios ) ) / ( max( binding_ratios ) - min( binding_ratios ) ) if score is not None else None for score in binding_ratios ]
# now get the rank order of both lists because it is essentially impossible for the normalized values to line up perfectly
# won't work right if there is a None in the binding_ratios!
'''
rank_order_dinterface_scores = get_rank_order_of_list( dinterface_scores )
rank_order_binding_ratios = get_rank_order_of_list( binding_ratios )
'''
# reverse the rank order for dinterface_scores. Adding 1 because reversing it would make it go from 0 to n-1 even though it started at 1 to n
'''
reversed_rank_order_dinterface_scores = [ ( max( rank_order_dinterface_scores ) - score ) + 1 for score in rank_order_dinterface_scores ]
'''
# should probably normalize again
#normalized_rank_order_dinterface_scores = [ ( score - min( rank_order_dinterface_scores ) ) / ( max( rank_order_dinterface_scores ) - min( rank_order_dinterface_scores ) ) for score in rank_order_dinterface_scores ]
#normalized_rank_order_binding_ratios = [ ( score - min( rank_order_binding_ratios ) ) / ( max( rank_order_binding_ratios ) - min( rank_order_binding_ratios ) ) if score is not None else None for score in rank_order_binding_ratios ]
# reverse the normalization of the dintf_score because a low dintf should correspond to a high binding ratio
#reversed_normalized_rank_order_dinterface_scores = [ 1 - score for score in normalized_rank_order_dinterface_scores ]


try:
    import pandas as pd
    pandas_on = True
except ImportError:
    pandas_on = False
    pass

if pandas_on:
    df = pd.DataFrame()
    df[ "mutations" ] = mutation_names
    df[ "dtotal_score" ] = dtotal_scores
    df[ "dintf_score" ] = dinterface_scores
    df[ "dintf_sasa" ] = dintf_sasa
    df[ "binding_ratio" ] = binding_ratios
    #df[ "rank_dintf_score" ] = reversed_rank_order_dinterface_scores
    #df[ "rank_binding_ratio" ] = rank_order_binding_ratios
    df[ "hits" ] = hits
    df = df.sort_values( "dintf_score" )
else:
    data_dict = {}
    data_dict[ "mutations" ] = mutation_names
    data_dict[ "dtotal_score" ] = dtotal_scores
    data_dict[ "dintf_score" ] = dinterface_scores
    data_dict[ "dintf_sasa" ] = dintf_sasa
    data_dict[ "binding_ratio" ] = binding_ratios
    #data_dict[ "rank_dintf_score" ] = reversed_rank_order_dinterface_scores
    #data_dict[ "rank_binding_ratio" ] = rank_order_binding_ratios
    data_dict[ "hits" ] = hits
        
    '''
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
    #jd.output_decoy( working_pose )
    # zip up the decoy pose, if desired
    #if input_args.zip_decoy:
    #    os.popen( "gzip %s" %working_pose.pdb_info().name() )

    # increment the decoy number counter
    cur_decoy_num += 1

    #for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
    #    print "Res num", res_num, "Reset phi:", GlycanModelProtocol.reset_pose.phi( res_num ), "end phi:", working_pose.phi( res_num )
    #    print "Res num", res_num, "Reset psi:", GlycanModelProtocol.reset_pose.psi( res_num ), "end psi:", working_pose.psi( res_num )
    #    print "Res num", res_num, "Reset omega:", GlycanModelProtocol.reset_pose.omega( res_num ), "end omega:", working_pose.omega( res_num )
    #    print

# move the lowest E pack and minimized native structure into the lowest_E_structs dir
#fasc_filename = decoy_name + ".fasc"
#lowest_E_native_filename = get_lowest_E_from_fasc( fasc_filename, GlycanModelProtocol.lowest_E_structs_dir, 10 )
'''
