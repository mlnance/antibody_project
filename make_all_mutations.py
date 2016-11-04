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
#pmm.apply( native_pose )





#########################
#### JOB DISTRIBUTOR ####
#########################

# imports
from antibody_functions import read_mutation_file
from rosetta.core.scoring import fa_rep

# get the list of mutations to make
all_mutations = read_mutation_file( input_args.mutation_file )

# create an appropriate decoy_name using the mutation_file name
mutation_filename = input_args.mutation_file.split( '/' )[-1].split( '.' )[0]
decoy_name = main_structure_dir + mutation_filename

sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )
orig_fa_rep = sf.get_weight( fa_rep )



##########################
#### PYJOBDISTRIBUTOR ####
##########################

# imports
#from rosetta import PyJobDistributor

# create and use the PyJobDistributor object
#jd = PyJobDistributor( decoy_name, 1, sf )
#jd.native_pose = native_pose
#cur_decoy_num = 1

#while not jd.job_complete:
# name to use when dumping the decoy. Should include full path
#working_pose.pdb_info().name( jd.current_name )



# imports
from antibody_functions import mutate_residue, calc_interface_sasa, \
    get_interface_score
from rosetta import standard_packer_task, RotamerTrialsMover
from rosetta import MoveMap, MinMover, score_type_from_name, Pose
# data
mutation_names = []
dtotal_scores = []
nat_tot_score = sf( native_pose )
dinterface_scores = []
nat_intf_score = get_interface_score( 2, sf, native_pose )
dintf_sasa = []
nat_intf_sasa = calc_interface_sasa( native_pose, 2 ) 
mutant_dict = {}
for mutation_set in all_mutations:
    # get a fresh pose object
    mutant_pose = native_pose.clone()
    mutant_pose.pdb_info().name( main_structure_dir + mutation_set + ".pdb" )
    print "%s..." %mutation_set


    # make each mutation one at a time
    ###########################
    #### make the mutation ####
    ###########################
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
            continue

    ###########################
    #### gradient pack/min ####
    ###########################
    min_pose = Pose()
    for ii in range( 5 ):
        working_pose = mutant_pose.clone()

        for jj in range( 2 ):
            # packing
            task = standard_packer_task( working_pose )
            task.or_include_current( True )
            task.restrict_to_repacking()
            rtm = RotamerTrialsMover( sf, task )
            rtm.apply( working_pose )

            # minimizing
            mm = MoveMap()
            mm.set_bb( True )
            mm.set_chi( True )
            min_mover = MinMover( movemap_in = mm,
                                  scorefxn_in = sf,
                                  min_type_in = "dfpmin_strong_wolfe", 
                                  tolerance_in = 0.01,
                                  use_nb_list_in = True )
            min_mover.max_iter( 2500 )
            min_mover.apply( working_pose )
        '''
        for jj in range( 3 ):
            if jj == 0:
                sf.set_weight( fa_rep, sf.get_weight( fa_rep ) * 0.1 )
            elif jj == 1:
                sf.set_weight( fa_rep, sf.get_weight( fa_rep ) * 0.33 )
            else:
                sf.set_weight( fa_rep, orig_fa_rep )
            min_mover = MinMover( movemap_in = mm,
                                  scorefxn_in = sf,
                                  min_type_in = "dfpmin_strong_wolfe", 
                                  tolerance_in = 0.01,
                                  use_nb_list_in = True )
            min_mover.max_iter( 2500 )
            min_mover.apply( working_pose )
            print sf( working_pose ), jj
        '''
        if sf( working_pose ) < sf( min_pose ):
            min_pose.assign( working_pose )
    min_pose.dump_pdb( mutant_pose.pdb_info().name() )
    mutant_dict[ mutation_set ] = min_pose

    ##############
    #### data ####
    ##############
    mutation_names.append( mutation_set )
    dtotal_scores.append( sf( min_pose ) - nat_tot_score )
    dinterface_scores.append( get_interface_score( 2, sf, min_pose ) - nat_intf_score )
    dintf_sasa.append( calc_interface_sasa( min_pose, 2 ) - nat_intf_sasa )

                
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
else:
    data_dict = {}
    data_dict[ "mutations" ] = mutation_names
    data_dict[ "dtotal_score" ] = dtotal_scores
    data_dict[ "dintf_score" ] = dinterface_scores
    data_dict[ "dintf_sasa" ] = dintf_sasa
        
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
