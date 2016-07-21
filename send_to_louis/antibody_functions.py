#!/usr/bin/python
__author__ = 'morganlnance'

# TODO: move imports to the functions where they actually are needed



#####################
#### ALL IMPORTS ####
#####################

if __name__ == "__main__":
    print
    print "Importing modules..."
else:
    print
    print "Importing modules from %s..." %__name__

# bread and butter Rosetta imports
from rosetta import init, pose_from_file, \
    standard_packer_task, change_cys_state, \
    Pose, MoveMap, RotamerTrialsMover, MinMover, \
    PyMOL_Mover, AtomID, aa_from_oneletter_code, \
    FoldTree, get_fa_scorefxn, Vector1
from rosetta.utility import vector1_bool
from rosetta.numeric import xyzVector_Real
from rosetta.core.chemical import VariantType
from rosetta.core.pose import remove_variant_type_from_pose_residue
from rosetta.core.scoring.constraints import AtomPairConstraint, AngleConstraint
from toolbox import mutate_residue, get_hbonds

# for sugar work
#from rosetta.protocols.carbohydrates import GlycanRelaxMover
#from rosetta.protocols.carbohydrates import LinkageConformerMover

# for loop work
from rosetta import Loop, Loops, add_single_cutpoint_variant
from rosetta.protocols.loops.loop_closure.ccd import CCDLoopClosureMover

# for extra scoring functionality
from rosetta.core.scoring import score_type_from_name, CA_rmsd
from rosetta.core.scoring.func import HarmonicFunc, CircularHarmonicFunc
from rosetta.protocols.analysis import InterfaceAnalyzerMover as IAM

# import extras
import os
import sys
try:
    import pandas as pd
    pandas_on = True
except ImportError:
    pandas_on = False
    print "Skipping Pandas import - consider downloading it! Who doesn't love Pandas??"
    pass



# create global pymol object
# TODO-add visualization options to each relevant function using this global PyMOL_Mover
pmm = PyMOL_Mover()

# global data path  --  SPECIFIC TO MORGAN'S COMPUTER  --  CHANGE ON YOURS
data_dir = "/Users/Research/pyrosetta_git_repo/mutational_data/"

# global variables
AA_list = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
all_letters_list = [ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' ]
CUTOFF_DISTANCE = 5.0  # used when calculating the number of residue contacts at the interface
PACK_RADIUS = 20.0  # used when making mutations to structures, repacks in this area
PROBE_RADIUS = 1.4  # used for calculating total SASA
kT = 0.8  # used in MonteCarlo and small and shear movers



## Pose numbering information ONLY relevant to PDB 3ay4
# native
native_Fc_chain_A = range( 1, 215 + 1 )
native_Fc_glycan_A = range( 216, 223 + 1 )
native_Fc_chain_B = range( 224, 341 + 1 )
native_Fc_glycan_B = range( 440, 447 + 1 )
native_FcR_protein = range( 448, 607 + 1 )
native_FcR_main_glycan = range( 608, 615 + 1 )
native_FcR_three_mer = range( 616, 618 + 1 )
native_FcR_glycan = range( 608, 618 + 1 )
native_Fc_protein = []
native_Fc_protein.extend( native_Fc_chain_A )
native_Fc_protein.extend( native_Fc_chain_B )
native_Fc_glycan = []
native_Fc_glycan.extend( native_Fc_glycan_A )
native_Fc_glycan.extend( native_Fc_glycan_B )
native_order = range( 1, 618 + 1 )
native_Fc_protein_chains = [ 'A', 'B' ]
native_FcR_protein_chains = [ 'C' ]
native_Fc_glycan_chains = [ 'D', 'E', 'F', 'G' ]
native_FcR_glycan_chains = [ 'H', 'I', 'J', 'K' ]

# glycosylated decoy
decoy_Fc_chain_A = range( 1, 215 + 1 )
decoy_Fc_glycan_A = range( 603, 610 + 1)
decoy_Fc_chain_B = range( 216, 431 + 1 )
decoy_Fc_glycan_B = range( 611, 618 + 1)
decoy_FcR_protein = range( 432, 591 + 1)
decoy_FcR_main_glycan = range( 592, 598 + 1 )
decoy_FcR_three_mer = range( 599, 602 + 1 )
decoy_FcR_glycan = range( 592, 602 + 1 )
decoy_Fc_protein = []
decoy_Fc_protein.extend( decoy_Fc_chain_A )
decoy_Fc_protein.extend( decoy_Fc_chain_B )
decoy_Fc_glycan = []
decoy_Fc_glycan.extend( decoy_Fc_glycan_A )
decoy_Fc_glycan.extend( decoy_Fc_glycan_B )
decoy_order = []
decoy_order.extend( decoy_Fc_chain_A )
decoy_order.extend( decoy_Fc_glycan_A )
decoy_order.extend( decoy_Fc_chain_B )
decoy_order.extend( decoy_Fc_glycan_B )
decoy_order.extend( decoy_FcR_protein )
decoy_order.extend( decoy_FcR_main_glycan )
decoy_order.extend( decoy_FcR_three_mer )
decoy_Fc_protein_chains = [ 'A', 'B' ]
decoy_FcR_protein_chains = [ 'C' ]
decoy_Fc_glycan_chains = [ 'H', 'I', 'J', 'K' ]
decoy_FcR_glycan_chains = [ 'D', 'E', 'F', 'G' ]

# make an appropriate dictionary map
native_to_decoy_res_map = {}
for ii in range( len( native_order ) ):
    native_to_decoy_res_map[ native_order[ii] ] = decoy_order[ii]
decoy_to_native_res_map = {}
for ii in range( len( decoy_order ) ):
    decoy_to_native_res_map[ decoy_order[ii] ] = native_order[ii]




# adds these two functions to FoldTree so that now FoldTree has the functions and can use elsewhere
def _new_loop( loop, loop_buffer = 1):
    begin = loop.start() - loop_buffer
    end = loop.stop() + loop_buffer
    cut = loop.cut()
    new_jump( begin, end, cut )


def _new_loops( loops ):
    for ii in range( 1, loops.num_loop() + 1 ):
        loop = loops[ ii ]
        new_loop( loop )

FoldTree.new_loop = _new_loop
FoldTree.new_loops = _new_loops



##########################
#### WORKER FUNCTIONS ####
##########################

def initialize_rosetta():
    # called when this file is imported rather than ran directly
    print "Initializing Rosetta with sugar flags"

    # makes Rosetta quiet and sugar I/O ready
    #init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records" )
    init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -write_pdb_link_records" )



def load_pose( pose_filename, is_3ay4 = False ):
    """
    Load pose from a filename
    :param pose_filename: str( /path/to/pose/filename )
    :param is_3ay4: bool( if this is 3ay4 native or decoy, attached chain and residue information to Pose object )
    :return: a Rosetta Pose
    """
    # imports
    from rosetta import Pose, pose_from_file, FoldTree
    
    # create Pose object from filename
#    print "Loading pose"
    pose = Pose()
    pose_from_file( pose, pose_filename )
    
    # clean up the name of the pose
    pose_name = pose.pdb_info().name()
    pose_name = pose_name.split( '/' )[-1]
    pose.pdb_info().name( pose_name )

    # store the original FoldTree and add empty loops for use later
    pose.orig_fold_tree = FoldTree( pose.fold_tree() )
    pose.loops = None
    pose.loops_file = None

    # store residue and chain designation information for both native and decoy for 3ay4
    if is_3ay4:
        # native information
        pose.native_Fc_chain_A = native_Fc_chain_A
        pose.native_Fc_glycan_A = native_Fc_glycan_A
        pose.native_Fc_chain_B = native_Fc_chain_B
        pose.native_Fc_glycan_B = native_Fc_glycan_B
        pose.native_FcR_protein = native_FcR_protein
        pose.native_FcR_main_glycan = native_FcR_main_glycan
        pose.native_FcR_three_mer = native_FcR_three_mer
        pose.native_FcR_glycan = native_FcR_glycan
        pose.native_Fc_protein = native_Fc_protein
        pose.native_Fc_glycan = native_Fc_glycan
        pose.native_order = native_order
        pose.native_Fc_protein_chains = native_Fc_protein_chains
        pose.native_FcR_protein_chains = native_FcR_protein_chains
        pose.native_Fc_glycan_chains = native_Fc_glycan_chains
        pose.native_FcR_glycan_chains = native_FcR_glycan_chains

        # decoy information
        pose.decoy_Fc_chain_A = decoy_Fc_chain_A
        pose.decoy_Fc_glycan_A = decoy_Fc_glycan_A
        pose.decoy_Fc_chain_B = decoy_Fc_chain_B
        pose.decoy_Fc_glycan_B = decoy_Fc_glycan_B
        pose.decoy_FcR_protein = decoy_FcR_protein
        pose.decoy_FcR_main_glycan = decoy_FcR_main_glycan
        pose.decoy_FcR_three_mer = decoy_FcR_three_mer
        pose.decoy_FcR_glycan = decoy_FcR_glycan
        pose.decoy_Fc_protein = decoy_Fc_protein
        pose.decoy_Fc_glycan = decoy_Fc_glycan
        pose.decoy_order = decoy_order
        pose.decoy_Fc_protein_chains = decoy_Fc_protein_chains
        pose.decoy_FcR_protein_chains = decoy_FcR_protein_chains
        pose.decoy_Fc_glycan_chains = decoy_Fc_glycan_chains
        pose.decoy_FcR_glycan_chains = decoy_FcR_glycan_chains


    return pose



def get_fa_scorefxn_with_given_weights( weights_dict, verbose = False ):
    """
    Return an sf from get_fa_scoretype but with adjusted weights <scoretypes> with given <weights>
    If <input_scoretype> is not already part of the <sf>, this function will add it to <sf> with a weight of <weight>, and then get the score
    Will exit if the string( <input_scoretype> ) is not a valid ScoreType
    :param weights_dict: dict( ScoreType or str of ScoreType name : int( or float( weight ) ) )
    :param verbose: bool( print the final weights of the returned ScoreFunction? ) Default = False
    "return: ScoreFunction( fa_scorefxn with adjusted weights of given scoretypes )
    """
    # argument check - check the passed argument is a dict
    if not isinstance( weights_dict, dict ):
        print "You didn't give me a dictionary for your input. I need a dict of ScoreType (or name) : weight. Exiting."
        sys.exit()
        
    # get a standard fa_scorefxn to start with
    sf = get_fa_scorefxn()
    
    # make a dummy ScoreType to use for isinstance() checking
    fa_dun = score_type_from_name( "fa_dun" )

    # for each entry of the dictionary, change the weight
    for scoretype_name in weights_dict:
        # if the key is a string
        if isinstance( scoretype_name, str ):
            try:
                scoretype = score_type_from_name( scoretype_name )
            except:
                print "\nThe string name: '%s' does not appear to be a valid ScoreType. Exiting" %scoretype_name
                sys.exit()
            
            # get the corresponding weight
            weight = weights_dict[ scoretype_name ]

            # set the weight
            sf.set_weight( scoretype, weight )

        # if the argument is a ScoreType object
        elif isinstance( scoretype_name, type( fa_dun ) ):
            # adjust the weight in the scorefxn using the corresponding weight given
            sf.set_weight( scoretype_name, weights_dict[ scoretype_name ] )

        # else, I don't know what they gave me as a scoretype
        else: 
            print "I'm not sure what '%s' is from your ScoreType key in your <weights_dict> argument. Exiting" %scoretype_name
            sys.exit()

    # if verbose, print out the weights of the new ScoreFunction
    if verbose:
        print "\nNew score weights sf:\n%s\n" %( "\n".join( [ "%s: %s" %( str( name ), sf.get_weight( name ) ) for name in sf.get_nonzero_weighted_scoretypes() ] ) )

    # return the newly weighted fa_scorefxn
    return sf



def make_fa_scorefxn_from_file( scorefxn_file, verbose = False ):
    """
    Return an sf from get_fa_scoretype but with adjusted weights read from the passed <scorefxn_file>
    If a ScoreType is not already part of the sf, this function will add it to sf with the specified weight
    :param scorefxn_file: str( path/to/scorefxn.sf file with ScoreTypes and weights )
    :param verbose: bool( print the final weights of the returned ScoreFunction? ) Default = False
    "return: ScoreFunction( fa_scorefxn with adjusted weights of given scoretypes )
    """
    # make a standard fa_scorefxn
    sf = get_fa_scorefxn()

    # edit the ScoreFunction
    try:
        with open( scorefxn_file, "rb" ) as fh:
            score_lines = fh.readlines()
    except:
        print "Something seems to be wrong with your scorefxn_file ( %s ). Please check your input" %scorefxn_file
        sys.exit()
    for score_line in score_lines:
        # ignore new line characters
        score_line = score_line.rstrip()

        # get the score type and weight from each line, which should be space delimited
        if score_line != '':
            # ignore commented lines
            if not score_line.startswith( '#' ):
                # check if the user followed the right format
                try:
                    score_line_split = score_line.split( ' ' )
                    score_type = str( score_line_split[0] )

                    # if user wants a multiplier
                    if score_line_split[1] == '*':
                        score_weight = sf.get_weight( score_type_from_name( score_type ) ) * float( score_line_split[2] )

                    # else user wants a specific value
                    else:
                        score_weight = float( score_line_split[1] )
                except:
                    print "\nIt seems that your scorefxn_file did not follow the proper format. Please follow 'ScoreType Weight' or 'ScoreType * Multiplier'"
                    sys.exit()
        try:
            sf.set_weight( score_type_from_name( score_type ), score_weight )
        except:
            print "It seems like you did not pass a valid ScoreType (or some other issue). Check your scorefxn_file (%s)" %input_args.scorefxn_file
            sys.exit()

    # if verbose, print the final weights
    if verbose:
        print "\nScore weights used in this sf:\n%s\n" %( "\n".join( [ "%s: %s" %( str( name ), sf.get_weight( name ) ) for name in sf.get_nonzero_weighted_scoretypes() ] ) )

    return sf




def show_score_breakdown( sf, pose ):
    """
    Shows the breakdown of the <pose>'s total score by printing the score of each nonzero weighted ScoreType in <sf>
    :param sf: ScoreFunction
    :param pose: Pose
    """
    # print out each score
    print "\n".join( [ "%s: %s" %( score_type, round( sf.score_by_scoretype( pose, score_type ), 3 ) ) for score_type in sf.get_nonzero_weighted_scoretypes() ] )
    print



def show_sf_weights_breakdown( scorefxn ):
    """
    Shows the breakdown of nonzeroweighted ScoreTypes in <scorefxn>
    :param sf: ScoreFunction
    :param pose: Pose
    """
    # print out each score
    print
    print "\n".join( [ "%s: %s" %( score_type, scorefxn.get_weight( score_type ) ) for score_type in scorefxn.get_nonzero_weighted_scoretypes() ] )
    print



def get_score_by_scoretype( sf, input_scoretype, pose, weight = 1.0, verbose = False ):
    """
    Return the specified <input_scoretype> value using <sf> on the <pose>
    If <input_scoretype> is not already part of the <sf>, this function will add it to <sf> with a weight of <weight>, and then get the score
    Will exit if the string( <input_scoretype> ) is not a valid ScoreType
    :param sf: ScoreFunction
    :param input_scoretype: str( name of a valid ScoreType ) or ScoreType object
    :param pose: Pose
    :param weight: int( or float( weight of the desired <input_scoretype> if it's not already a part of <sf> ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( <input_scoretype> score of <pose> )
    """
    if verbose:
        print "Getting the total score of a specific ScoreType for the pose"

    # make a dummy ScoreType to use for isinstance() checking
    fa_dun = score_type_from_name( "fa_dun" )

    # check to make sure the string <input_scoretype> is a valid ScoreType
    if isinstance( input_scoretype, str ):
        try:
            scoretype = score_type_from_name( input_scoretype )
        except:
            print
            print
            print "The string name:", input_scoretype, "does not appear to be a valid ScoreType. Exiting"
            sys.exit()

    # if the input was a ScoreType, it must be valid. Turn its name into a string version for use
    elif isinstance( input_scoretype, type( fa_dun ) ):
        scoretype = input_scoretype
        input_scoretype = str( input_scoretype )

    # else, I don't know what they gave me
    else:
        print "I'm not sure what '%s'is. Exiting" %scoretype
        sys.exit()

    # tell user it was successful in turning str( <input_scoretype> ) to a ScoreType object
    if verbose:
        print "Successfully turned", input_scoretype, "into a valid ScoreType object"

    # check to see if <sf> already has the given <input_scoretype>
    has_scoretype = False
    current_scoretypes = sf.get_nonzero_weighted_scoretypes()
    for score_type in current_scoretypes:
        if str( score_type ) == input_scoretype:
            has_scoretype = True
            break

    # if has_scoretype is False, add the <input_scoretype> to the <sf> with weight <weight>
    if not has_scoretype:
        print "It appears the ScoreType is not already included in your ScoreFunction, adding", input_scoretype, "with a weight of", weight
        sf.set_weight( scoretype, weight )

    # score the <pose> with the specified <input_scoretype> and return the value
    score = sf.score_by_scoretype( pose, scoretype )

    return score



def get_residue_score_by_scoretype( sf, input_scoretype, seq_pos, pose, weight = 1.0, verbose = False ):
    """
    Given a <pose>, a specific <seq_pos> for a residue, a <input_scoretype> of type string or ScoreType, return the score of that residue using the ScoreFunction <sf>
    If the <sf> does not already have <input_scoretype> as a non-zero weighted ScoreType, then this function will add it with a weight of <weight> and then score the residue
    Will exit and print an error message if the string <input_scoretype> is invalid
    :param sf: ScoreType
    :param input_scoretype: str( a valid ScoreType name ) or a ScoreType object
    :param seq_pos: int( the residue of interest )
    :param pose: Pose
    :param weight: int( or float( the weight to be used if <input_scoretype> is not already in <sf> ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( the <input_scoretype> score of residue <seq_pos> )
    """
    if verbose:
        print "Getting and comparing each residue between the two poses based on ScoreType"

    # make a dummy ScoreType to use for isinstance() checking
    fa_dun = score_type_from_name( "fa_dun" )

    # check to make sure the string <input_scoretype> is a valid ScoreType
    if isinstance( input_scoretype, str ):
        try:
            scoretype = score_type_from_name( input_scoretype )
        except:
            print
            print
            print "The string name:", input_scoretype, "does not appear to be a valid ScoreType. Exiting"
            sys.exit()

    # if the input was a ScoreType, it must be valid. Turn its name into a string version for use
    elif isinstance( input_scoretype, type( fa_dun ) ):
        scoretype = input_scoretype
        input_scoretype = str( input_scoretype )

    # else, I don't know what they gave me
    else:
        print "I'm not sure what", input_scoretype, "is. Exiting"
        sys.exit()

    # tell user it was successful going from str( <input_scoretype> ) to ScoreType
    if verbose:
        print "Successfully turned", input_scoretype, "into a valid ScoreType object"

    # check to see if <sf> already has the given <input_scoretype>
    has_scoretype = False
    current_scoretypes = sf.get_nonzero_weighted_scoretypes()
    for score_type in current_scoretypes:
        if str( score_type ) == input_scoretype:
            has_scoretype = True
            break

    # if has_scoretype is False, add the <input_scoretype> to the <sf> with weight <weight>
    if not has_scoretype:
        if verbose:
            print "It appears the ScoreType is not already included in your ScoreFunction, adding", input_scoretype, "it with a weight of", weight
        sf.set_weight( scoretype, weight )

    # get access to <pose>'s energies() object
    sf( pose )

    # get the <input_scoretype> score of the residue at <seq_pos>
    score_of_res_of_scoretype = pose.energies().residue_total_energies( seq_pos ).get( scoretype )
    
    
    # return the score of the residue given the ScoreType given
    return score_of_res_of_scoretype



def get_sugar_bb_only_sf( weight = 1 ):
    """
    Creates and returns a ScoreFunction with only the sugar_bb term as a non-zero. Can specify weight with <weight> param
    :param weight: int( or float( weight of the sugar_bb term in the sf ) ). Default = 1
    :return: ScoreFunction
    """
    # instantiate a fa ScoreFunction
    sugar_sf = get_fa_scorefxn()
    
    # get a list of all the ScoreTypes in the fa_scorefxn
    scoretypes = sugar_sf.get_nonzero_weighted_scoretypes()

    # set all ScoreTypes to 0
    for scoretype in scoretypes:
        sugar_sf.set_weight( scoretype, 0 )
            
    # set sugar_bb ScoreType to given <weight>
    # doing this outside the loop in case "sugar_bb" stops becoming a part of fa_scorefxn
    sugar_bb = score_type_from_name( "sugar_bb" )
    sugar_sf.set_weight( sugar_bb, weight )
    
    
    # return the sugar_sf that now only has the sugar_bb ScoreType as non-zero
    return sugar_sf



def apply_sugar_constraints_to_sf( sf, pose, weight = 1.0, verbose = False ):
    """
    Applies bond distance and bond angle constraints to sugar branch points within <pose> with a weight of <weight>
    Uses a HarmonicFun for the bond distance, and a CircularHarmonicFun for the bond angle
    :param sf: ScoreFunction
    :param pose: Pose
    :param weight: int( or float( what do you want the weight to be for the bond distance and bond angle constraints? ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: ScoreFunction( sf including sugar constraints )
    """
    # get list of chemical edges ( marked as -2 in FoldTree f )
    edges = pose.fold_tree().get_chemical_edges()

    if verbose:
        print "Constraining the bond distance and bond angles for", len( edges ), "sugar residues"

    # need atom IDs from residue numbers
    for ii in range( 1, len( edges ) + 1 ):
        # extract info from edge
        start_res_num = edges[ii].start()
        start_atom_name = edges[ii].start_atom()
        stop_res_num = edges[ii].stop()
        stop_atom_name = edges[ii].stop_atom()

        # get atom id's of upper and lower branch point atoms
        upper_branch_atom_index = pose.residue( start_res_num ).atom_index( start_atom_name )
        lower_branch_atom_index = pose.residue( stop_res_num ).atom_index( stop_atom_name )

        # get atom index numbers of the first adjacent heavy atoms to each branch point atom
        upper_branch_angle_atom_index = pose.residue( start_res_num ).first_adjacent_heavy_atom( upper_branch_atom_index )
        lower_branch_angle_atom_index = pose.residue( stop_res_num ).first_adjacent_heavy_atom( lower_branch_atom_index )

        # make AtomID objects
        upper_branch_atom = AtomID( upper_branch_atom_index, start_res_num )
        lower_branch_atom = AtomID( lower_branch_atom_index, stop_res_num )
        upper_branch_angle_atom = AtomID( upper_branch_angle_atom_index, start_res_num )
        lower_branch_angle_atom = AtomID( lower_branch_angle_atom_index, stop_res_num )

        # get bond length of the atoms to constrain
        dist_upper_lower = pose.conformation().bond_length( upper_branch_atom, lower_branch_atom )

        # make and apply harmonic function
        harmonic_func = HarmonicFunc( dist_upper_lower, 0.5 )
        atm_pair_constraint = AtomPairConstraint( upper_branch_atom, lower_branch_atom, harmonic_func )
        pose.add_constraint( atm_pair_constraint )

        # get bond angle of the atoms to constrain protein C, N - C sugar
        ang_upper_angle_upper_lower = pose.conformation().bond_angle( upper_branch_angle_atom, upper_branch_atom, lower_branch_atom )

        # make and apply circular harmonic function
        circ_harmonic_func = CircularHarmonicFunc( ang_upper_angle_upper_lower, 5 )
        ang_constraint = AngleConstraint( upper_branch_angle_atom, upper_branch_atom, lower_branch_atom, circ_harmonic_func )
        pose.add_constraint( ang_constraint )

        # get bond angle of the atoms to constrain protein N - C, O sugar
        ang_upper_lower_lower_angle = pose.conformation().bond_angle( upper_branch_atom, lower_branch_atom, lower_branch_angle_atom )

        # make and apply circular harmonic function
        circ_harmonic_func = CircularHarmonicFunc( ang_upper_lower_lower_angle, 5 )
        ang_constraint = AngleConstraint( upper_branch_atom, lower_branch_atom, lower_branch_angle_atom, circ_harmonic_func )
        pose.add_constraint( ang_constraint )

    # add constraints to score function weights
    sf.set_weight( score_type_from_name( "atom_pair_constraint" ), weight )
    sf.set_weight( score_type_from_name( "angle_constraint" ), weight )

    return sf



def SugarSmallMover( seqpos, in_pose, angle_max, set_phi = True, set_psi = True, set_omega = True ):
    """
    Randomly resets the phi, psi, and omega values of the sugar residue <seqpos> in <pose> to old_value +/- angle_max/2
    Emulates the SmallMover but with the additional omega mover
    :param seqpos: int( the pose number for the residue )
    :param in_pose: Pose
    :param angle_max: int( or float( the max angle around which the phi/psi/omega could move ) )
    :param set_phi: bool( do you want to change the phi angle? ) Default = True
    :param set_psi: bool( do you want to change the psi angle? ) Default = True
    :param set_omega: bool( do you want to change the omega angle? ) Default = True
    :return: Pose
    """
    # imports
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg

    # copy the input pose
    pose = Pose()
    pose.assign( in_pose )

    # from BackboneMover.cc file for SmallMover
    big_angle = angle_max
    small_angle = big_angle / 2.0

    # get current phi, psi, and omega
    old_phi = pose.phi( seqpos )
    old_psi = pose.psi( seqpos )
    old_omega = pose.omega( seqpos )

    # get random values for phi, psi, and omega
    new_phi = periodic_range( old_phi - small_angle + rg().uniform() * big_angle, 360.0 )
    new_psi = periodic_range( old_psi - small_angle + rg().uniform() * big_angle, 360.0 )
    new_omega = periodic_range( old_omega - small_angle + rg().uniform() * big_angle, 360.0 )

    # set the new values
    if set_phi:
        pose.set_phi( seqpos, new_phi )
    if set_psi:
        pose.set_psi( seqpos, new_psi )
    if set_omega:
        pose.set_omega( seqpos, new_omega )

    return pose


def make_pack_rotamers_mover( sf, input_pose, apply_sf_sugar_constraints = True, pack_branch_points = True, residue_range = None, use_pack_radius = False, pack_radius = PACK_RADIUS, verbose = False ):
    """
    Returns a standard pack_rotamers_mover restricted to repacking and allows for sampling of current residue conformations for the <pose>
    IMPORTANT: DON'T USE for a Pose you JUST mutated  --  it's not setup to handle mutations
    If you want a packer task for a specific set of residues, set <residue_range> to a list of the residue sequence positions of interest
    If you have one or more residues of interest where you want to pack within a radius around each residue, set <use_pack_radius> to True, give a <pack_radius> (or use default), and set <residue_range> to the residue sequence positions of interest
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param pack_branch_points: bool( allow packing at branch points? ). Default = True
    :param residue_range: list( of int( valid residue sequence positions ) )
    :param use_pack_radius: bool( Do you want to pack residues within a certain <pack_radius> around residues specified in <residue_range>? ). Default = False
    :param pack_radius: int( or float( the radius around which you want to pack additional residues around the residues from <residue_range>. MUST have <set use_pack_radius> to True to do this ). Default = PACK_RADIUS = 20.0 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: a pack_rotamers_mover object
    """
    if verbose:
        print "Making a pack rotamers mover"
        
    # copy a fresh pose
    pose = Pose()
    pose.assign( input_pose )

    # check to make sure <residue_range> is a list of valid residue numbers
    if residue_range is not None:
        if isinstance( residue_range, list ):
            # quick check to make sure the numbers are valid residue sequence positions
            for res_num in residue_range:
                if not 0 < res_num <= pose.n_residue():
                    print
                    print "It seems that", res_num, "is not a valid sequence position for the Pose you gave me. Exiting"
                    sys.exit()
        else:
            # not sure what type of thing they gave me, print an error and exit
            print
            print "I'm not sure what you gave me. <residue_range> is of type", type( residue_range ), "and I needed a list. Even if it's just one residue (sorry). Exiting"
            sys.exit()

    # make the packer task
    task = standard_packer_task( pose )
    task.or_include_current( True )
    task.restrict_to_repacking()

    # if specified, turn off packing for each branch point residue
    if pack_branch_points is False:
        if verbose:
            print "  Turning off packing for branch points"
        for res_num in range( 1, pose.n_residue() + 1 ):
            if pose.residue( res_num ).is_branch_point():
                task.nonconst_residue_task( res_num ).prevent_repacking()
                
    # if <use_pack_radius> is True
    if use_pack_radius:
        # container for the centers of each residue from <residue_range>
        centers_of_res_from_residue_range = []

        # fill up the centers container
        for res_num in residue_range:
            # add all the centers of the residues from <residue_range> into the center container  --  used to find residues that should not be repacked
            center = pose.residue( res_num ).nbr_atom_xyz()
            centers_of_res_from_residue_range.append( center )

        # container for residues outside the <pack_radius> range  --  to be prevented from repacking
        res_nums_outside_pack_radius = []

        # make a counter to print out the number of residues that will be packed  --  if user turned verbose to True
        counter = 0

        # find all residues outside <pack_radius> of each residue from <residue_range>  --  they will be prevented from repacking
        for res_num in range( 1, pose.n_residue() + 1 ):
            # don't need to do this for the residues from <residue_range> since they should be packed
            if res_num not in residue_range:
                # get the center of the residue
                center = pose.residue( res_num ).nbr_atom_xyz()

                # boolean used to ensure that the residue is outside the <pack_radius> for every residue from <residue_range>
                outside_pack_radius = True

                # loop over each center of the residues from <residue_range>
                for center_res_range in centers_of_res_from_residue_range:
                    if center.distance( center_res_range ) <= pack_radius:
                        # if its within the distance, change outside_pack_radius to False and break from this loop
                        outside_pack_radius = False
                        break

                # if the center of the residue in question is outside the <pack_radius> from all of the residues from <residue_range>, add it to the residue container to have its repacking be prevented
                if outside_pack_radius:
                    res_nums_outside_pack_radius.append( res_num )

                    # up the counter
                    counter += 1

        if verbose:
            print "  Preventing", counter, "residues from repacking in this packer task"

        # turn off repacking for the residues not in the <pack_radius>
        for res_num in res_nums_outside_pack_radius:
            task.nonconst_residue_task( res_num ).prevent_repacking()
            
    # else, if no <pack_radius> given, but a <residue_range> was given, prevent repacking for the other residues
    else:
        if residue_range is not None:
            # prevent repacking for every residue NOT given by the user through <residue_range>
            for res_num in range( 1, pose.n_residue() + 1 ):
                if res_num not in residue_range:
                    task.nonconst_residue_task( res_num ).prevent_repacking()

    # apply sugar branch point constraints to sf, if desired
    if apply_sf_sugar_constraints:
        apply_sugar_constraints_to_sf( sf, pose )

    # make the pack rotamers mover and return it
    pack_rotamers_mover = RotamerTrialsMover( sf, task )

    return pack_rotamers_mover



def make_min_mover( sf, input_pose, apply_sf_sugar_constraints = True, jumps = None, allow_sugar_chi = False, minimization_type = "dfpmin_strong_wolfe", verbose = False ):
    """
    Returns a min_mover object suitable for glycosylated proteins (turns off carbohydrate chi for now)
    IMPORTANT: If using a specific type of minimization, <minimization_type>, it DOES NOT check beforehand if the type you gave is valid, so any string actually will work. It will break when used
    See 'https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview' for minimization type options
    Prints error message and exits if there was a problem
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param jumps: list( Jump numbers of Jump(s) to be minimized ). Default = None (ie. all Jumps) (Give empty list for no Jumps)
    :param allow_sugar_chi: bool( allow the chi angles of sugars to be minimized ). Default = False
    :param minimization_type: str( the type of minimization you want to use ). Default = "dfpmin_strong_wolfe"
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: min_mover object
    """
    if verbose:
        print "Making a min mover"
        
    # copy a fresh pose
    pose = Pose()
    pose.assign( input_pose )

    # instantiate a MoveMap
    mm = MoveMap()

    # set all Jumps to True if <jumps> is None
    if jumps is None:
        if verbose:
            print "  Setting all Jumps to be minimized"
        for ii in range( 1, pose.num_jump() + 1 ):
            mm.set_jump( ii, True )

    # otherwise, set the user-specified Jumps to True
    else:
        # make sure <jumps> is of type list
        if isinstance( jumps, list ):
            for ii in jumps:
                if verbose:
                    print "  Setting only", jumps, "to be minimized"
                # ensure <jumps> only contains valid Jumps for the <pose>
                if 0 < ii <= pose.num_jump():
                    mm.set_jump( ii, True )
                else:
                    print
                    print ii, "doesn't appear to be a valid Jump for your pose, exiting"
                    sys.exit()
        else:
            print
            print "You didn't give me a list for the jumps param - I need a list of valid jump numbers. Exiting"
            sys.exit()

    # set backbone angles to be minimized for both protein and sugar residues
    mm.set_bb( True )

    # turn on chi angle minimization for all amino acids
    # turn on or off chi angle minimization for sugars based on user input
    if allow_sugar_chi:
        if verbose:
            print "  Turning on chi angle mobility for all residues"
        for residue in pose:
            mm.set_chi( residue.seqpos(), True )
    else:
        if verbose:
            print "  Turning on chi angle mobility only for protein residues"
        for residue in pose:
            if not residue.is_carbohydrate():
                mm.set_chi( residue.seqpos(), True )

    # apply sugar branch point constraints to sf, if desired
    if apply_sf_sugar_constraints:
        apply_sugar_constraints_to_sf( sf, pose )

    # create a MinMover with options
    min_mover = MinMover( mm, sf, minimization_type, 0.01, True )

    return min_mover



def make_movemap_for_loop( loop, allow_bb_movement = True, allow_chi_movement = True, verbose = False ):
    """
    Given a loop, creates and returns a MoveMap allowing for bb and chi angle freedom, with a +/- 1 buffer
    :param loop: Loop( loop of interest )
    :param allow_bb_movement: bool( Do you want to allow BackBone movement for the loop? ). Default = True
    :param allow_chi_movement: bool( Do you want to allow chi angle movement for the loop? ). Default = True
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for Loop
    """
    if verbose:
        print "Making a MoveMap for a Loop"
    
    # instantiate a MoveMap
    mm = MoveMap()

    # set only the Loop area to true
    # if BB is True
    if allow_bb_movement:
        mm.set_bb_true_range( loop.start() + 1, loop.stop() - 1 )

    # if chi is True
    if allow_chi_movement:
        mm.set_chi_true_range( loop.start() + 1, loop.stop() - 1 )

    return mm



def make_movemap_for_jumps( JUMP_NUMbers, verbose = False ):
    """
    Given a jump number, creates and returns a MoveMap allowing for movement of only the Jump(s) given in the list
    :param JUMP_NUM: list( int( valid Jump number ) )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for Jump(s)
    """
    # instantiate a MoveMap
    mm = MoveMap()

    # print different verbose statements and make the MoveMap differently depending on if it's a single number, or a list
    if isinstance( JUMP_NUMbers, int ):
        # this is only one Jump
        if verbose:
            print "Making a MoveMap for Jump number", JUMP_NUMbers

        # set only the given Jump to True
        mm.set_jump( JUMP_NUMbers, True )

    elif isinstance( JUMP_NUMbers, list ):
        # this is a list of Jump numbers
        if verbose:
            print "Making a MoveMap for Jump numbers", JUMP_NUMbers

        # set only the given Jumps to True
        for JUMP_NUM in JUMP_NUMbers:
            mm.set_jump( JUMP_NUM, True )

    else:
        # I don't know what they gave me  -  I need an integer or a list
        print
        print "I'm not sure what ", JUMP_NUMbers, "is. I need an integer or a list of integers. Exiting"
        sys.exit()

    return mm



def make_movemap_for_sugars( pose, allow_bb_movement = True, allow_chi_movement = False, verbose = False ):
    """
    Given a <pose> object, return a MoveMap allowing for bb and/or chi movement of only sugar residues.
    :param pose: Pose( a Pose object )
    :param allow_bb_movement: bool( Do you want to allow sugarback bone movement? ). Default = True
    :param allow_chi_movement: bool( Do you want to allow sugarchi angle movement? ). Default = False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for all sugar residues in <pose>
    """
    if verbose:
        print "Making a MoveMap for a Loop"
    
    # instantiate a MoveMap
    mm = MoveMap()

    ## set movement for only sugar residues
    # if BB is True
    if allow_bb_movement:
        for res in pose:
            if res.is_carbohydrate():
                mm.set_bb( res.seqpos(), True )

    # if chi is True
    if allow_chi_movement:
        for res in pose:
            if res.is_carbohydrate():
                mm.set_chi( res.seqpos(), True )
    
    # return the sugar MoveMap
    return mm



def make_movemap_for_range( seqpos, allow_bb_movement = True, allow_chi_movement = True, verbose = False ):
    """
    Given a list of <seqpos>, return a MoveMap allowing for bb and/or chi movement of only the passed residues.
    :param seqpos: list( int( a list of sequence positions for your pose ) )
    :param allow_bb_movement: bool( Do you want to allow sugarback bone movement? ). Default = True
    :param allow_chi_movement: bool( Do you want to allow sugarchi angle movement? ). Default = False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for all sugar residues in <pose>
    """
    # check to make sure <seqpos> is a list
    if not isinstance( seqpos, list ):
        print "You didn't pass me a list for your <seqpos> argument. That's what I need to make your MoveMap. Exiting."
        sys.exit()
        
    if verbose:
        print "Making a MoveMap for the following residues:", seqpos
    
    # instantiate a MoveMap
    mm = MoveMap()

    # if BB is True
    if allow_bb_movement:
        for num in seqpos:
            mm.set_bb( num, True )

    # if chi is True
    if allow_chi_movement:
        for num in seqpos:
            mm.set_chi( num, True )
    
    # return the MoveMap
    return mm



# TODO: Do I actually need this function?
def do_min_with_this_mm( mm, sf, pose, apply_sf_sugar_constraints = True, minimization_type = "dfpmin_strong_wolfe", verbose = False ):
    """
    Minimizes a given Pose using dfpmin_strong_wolfe and the user-supplied ScoreFunction <sf> and MoveMap <mm>
    IMPORTANT: If using a specific type of minimization, <minimization_type>, it DOES NOT check beforehand if the type you gave is valid, so any string actually will work. It will break when used
    See 'https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview' for minimization type options
    :param mm: finished MoveMap
    :param sf: ScoreFunction
    :param pose:  Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param minimization_type: str( the type of minimization you want to use ). Default = "dfpmin_strong_wolfe"
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: minimized Pose
    """
    # apply sugar branch point constraints to sf, if desired
    if apply_sf_sugar_constraints:
        apply_sugar_constraints_to_sf( sf, pose )

    # create a MinMover with options
    min_mover = MinMover( mm, sf, minimization_type, 0.01, True )

    if verbose:
        print "Minimizing the Pose with the following min_mover"
        print min_mover

    # apply the min_mover
    min_mover.apply( pose )

    return pose



def do_pack_min( sf, input_pose, apply_sf_sugar_constraints = True, residue_range = None, jumps = None, pack_branch_points = True, use_pack_radius = False, pack_radius = PACK_RADIUS, allow_sugar_chi = False, minimization_type = "dfpmin_strong_wolfe", verbose = False, pmm = None ):
    """
    Makes and applies a packer task and basic min mover to <pose> using the supplied ScoreFunction <sf>
    IMPORTANT: If using a specific type of minimization, <minimization_type>, it DOES NOT check beforehand if the type you gave is valid, so any string actually will work. It will break when used
    See 'https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview' for minimization type options
    min_mover will set SUGAR BB freedom only - not CHI
    Will minimize over all Jumps, unless <jumps> is set to a list of Jump numbers to minimize instead
    If you want a packer task for a specific set of residues, set <residue_range> to a list of the residue sequence positions of interest
    If you have one or more residues of interest where you want to pack within a certain radius around each residue, set <use_pack_radius> to True, give a <pack_radius> (or use default value), and set <residue_range> to the residue sequence positions of interest
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param residue_range: list( of int( valid residue sequence positions ) )
    :param jumps = list( valid Jump numbers to be minimized if not all Jumps should be minimized). Default = None (ie. all jumps minimized) (Give empty list for no Jumps)
    :param pack_branch_points: bool( allow packing at branch points? ). Default = True
    :param use_pack_radius: bool( Do you want to pack residues within a certain <pack_radius> around residues specified in <residue_range>? ). Default = False
    :param pack_radius: int( or float( the radius around which you want to pack additional residues around the residues from <residue_range>. MUST have <set use_pack_radius> to True to do this ). Default = PACK_RADIUS = 20.0 Angstroms
    :param allow_sugar_chi: bool( allow chi minimization for sugars? ). Default = False
    :param minimization_type: str( the type of minimization you want to use ). Default = "dfpmin_strong_wolfe"
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: packed and minimized Pose
    """
    if verbose:
        print "Packing and minimizing the Pose"
        
    # get a fresh copy of the input_pose
    pose = Pose()
    pose.assign( input_pose )
        
    # apply the PyMOL_Mover if passed
    if pmm is not None:
        # used in the function elsewhere to ensure the program doesn't use a broken mover
        pmm_worked = False
        try:
            pmm.apply( pose )
            pmm_worked = True
        except:
            print "Something was wrong with your PyMOL_Mover -- continuing without watching"
            pass

    # make and apply the pack_rotamers_mover
    pack_rotamers_mover = make_pack_rotamers_mover( sf, 
                                                    pose, 
                                                    apply_sf_sugar_constraints = apply_sf_sugar_constraints, 
                                                    pack_branch_points = pack_branch_points, 
                                                    residue_range = residue_range, 
                                                    use_pack_radius = use_pack_radius, 
                                                    pack_radius = pack_radius, 
                                                    verbose = verbose )
    pack_rotamers_mover.apply( pose )
    if pmm is not None and pmm_worked:
        pmm.apply( pose )

    # make and apply the min_mover
    min_mover = make_min_mover( sf, 
                                pose, 
                                apply_sf_sugar_constraints = apply_sf_sugar_constraints, 
                                jumps = jumps, 
                                allow_sugar_chi = allow_sugar_chi,
                                minimization_type = minimization_type, 
                                verbose = verbose )
    min_mover.apply( pose )
    if pmm is not None and pmm_worked:
        pmm.apply( pose )

    return pose



# used in functions within Protocols code block
def ramp_score_weight( sf, score_type_str, target_weight, current_step, total_steps, verbose = False ):
    """
    Ramps up or down a given ScoreType toward its target score. Uses percentage of completion to calculate amount of steps
    Current and Total steps -- Say you're doing 100 rounds, if you're on round 57, current_step = 57, total_steps = 100
    :param sf: ScoreFunction
    :param score_type_str: str( valid ScoreType )
    :param target_weight: int( or float( value of your target ScoreType weight )
    :param current_step: int( current step of how many rounds you're doing )
    :param total_steps: int( number of rounds you're doing )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: the ScoreFunction with the adjusted weight
    """
    # get the ScoreType object given the str value passed in
    try:
        score_type = score_type_from_name( score_type_str )
    except RuntimeError:
        print
        print
        print score_type_str, "is not a valid score type, exiting"
        sys.exit()
    except:
        print
        print
        print "I'm not sure what happened. Check your input to ramp_score_weight. Exiting"
        sys.exit()

    # need to be able to do this a minimum of 10 times, otherwise the math will break
    if total_steps < 10:
        print
        print "You need to run a loop at least more than 10 times - otherwise the math in here won't work"
        sys.exit()

    # adjust the total_steps so that it's actually 10% less than what it actually is
    # this is so that the final 10% of the simulation will run with the target score weights
    total_steps = int( total_steps * 0.9 )

    # get the current weight of the score type
    current_weight = sf.get_weight( score_type )

    # make appropriate weight adjustments
    if current_weight > target_weight:  # we're decreasing the weight
        amount_left = target_weight - current_weight
        moves_left = float( total_steps - current_step )
        add_weight = amount_left / moves_left  # should be negative

        if verbose:
            print "Changing", score_type_str, "from", current_weight, "to", current_weight + add_weight

        # set the new weight
        sf.set_weight( score_type, current_weight + add_weight )

    elif target_weight > current_weight:  # we're increasing the weight
        amount_left = target_weight - current_weight
        moves_left = float( total_steps - current_step )
        add_weight = amount_left / moves_left  # should be positive

        if verbose:
            print "Changing", score_type_str, "from", current_weight, "to", current_weight + add_weight

        # set the new weight
        sf.set_weight( score_type, current_weight + add_weight )

    else:
        # we're at the target weight, so keep it there
        if verbose:
            print score_type_str, "is currently at the target weight of", target_weight

        # set to target weight
        sf.set_weight( score_type, target_weight )

    return sf



def CCD_loop_closure( loop, pose ):
    """
    Closes the given <loop> in <pose> using CCD
    :param loop: Loop( loop of interest )
    :param pose: Pose
    :return: Pose( <pose> with closed <loop> )
    """
    # make a MoveMap allowing only the loop to be flexible
    mm = make_movemap_for_loop( loop )

    # close the loop using CCD
    ccd = CCDLoopClosureMover( loop, mm )
    ccd.apply( pose )

    return pose



def calc_distance( vec1, vec2 ):
    """
    Calculates the distances between two points in 3D space
    :param vec1: list( x, y, z coordinates ) of point 1
    :param vec2: list( x, y, z coordinates ) of point 2
    :return: float( distance )
    """
    from math import sqrt, pow

    # assert that two lists were passed for coordinates
    if not isinstance( vec1, list ):
        print
        print "vec1 isn't a list - give me a list of xyz coordinates. Exiting"
        sys.exit()
    if not isinstance( vec2, list ):
        print
        print "vec2 isn't a list - give me a list of xyz coordinates. Exiting"
        sys.exit()

    # takes in two lists of xyz coordinates
    x1 = vec1[0]
    y1 = vec1[1]
    z1 = vec1[2]
    x2 = vec2[0]
    y2 = vec2[1]
    z2 = vec2[2]

    # calculate the distance
    dist = sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) + pow( z2 - z1, 2 ) )

    return dist



def compare_pose_energy_per_residue( sf, pose1, pose2, diff_cutoff = 1.0, detailed_analysis = False, compare_using_this_scoretype = None, weight = 1.0, verbose = False ):
    """
    Uses the ScoreFunction <sf> to compare the individual energy per residue of <pose1> and returns the count of residues that have a difference of greater than +/- 1 of <diff_cutoff> compared to the corresponding residue in <pose2>
    If you are collecting this information for use in a detailed analysis, consider setting <detailed_analysis> to True and this function will return the residue objects of <pose2> that satisfied the energy <diff_cutoff>
    If <compare_using_this_scoretype> is not already a ScoreType that is part of the <sf>, this function will add it to <sf> with a weight of <weight>
    :param sf: ScoreFunction
    :param pose1: Pose( original pose )
    :param pose2: Pose( new pose )
    :param diff_cutoff: int( or float( energy cutoff in REUs ). Positive number. Default = 1.0 REU
    :param detailed_analysis: bool( changes return values if you want the residues that satisfied the energy <diff_cutoff>). Default = False
    :param compare_using_this_scoretype: str( what specific ScoreType do you want to use for comparison? ). Default = None = total_energy
    :param weight: int( or float( weight of the desired <str_scoretype> if it's not already a part of <sf> ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: detailed_analysis = False, int( count of residues that satisfied the energy <diff_cutoff> )
    :return: detailed_analysis = True, int( count of residues that satisfied the energy <diff_cutoff> ), list( residues that satisfied the energy <diff_cutoff> of type Residue() )
    """
    # ensure that a number was given as a cutoff value
    try:
        diff_cutoff = abs( diff_cutoff )
        neg_diff_cutoff = diff_cutoff * -1
    except TypeError:
        print
        print
        print "Did you mean to give me a number for the diff_cutoff? Exiting"
        sys.exit()
    except:
        print
        print
        print "Something I don't understand happened. Raising the error."
        raise

    # ensure that the poses have the same number of amino acids (sequence identity doesn't matter)
    if pose1.total_residue() != pose2.total_residue():
        print "Your poses don't have the same number of residues, this comparison function won't work for this case. Exiting"
        print pose1
        print "vs"
        print pose2
        sys.exit()

    # otherwise continue on
    else:
        if verbose:
            if compare_using_this_scoretype is not None:
                print "Comparing residue v. residue energy using the specific ScoreType", compare_using_this_scoretype, "with a difference in energy cutoff value of", diff_cutoff
            else:
                print "Comparing residue v. residue total energy with a difference in energy cutoff value of", diff_cutoff

        # getting score only to allow access to the Pose's energies object
        sf( pose1 )
        sf( pose2 )

        # instantiate list that will hold Residue objects of residues that satisfied the energy <diff_cutoff>
        pose2_residues_with_greater_than_cutoff = []

        # number of residues in <pose2> that satisfied the energy <diff_cutoff>
        count = 0
        for ii in range( 1, pose1.total_residue() + 1 ):
            # if a ScoreType was given to use for comparison, calculate that energy
            if compare_using_this_scoretype is not None:
                # get ScoreType energy of the two residues
                pose1_res_E = get_residue_score_by_scoretype( sf, compare_using_this_scoretype, ii, pose1, weight = weight )
                pose2_res_E = get_residue_score_by_scoretype( sf, compare_using_this_scoretype, ii, pose2, weight = weight )

            # else, calculate total_energy and use that for comparison
            else:
                # get the total energy of each residue
                pose1_res_E = pose1.energies().residue_total_energy( ii )
                pose2_res_E = pose2.energies().residue_total_energy( ii )

            # get the energy difference, whether from total_energy or from a specific ScoreType
            E_diff = pose2_res_E - pose1_res_E

            # if E_diff is outside the +/- diff_cutoff
            if not neg_diff_cutoff < E_diff < diff_cutoff:
                pose2_residues_with_greater_than_cutoff.append( pose2.residue( ii ) )
                count += 1
                ##print "Pose1", pose1_res_E, "vs Pose2", pose2_res_E, "at position", ii

        if detailed_analysis:
            return count, pose2_residues_with_greater_than_cutoff
        else:
            return count



def compare_these_poses_by_score( sf, pose1, pose2, compare_using_this_scoretype = None, verbose = False ):
    """
    Calculates the energy of <pose1> and compares it to the energy of <pose2>, returning the Pose with the lower score
    :param sf: ScoreFunction
    :param pose1: Pose( pose1 )
    :param pose2: Pose( pose2 )
    :param compare_using_this_scoretype: str( what specific ScoreType do you want to use for comparison? ). Default = None = total_energy
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: the Pose with the lower energy of the two inputs
    """
    # return the better scored structure based on total_energy
    if compare_using_this_scoretype is None:
        if verbose:
            print "Returning the better scored Pose based on total energy"

        # return whichever Pose has the lower total energy score
        if sf( pose1 ) < sf( pose2 ):
            return pose1
        else:
            return pose2

    # otherwise, score by the given ScoreType and return the better scored Pose
    else:
        if verbose:
            print "Returning the better scored Pose based on", compare_using_this_scoretype, "energy"

        # get the specific ScoreType energy for the Poses
        pose1_E = get_score_by_scoretype( sf, compare_using_this_scoretype, pose1 )
        pose2_E = get_score_by_scoretype( sf, compare_using_this_scoretype, pose2 )

        # check to make sure neither energy is None (which means there was an error). If it is, exit this function
        if pose1_E is None or pose1_E is None:
            # don't need to print an error message because get_score_by_scoretype already does that
            sys.exit()

        # else return the better scored pose based on ScoreType
        else:
            if pose1_E < pose2_E:
                return pose1
            else:
                return pose2



# TODO-see how many inner and outer trials are actually necessary to hit a decent minimum
def get_best_structure_based_on_score( sf, pose, apply_sf_sugar_constraints = True, outer_trials = 3, inner_trials = 3, compare_using_this_scoretype = None, dump_best_pose = False, dump_pose_name = None, dump_dir = None, verbose = False, pmm = None ):
    """
    Packs and minimizes a <pose> <inner_trials> times, then uses the best Pose based on total score to pack and minimize again <outer_trials> times using the ScoreFunction <sf>
    :param sf: ScoreFunction
    :param pose: Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param outer_trials: int( number of times to run <inner_trials> ). Default = 3
    :param inner_trials: int( number of times to pack and minimize before calling that the temporary "best" structure ). Default = 3
    :param compare_using_this_scoretype: str( what specific ScoreType do you want to use for comparison? ). Default = None = total_energy
    :param dump_best_pose: bool( after finding the lowest energy Pose, dump the structure into the current directory (or to <dump_dir>). Default = False
    :param dump_pose_name: bool( filename of the Pose to be dumped ). Default = None (ie. "Best_" + current Pose name)
    :param dump_dir: str( path/to/dir/where/pose/will/be/dumped ). Default = None (ie. current working directory)
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: Pose( the Pose with the lowest total score after the trials of packing and minimization )
    """
    # inform user of options depending on how much they want to hear
    if verbose:
        print "\nUsing", outer_trials, "outer loops and", inner_trials, "inner loops to find a low-energy structure based on total score for pose", pose.pdb_info().name().split( '/' )[-1]
    else:
        print "Getting a base structure for your pose using packing and minimization..."
        
    # apply the PyMOL_Mover if passed
    if pmm is not None:
        # used in the function elsewhere to ensure the program doesn't use a broken mover
        pmm_worked = False
        try:
            if verbose:
                print "Watching your pose now"
            pmm.apply( pose )
            pmm_worked = True
        except:
            print "Something was wrong with your PyMOL_Mover -- continuing without watching"
            pass

    # apply sugar branch point constraints to sf, if desired
    if apply_sf_sugar_constraints:
        apply_sugar_constraints_to_sf( sf, pose )

    # make a dummy best pose - will be replaced as a new best pose is found
    best_pose = Pose()
    best_pose.assign( pose )

    # pack and minimize a temp pose
    for ii in range( outer_trials ):
        pose.assign( best_pose )
        
        # print current outer trial
        if verbose:
            print "Outer trial", ii + 1
            
        # watch pose in PyMOL
        if pmm is not None and pmm_worked:
            pmm.apply( pose )

        # now cycle through the inner pack_min loops
        for jj in range( inner_trials ):
            # print current inner trial
            if verbose:
                print "  Inner trial", jj + 1
                
            # this makes a temporary Pose which eventually gets replaced with an updated best Pose, make the minimization trajectory a little more directed
            temp_pose = Pose()
            temp_pose.assign( pose )
            if pmm is not None and pmm_worked:
                pmm.apply( temp_pose )
                
            # pack and minimize
            temp_pose.assign( do_pack_min( sf, temp_pose, apply_sf_sugar_constraints = apply_sf_sugar_constraints ) )
            if pmm is not None and pmm_worked:
                pmm.apply( temp_pose )
            best_pose.assign( compare_these_poses_by_score( sf, best_pose, temp_pose, compare_using_this_scoretype, verbose = verbose ) )


    # dump the Pose, if desired
    if dump_best_pose:
        # set up the correct Pose filename if none given
        if dump_pose_name is None:
            # default Pose name
            dump_pose_name = "Best_" + pose.pdb_info().name()

        # set up the correct dumping directory
        if dump_dir is None:
            # dump directory is the current working directory since it was not specified at input
            dump_dir = os.getcwd()

            if verbose:
                print "Dumping the Pose now named", dump_pose_name, "into", dump_dir

        else:
            # check if given directory is a valid path, else, dump to current working directory instead and inform user
            if not os.path.isdir( dump_dir ):
                print dump_dir, "is not a valid directory. Dumping the Pose now named", dump_pose_name, "to current working directory (", dump_dir, ") instead"
                dump_dir = os.getcwd()

            # the specified dump directory is valid, so the Pose will be dumped there
            else:
                if verbose:
                    print "Dumping the Pose now named", dump_pose_name, "into the directory", dump_dir

        # dump the Pose by navigating to dump directory
        cur_dir = os.getcwd()
        os.chdir( dump_dir )

        # dump the pose and return to the current working directory
        best_pose.dump_pdb( dump_pose_name )
        os.chdir( cur_dir )

    return best_pose



####################
#### LOOPS CODE ####
####################


def make_loop( start, stop, cut = None ):
    """
    Makes a loop starting at <start>, ending at <end>, and with a cutpoint at either <cut>, or the midpoint between <start> and <stop>
    :param start: int( start of loop )
    :param stop: int( end of loop )
    :param cut: int( cutpoint ). Default = None (ie. the midpoint between <start> and <stop>)
    :return: a Loop object
    """
    # if a cut was specified on input
    if cut is not None:
        loop = Loop( start, stop, cut )

    # otherwise, make the cut the midpoint between start and stop
    else:
        loop = Loop( start, stop, int( ( start + stop ) / 2 ) )

    return loop



def make_loops( starts, stops, cuts = None ):
    """
    Makes a Loops object out of a list of <starts> and <stops> with an optional list of <cuts>
    :param starts: list( int( start of loops ) )
    :param stops: list( int( end of loops ) )
    :param cuts: list( int( cutpoints ) ). Default = None (ie. the midpoint between <start> and <stop>)
    :return: a Loops object
    """
    # check that the arguments passed are lists - need lists for code to work
    if not isinstance( starts, list ):
        print
        print starts, "isn't a list. I need a list of int( start points ). Exiting"
        sys.exit()
    if not isinstance( stops, list ):
        print
        print stops, "isn't a list. I need a list of int( stop points ). Exiting"
        sys.exit()
    if cuts is not None:
        if not isinstance( cuts, list ):
            print
            print cuts, "isn't a list. I need a list of int( cutpoints ). Exiting"
            sys.exit()

    # check to make sure there are as many starts as stops and as cuts (if passed)
    if cuts is None:
        if len( starts ) != len( stops ):
            print
            print "You didn't give me the right amount of starts and ends. Exiting"
            sys.exit()
    else:
        if len( starts ) != len( stops ) != len( cuts ):
            print
            print "You didn't give me the right amount of starts/ends/cuts. Exiting"
            sys.exit()

    # now make the Loops object
    loops = Loops()
    for ii in range( len( starts ) ):
        # if cut points were specified at input, use those
        if cuts is not None:
            loop = make_loop( starts[ii], stops[ii], cuts[ii] )

        # otherwise, make the cutpoints the midpoint between the starts and stops
        else:
            loop = make_loop( starts[ii], stops[ii] )
        loops.add_loop( loop )

    return loops



def make_loops_from_file( loops_file, pose, PDB_numbering = False, anchor_loops = False, verbose = False ):
    """
    Reads the <loops_file> and creates the corresponding <loops>. Returns None if there was an error.
    IMPORTANT: can go from PDB --> Pose numbering ONLY IF in your LOOPs file you include the chain LETTER after the stop or cut value. Set PDB_numbering to True for this feature
    :param loops_file: str( /path/to/LOOPs/file )
    :param pose: Pose
    :param PDB_numbering: bool( is this file numbered using PDB numbering? ) Default = False
    :param anchor_loops: bool( do you want to add +2-residue anchors to the ends of your loops? ). Default = True
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: a Loops object
    """
    # makes the Loops object that will hold all the Loops created
    loops = Loops()

    try:
        # make sure the <loops_file> is a valid file and can be opened
        with open( loops_file, 'rb' ) as fh:
            lines = fh.readlines()
    except IOError:
        # if it's not a valid path, print an error message and exit
        print
        print
        print loops_file, "does not appear to be a valid path. Exiting"
        sys.exit()
    except:
        # if something else happened, print an error message and exit
        print
        print
        print "Something I don't understand happened. Raising the error."
        raise

    if verbose:
        if PDB_numbering:
            print "Making", len( lines ), "Loops from the given LOOPs file and changing from PDB numbering to Pose numbering\n"
        else:
            print "Making", len( lines ), "Loops from the given LOOPs file\n"

    # create the Loops from the file
    for line in lines:
        if line != '':
            # strip out new line characters and white space  -  used to make string parsing more efficient for Loop creation
            line = line.rstrip( '\n' )
            white_space_columns = line.split( ' ' )

            columns = []
            # splits up the line by columns, ignoring anything that is a Blank Space (Taylor Swift approves!)
            for ii in white_space_columns:
                if ii != '':
                    columns.append( ii )

            # the start and stop point should always be at [1] and [2]  -  "EDGE" is always at [0]
            start = int( columns[1] )
            stop = int( columns[2] )
            
            # if anchor_loops is True, add subtract 2 from start and add 2 to stop
            if anchor_loops:
                start -= 2
                stop += 2

            # the first try and the following three excepts are relevant to finding or making the cutpoint
            # in the try case, we are trying to see if an integer cutpoint is specified at columns[3]
            try:
                # an integer cutpoint should be here
                cut = int( columns[3] )

                # check to make sure cut is between the given start and stop, otherwise print an error message and exit
                if not start < cut < stop:
##                if not ( start < cut < stop or stop < cut < start ):
                    print
                    print "The cut you gave me,", cut, ",is not between your start and stop points:", start, "and", stop, " -  Check your file. Exiting"
                    sys.exit()

            # in this case, the cutpoint doesn't exist, but the chain ID does  -  so make cutpoint instead
            # trying to make an integer out of a string results in a ValueError, so this will work for this case, but not for making a string out of an integer
            except ValueError:
                # no cut given, so make one by getting the midpoint of the start and stop points
                cut = int( ( start + stop ) / 2 )

                # ensure a chain ID does indeed exist here by checking if it's a letter
                chain = str( columns[3] )
                if chain not in all_letters_list:
                    print
                    print
                    print "I'm not sure what this is:", "'%s'" %chain, " -  I was expecting a letter from A through Z. Check your file. sys.exit()"
                    sys.exit()

            # in this case, the cutpoint doesn't exist, and there is no chain ID there either  -  so make the cutpoint, or exit if there was supposed to be a chain ID (if <PDB_numbering> is True)
            except IndexError:
                # no cut given, so make one by getting the midpoint of the start and stop points
                cut = int( ( start + stop ) / 2 )

                # chain ID doesn't exist in this case, exit if user set <PDB_numbering> to True. They need to have added a chain here
                if PDB_numbering:
                    print
                    print
                    print "You set PDB_numbering to True, but did not add a chain ID after the stop residue. Exiting"
                    sys.exit()

            # otherwise, something weird happened
            except:
                print
                print
                print "Something I don't understand happened. Raising the error."
                raise

############

            # now, look for a chain ID, if <PDB_numbering> is True (otherwise, it wouldn't be there anyway)
            if PDB_numbering:
                # try to see if anything exists at columns[4]  -  don't need to check columns[3] because that was done above
                try:
                    chain = columns[4]

                # if there is nothing there, exit because <PDB_numbering> is True, thus we need a chain ID
                except IndexError:
                    print
                    print
                    print "You set PDB_numbering to True, but did not add a chain ID after the cut residue. Exiting"
                    sys.exit()

                # if something weird happened, exit
                except:
                    print
                    print
                    print "Something I don't understand happened. Raising the error."
                    raise

                # try except won't work for there being a ValueError because it is allowed to make a string out of an integer, so we need another way to ensure that the chain ID given is a valid letter
                # check to make sure that columns[4] is a string  -  exit if it is an integer or anything else
                if isinstance( chain, str ):
                    # check to make sure that the chain is actually only a letter
                    if chain not in all_letters_list:
                        print
                        print "I'm not sure what this is:", "'%s'" %chain, " -  I was expecting a letter from A through Z. Check your file. Exiting"
                        sys.exit()

                # if it's an integer, exit since we are expecting a string
                elif isinstance( chain, int ):
                    print
                    print "I was not expecting an integer here, I need a string for the chain ID. Exiting"
                    sys.exit()

                # otherwise, I'm not sure what is in columns[4], so exit
                else:
                    print
                    print "I'm not sure what this is", columns[4], "I need a string for the chain ID. Exiting"
                    sys.exit()


            # if <PDB_numbering> is True, go from PDB to Pose numbering using the positions and chain IDs from file
            if PDB_numbering:
                start = pose.pdb_info().pdb2pose( chain, start )
                stop = pose.pdb_info().pdb2pose( chain, stop )
                cut = pose.pdb_info().pdb2pose( chain, cut )

            # add each Loop to the Loops object
            loop = Loop( start, stop, cut )
            loops.add_loop( loop )

    # print loop to user if verbose is True
    if verbose:
        print "Returning the following Loops object\n\t", loops
        
    return loops



def add_cutpoint_variants( loops, pose ):
    """
    Adds cutpoint_variants to the <pose> for each loop given in <loops>
    :param loops: a Loops object (can be one loop)
    :param pose: Pose
    :return: Pose( <pose> with <loops> that are cutpoint_variants )
    """

    # for every Loop in the Loops object, add it as a cutpoint variant to the <pose>
    for loop_num in range( 1, loops.num_loop() + 1 ):
        add_single_cutpoint_variant( pose, loops[ loop_num ] )

    return pose



########################
#### FOLD TREE CODE ####
########################

def foldtree_to_string_and_stripped( pose, verbose = False ):
    """
    Turns the FoldTree of <pose> into a list of strings only containing the relevant numbers to each EDGE
    Removes the words FOLD_TREE and EDGE and strips whitespace and splits the rest into chunks
    Example output:
    :param pose: Pose
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: list( relevant information from each EDGE)
    """
    if verbose:
        print "Turning the Pose's FoldTree into a string by splitting on 'EDGE' and stripping off whitespace"

    # get string version of FoldTree
    ft = pose.fold_tree().to_string()
    ft = ft.rstrip()

    # split FoldTree on "EDGE"
    ft_split_on_edge = ft.split( "EDGE" )

    # split on white space to get access to numbers and remove the FOLD_TREE tag at the beginning
    ft_split_on_white_space = []
    for edge in ft_split_on_edge:
        # strips off white space from the ends
        edge = edge.strip()
        # ignores the FOLD_TREE chunk
        if edge != "FOLD_TREE":
            # ignores any blank spaces
            if edge != '':
                ft_split_on_white_space.append( edge )

    return ft_split_on_white_space



def get_JUMP_NUM_from_seq_pos( seq_pos, pose, downstream = False, verbose = False ):
    """
    Return the associated up- or downstream jump number given a specific <seq_pos> of a <pose>
    Returns the upstream jump number by default, set <downstream> to True to retrieve upstream jump number of <seq_pos>
    :param seq_pos: int( sequence position of residue in question )
    :param pose: Pose
    :param downstream: bool( returns upstream residue by default, set <downstream> to True to return downstream jump num )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: int( up- or downstream jump number )
    """
    # first check to see that the <seq_pos> entered is a valid position
    if seq_pos <= 0 or seq_pos > pose.n_residue():
        print
        print "Invalid sequence position given:", seq_pos, "Exiting"
        sys.exit()

    if verbose:
        if downstream:
            print "Given the residue sequence number", seq_pos, "this will find the closest downstream Jump number"
        else:
            print "Given the residue sequence number", seq_pos, "this will find the closest upstream Jump number"

    # retrieve both upstream and downstream jump numbers relative to <seq_pos>
    upstream_jump = 0
    downstream_jump = 0
    last_jump_seen = 0

    # get string version without whitespace of FoldTree EDGE's
    ft_split_on_white_space = foldtree_to_string_and_stripped( pose )

    # split up and use each EDGE from the FoldTree
    for edge in ft_split_on_white_space:
        # split on whitespace to get only numbers
        edge_split = edge.split( ' ' )
        edge_split_on_white = []

        for e in edge_split:
            # ignore resulting white space and blank spaces
            if e != ' ':
                if e != '':
                    # append only the relevant numbers for this EDGE
                    edge_split_on_white.append( e )

        # get relevant EDGE numbers  -  not using start because the end residue is more relevant to use when looking for the closest Jump number
        end_res = int( edge_split_on_white[1] )
        edge_type = int( edge_split_on_white[2] )

        # if this is a jump EDGE  (if edge_type is a positive integer)
        if edge_type > 0:
            # store it as the last_jump_seen to be used in case the residue in question is beyond the last Jump
            last_jump_seen = edge_type

            # if the <seq_pos> is within this Jump EDGE  (or at least less than the end of it)
            if seq_pos <= end_res:
                # store both the upstream and downstream jump
                upstream_jump = int( edge_type )
                downstream_jump = upstream_jump - 1

                # return upstream or downstream jump depending on user input when it is found  -  no need to go through all of the FoldTree
                if downstream:
                    return downstream_jump
                else:
                    return upstream_jump

    # if the last_jump_seen != 0 but upstream_jump == 0, it means there are jumps, but the seq_pos is beyond the last j=Jump. So set the Jump to last_jump_seen
    if last_jump_seen > 0 and upstream_jump == 0:
        upstream_jump = last_jump_seen
        downstream_jump = last_jump_seen - 1

    # return upstream or downstream jump depending on user input
    if downstream:
        return downstream_jump
    else:
        return upstream_jump



def get_seq_pos_from_JUMP_NUM( JUMP_NUM, pose, start = False, verbose = False):
    """
    Return the associated start or stop sequence position corresponding to the <JUMP_NUM> in a <pose>
    By default, returns the stop position. Set <start> to True if you want the start position
    :param JUMP_NUM: int( Jump number )
    :param pose: Pose
    :param start: bool( do you want the start sequence position? Please set to True). Default = False (ie. returns stop position )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: int( sequence position of the start or stop residue )
    """
    # check to make sure a valid <JUMP_NUM> was passed
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print "Invalid jump number given:", JUMP_NUM, "Exiting"
        sys.exit()

    if verbose:
        if start:
            print "Getting the start residue closest to the Jump number", JUMP_NUM
        else:
            print "Getting the stop residue closest to the Jump number", JUMP_NUM

    # get string version without whitespace of FoldTree EDGE's
    ft_split_on_white_space = foldtree_to_string_and_stripped( pose )

    # split up and use each EDGE from the FoldTree
    for edge in ft_split_on_white_space:
        # split on whitespace to get only numbers
        edge_split = edge.split( ' ' )
        edge_split_on_white = []

        for e in edge_split:
            # ignore resulting white space and blank spaces
            if e != ' ':
                if e != '':
                    # append only the relevant numbers for this EDGE
                    edge_split_on_white.append( e )


        # find the start or stop residue of the given <JUMP_NUM>
        start_pos = int( edge_split_on_white[0] )
        stop_pos = int( edge_split_on_white[1] )
        edge_type = int( edge_split_on_white[2] )

        # if the edge_type is the <JUMP_NUM>
        if edge_type == JUMP_NUM:
            # return start or stop position, given user input
            if start:
                return start_pos
            else:
                return stop_pos



# LOOPs file --> appropriate fold tree
def hard_code_fold_tree( orig_loops, pose, verbose = False ):
    ft = FoldTree()
    
    # make a copy of the original loops as to not alter them for later use
    loops = Loops()
    for ii in range( 1, orig_loops.num_loop() + 1 ):
        loops.add_loop( orig_loops[ii].start(), orig_loops[ii].stop(), orig_loops[ii].cut() )
    if verbose:
        print "Hard coding FoldTree with the following Loops object:\n\t", loops
    
    # attach the Loops object to the <pose> to be used when tearing down the FoldTree with restore_original_fold_tree
    pose.loops = Loops()
    for ii in range( 1, loops.num_loop() + 1 ):
        pose.loops.add_loop( loops[ii] )
    
    # print original FoldTree, if verbose
    if verbose:
        print "Original FoldTree:\n\t", pose.fold_tree()
    
    # parse through the loops first
    starts_from_loops = []
    for ii in range( 1, loops.num_loop() + 1 ):
        starts_from_loops.append( int( loops[ii].start() ) )  # used for indexing
        
    hold_ft = []
    break_out = False
    JUMP_NUM = 1
    
    orig_ft = pose.fold_tree().to_string().split( "EDGE" )
    hold_ft = []
    for edge in orig_ft:
        edge = edge.strip()
        if edge != "FOLD_TREE":
            if edge != '':
                hold_ft.append( edge )
                        
    # outer loop should be the loop_starts in case there are two loops within one EDGE
    for edge in hold_ft:
        # split on whitespace
        temp_columns = edge.split( ' ' )
        columns = []
        for ii in temp_columns:
            if ii != '':
                columns.append( ii )
                    
        # check to see that this is a -1 EDGE ( this is the only EDGE I'd be altering anyway )
        if int( columns[2] ) == -1:
            loops_to_add = []
            for loop_start in starts_from_loops:
                # second number will be the end of the edge or jump, check to see that the new loop start isn't smaller than this sequence position
                if loop_start < int( columns[1] ):
                    loops_to_add.append( loop_start )
                        
            if len( loops_to_add ) != 0:
                if len( loops_to_add ) == 1:
                    for loop_start in loops_to_add:
                        # split this EDGE up to make room for the loop
                        # get the loop's index number within the loops object to get it's end and cut info
                        loop_index_number = starts_from_loops.index( loop_start ) + 1  # Loops are indexed starting at 1
                        loop = loops[ loop_index_number ]
                        
                        # add the edge from the original start to this loop's start point
                        ft.add_edge( int( columns[0] ), loop.start(), -1 )
                        
                        # add the edge from the start of the loop to its cutpoint
                        ft.add_edge( loop.start(), loop.cut(), -1 )
                        
                        # add the jump from the start to the end of the loop
                        ft.add_edge( loop.start(), loop.stop(), JUMP_NUM )
                        JUMP_NUM += 1
                        
                        # add the edge from the end of the loop to its upper cutpoint
                        ft.add_edge( loop.stop(), loop.cut() + 1, -1 )
                        
                        # finally, go from end of loop to the original end point
                        ft.add_edge( loop.stop(), int( columns[1] ), -1 )
                        
                        # remove the loop from the loops list
                        starts_from_loops.remove( loop_start )
                        loops.delete_loop( loop.start(), loop.stop() )
                else:
                    ii = 0
                    jj = 1
                    while ( ii + 1 ) != len( loops_to_add ):
                        loop_index_number = starts_from_loops.index( loops_to_add[ii] ) + 1  # Loops are indexed starting at 1
                        loop1 = loops[ loop_index_number ]
                        
                        loop_index_number = starts_from_loops.index( loops_to_add[jj] ) + 1  # Loops are indexed starting at 1
                        loop2 = loops[ loop_index_number ]
                        
                        if ii == 0:
                            # the start is from the current edge start to the new loop start
                            ft.add_edge( int( columns[0] ), loop1.start(), -1 )
                            
                            # add the edge from the start of the loop to its cutpoint
                            ft.add_edge( loop1.start(), loop1.cut(), -1 )
                            
                            # add the jump from the start to the end of the loop
                            ft.add_edge( loop1.start(), loop1.stop(), JUMP_NUM )
                            JUMP_NUM += 1
                            
                            # add the edge from the end of the loop to its upper cutpoint
                            ft.add_edge( loop1.stop(), loop1.cut() + 1, -1 )
                            
                            # finally, go from end of loop to the start of the next loop
                            ft.add_edge( loop1.stop(), loop2.start(), -1 )
                            
                            # remove the loop from the loops list
                            if ii != ( len( loops_to_add ) - 1 ):
                                starts_from_loops.remove( loop1.start() )
                                loops.delete_loop( loop1.start(), loop1.stop() )
                            
                            # up the loop indices
                            ii += 1
                            jj += 1
                                
                        else:
                            # add the edge from the start of the loop to its cutpoint
                            ft.add_edge( loop1.start(), loop1.cut(), -1 )
                            
                            # add the jump from the start to the end of the loop
                            ft.add_edge( loop1.start(), loop1.stop(), JUMP_NUM )
                            JUMP_NUM += 1
                            
                            # add the edge from the end of the loop to its upper cutpoint
                            ft.add_edge( loop1.stop(), loop1.cut() + 1, -1 )
                            
                            # finally, go from end of loop to the start of the next loop
                            ft.add_edge( loop1.stop(), loop2.start(), -1 )
                            
                            # remove the loop from the loops list
                            if ii != ( len( loops_to_add ) - 2 ):
                                starts_from_loops.remove( loop2.start() )
                                loops.delete_loop( loop2.start(), loop2.stop() )
                                
                            # up the loop indices
                            ii += 1
                            jj += 1
                            
                    # add the final loop
                    loop_start = loops_to_add[-1]
                    loop_index_number = starts_from_loops.index( loop_start ) + 1  # Loops are indexed starting at 1
                    loop2 = loops[ loop_index_number ]
                    
                    # add the edge from the start of the loop to its cutpoint
                    ft.add_edge( loop2.start(), loop2.cut(), -1 )
                    
                    # add the jump from the start to the end of the loop
                    ft.add_edge( loop2.start(), loop2.stop(), JUMP_NUM )
                    JUMP_NUM += 1
                    
                    # add the edge from the end of the loop to its upper cutpoint
                    ft.add_edge( loop2.stop(), loop2.cut() + 1, -1 )
                    
                    # finally, go from end of loop to the original end point
                    ft.add_edge( loop2.stop(), int( columns[1] ), -1 )
                    
                    # remove the loop from the loops list
                    starts_from_loops.remove( loop2.start() )
                    loops.delete_loop( loop2.start(), loop2.stop() )
                    
            else:
                ft.add_edge( int( columns[0] ), int( columns[1] ), -1 )
                    
        if int( columns[2] ) == -2:
            # chemical edge ( sugar )
            ft.add_edge( int( columns[0] ), int( columns[1] ), columns[3], columns[4] )
            
        if int( columns[2] ) > 0:
            # a jump
            ft.add_edge( int( columns[0] ), int( columns[1] ), JUMP_NUM )
            JUMP_NUM += 1
            
    if ft.check_fold_tree():
        if verbose:
            print "New FoldTree:\n\t", ft
        return ft
    else:
        print "There was a problem hard coding your fold tree"
        sys.exit()



def setup_new_fold_tree( loops_file, pose, PDB_numbering = False, anchor_loops = True, add_cutpoints = True, verbose = False ):
    """
    Creates a FoldTree given a <loops_file> and the <pose>'s current FoldTree and then gives it to <pose>
    If the Loops result in a invalid FoldTree, the program will result in an error and exit
    :param loops_file: str( /path/to/valid/LOOPs/file )
    :param pose: Pose
    :param PDB_numbering: bool( is this file numbered using PDB numbering? ) Default = False
    :param anchor_loops: bool( do you want to add +2-residue anchors to the ends of your loops? ). Default = True
    :param add_cutpoints bool( do you want to add cutpoints in the FoldTree for each of the loop's cutpoints? ). Default = True
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: Pose with the new FoldTree
    """
    # attach the Loops object to the <pose> to be used when tearing down the FoldTree with restore_original_fold_tree
    pose.loops = Loops()
    
    # check to see if user passed a Loop(s) object or the file name of a Loop(s) object
    if isinstance( loops_file, str ):
        # add the name of the loops file to the pose to be used when returning it to its normal FoldTree
        pose.loops_file = loops_file
        
        # make the Loops from LOOPs file
        loops = make_loops_from_file( loops_file, pose, PDB_numbering = PDB_numbering, anchor_loops = anchor_loops, verbose = verbose )
        for ii in range( 1, loops.num_loop() + 1 ):
            pose.loops.add_loop( loops[ii] )

    else:
        print
        print "I'm not quite sure what you passed to me. I need the name of the file containing your Loops, or an actual Loop or Loops object. Exiting."
        sys.exit()
    
    # apply the new FoldTree by altering the current one to include the new Loops
    pose.fold_tree( hard_code_fold_tree( loops, pose, verbose = verbose ) )

    # add cutpoint variants for each new Loop
    if add_cutpoints:
        if verbose:
            print "Pose's original set of cutpoints:\n\t", pose.fold_tree().cutpoints()
        pose.assign( add_cutpoint_variants( loops, pose ) )
        if verbose:
            print "Pose's new set of cutpoints:\n\t", pose.fold_tree().cutpoints()
    
    return pose



def restore_original_fold_tree( pose, verbose = False ):
    """
    Returns the Pose's FoldTree back to its original version and removes the added cutpoint variants
    The Pose's original FoldTree was attached to itself during the function load_pose
    :param pose: Pose
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: Pose with original FoldTree
    """
    # check to make sure there's actually a reason to restore the original FoldTree
    # if pose.loops_file is still None, that means it hasn't been changed  -  therefore nothing to do for this <pose>
    if pose.loops_file is None:
        print
        print "I see you haven't actually changed the fold tree. Is this true?"
        return pose

    # otherwise, restore it
    else:
        if verbose:
            print "Returning the Pose's FoldTree back to its original one by removing the cutpoint variants added and then applying the original FoldTree stored from from the load_pose function"

        ## IMPORTANT: have to get the loops again because the fold tree thinks the sugars are cutpoints - will mess it up
        loops = pose.loops
        if verbose:
            print "  Tearing down the following Loops:\n\t", loops

        # remove cutpoint variants by reading the LOOPs file attached to the <pose>
        for ii in range( 1, loops.size() + 1 ):
            # make the Loop from the Loops object
            loop = loops[ii]

            # get its lower and upper cutpoints
            lower_cut = loop.cut()
            upper_cut = lower_cut + 1

            # remove the variant types for the Loop from the LOOPs file
            remove_variant_type_from_pose_residue( pose, VariantType.CUTPOINT_LOWER, lower_cut )
            remove_variant_type_from_pose_residue( pose, VariantType.CUTPOINT_UPPER, upper_cut )

        # apply the original FoldTree back to the <pose>
        pose.fold_tree( pose.orig_fold_tree )
        
        if verbose:
            print "Restored Pose to following FoldTree:\n\t", pose.fold_tree()
            print "Pose now has the following cutpoints:\n\t", pose.fold_tree().cutpoints()

        # return the LOOPs file object in the <pose> back to None since it's now back to original FoldTree
        pose.loops_file = None
        pose.loops = None

        return pose



#####################################
#### MUTATIONAL WORKER FUNCTIONS ####
#####################################

def make_mutation_packer_task( amino_acid, seq_pos, sf, pose, pack_radius = PACK_RADIUS ):
    """
    Returns a packer task that can handle a <pose> with a SINGLE mutation at <seq_pos>
    :param amino_acid: str( ONE LETTER amino acid code for the SINGLE mutated residue )
    :param seq_pos: int( sequence position of the mutated residue )
    :param sf: ScoreFunction
    :param pose: Pose
    :param pack_radius: int( or float( the distance in Angstroms around the mutated residue you want to be packed ). Default = PACK_RADIUS = 20
    :return: packer_task ready to pack a Pose with a SINGLE mutation
    """
    # create a packer task handling a SINGLE mutated residue
    # tell packer task that only the current <amino_acid> should be allowed to be modified
    aa_bool = vector1_bool()
    mutant_amino_acid = aa_from_oneletter_code( amino_acid )
    for ii in range( 1, 20 + 1 ):
        aa_bool.append( ii == mutant_amino_acid )

    # get center of mutated residue to be used in distance calculation for pack radius
    center = pose.residue( seq_pos ).nbr_atom_xyz()

    # make the packer task
    task = standard_packer_task( pose )
    task.or_include_current( True )
    task.nonconst_residue_task( seq_pos ).restrict_absent_canonical_aas( aa_bool )
    task.restrict_to_repacking()

    # turn off repacking for residues not in the pack radius
    for ii in range( 1, pose.total_residue() + 1 ):
        if center.distance( pose.residue( ii ).nbr_atom_xyz() ) > pack_radius:
            task.nonconst_residue_task( ii ).prevent_repacking()
    pack_rotamers_mover = RotamerTrialsMover( sf, task )

    return pack_rotamers_mover



def do_mutation_pack( seq_pos, amino_acid, sf, mutated_pose, pack_radius = PACK_RADIUS ):
    """
    Returns a packed and minimized <pose> that has a SINGLE mutation of <amino_acid> at <seq_pos>
    :param seq_pos: int( sequence position of the mutated residue _
    :param amino_acid: str( ONE LETTER amino acid code for the SINGLE mutated residue )
    :param sf: ScoreFunction
    :param mutated_pose: Pose with the SINGLE mutation
    :param pack_radius: int( or float( the distance in Angstroms around the mutated residue you want to be packed ). Default = PACK_RADIUS = 20
    :return: Pose packed around the SINGLE mutation
    """
    pose = Pose()
    pose.assign( mutated_pose )
    
    # pack the <mutated_pose>
    pack_rotamers_mover = make_mutation_packer_task( amino_acid, seq_pos, sf, pose, pack_radius )
    pack_rotamers_mover.apply( pose )

    return pose



###################################
#### MAIN MUTANT POSE CREATION ####
###################################

def get_best_mutant_of_20( seq_pos, sf, pose, apply_sf_sugar_constraints = True, rounds = 1, pack_radius = PACK_RADIUS ):
    """
    For a single position <seq_pos>, mutated to all 20 amino acids, pack and minimize over the given number of <rounds>, and return the best mutant pose
    :param seq_pos: int( the residue position that will be mutated )
    :param sf: ScoreFunction
    :param pose: Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param rounds: int( the number of times to pack and minimize the <pose> after the mutation ). Default = 1
    :param pack_radius: int( or float( the distance in Angstroms around the mutated residue you want to be packed ). Default = PACK_RADIUS = 20
    :return: the mutated Pose of the lowest total score out of the twenty mutations
    """
    # apply sugar branch point constraints to sf, if desired
    if apply_sf_sugar_constraints:
        apply_sugar_constraints_to_sf( sf, pose )

    # talk to user
    orig_amino_acid = pose.residue( seq_pos ).name1()
    print "Mutating", seq_pos, orig_amino_acid

    for amino_acid in AA_list:
        print "Now mutating to", amino_acid, "..."

        # mutate to a new residue
        mutant_pose = mutate_residue( pose, seq_pos, amino_acid )

        # do one round of packmin to get rid of clashes and prepare for comparison
        mutant_pose = do_mutation_pack( seq_pos, amino_acid, sf, mutant_pose, pack_radius )

        # do mutation packmin as many times as specified, saving the best scored one
        for ii in range( rounds ):
            # assign a temp pose that will be separately packmined
            temp_pose = Pose()
            temp_pose.assign( mutant_pose )

            # packmin the temp pose
            temp_pose = do_mutation_pack( seq_pos, amino_acid, sf, temp_pose )

            # (re)-assign best pose to whichever packmin pose is better (temp or orig)
            best_mutant = compare_pose_energy_per_residue( sf, mutant_pose, temp_pose )

    return best_mutant



def make_all_mutations( sf, orig_pose_file, mutant_list_file, pack_around_mut = True, dump_pose = False, dump_dir = None ):
    # ensure the mutation list filename is accurate and can be opened
    try:
        mutant_list = []
        f = open( mutant_list_file, 'rb' )
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            if line != '' and line[0] != '#':
                mutant_list.append( line )
    except IOError:
        print
        print
        print mutant_list_file, "didn't work. Check your input"
        sys.exit()
    except:
        print
        print
        print "Unexpected error"
        raise
    
    # ensure the validity of the orig_pose_file path
    if os.path.isfile( orig_pose_file ):
        orig_pose = Pose()
        orig_pose.assign( load_pose( orig_pose_file ) )
    else:
        print orig_pose_file, "is not a valid file path"
        print "Exiting"
        sys.exit()
    
    # if no dump directory was given, use the current working directory
    if dump_dir is None:
        dump_dir = os.getcwd()

    # if one was given, ensure it is valid; otherwise, create that directory
    else:
        if not os.path.isdir( dump_dir ):
            os.mkdir( dump_dir )

    # for each mutation in list, run the mutation function given sym or asym designation
    for mut_line_full in mutant_list:
        mut_line = mut_line_full.split( ' ' )
        mut = mut_line[ 0 ]
        symmetry = mut_line[ 1 ]
        if symmetry == '' or symmetry == "sym":
            print "Symmetrical", mut
            make_my_new_symmetric_antibody( mut, sf, orig_pose, 
                                            pack_around_mut = pack_around_mut, 
                                            dump_pose = dump_pose, 
                                            dump_dir = dump_dir )
        elif symmetry == "asym":
            print "Asymmetrical", mut
            make_my_new_asymmetric_antibody( mut, sf, orig_pose, 
                                             pack_around_mut = pack_around_mut, 
                                             dump_pose = dump_pose, 
                                             dump_dir = dump_dir )
        else:
            print "'%s'" %symmetry, "isn't a valid a symmetrical designation. Please put 'sym' or 'asym'"
            print "Exiting"
            sys.exit()
    
    return True



def make_my_new_symmetric_antibody( mutation_string, sf, input_pose, apply_sf_sugar_constraints = True, pack_around_mut = True, dump_pose = False, dump_dir = None, verbose = False ):
    """
    Turn <input_pose> into a mutant with mutations specified by <mutation_string>. If specified, dump mutant into <dump_dir>
    IMPORTANT: Since the Fc region of antibodies are symmetric, this code automatically takes the mutation from chain A and mutates its equivalent in chain B
    Example format for mutation_string ( 1 mutation A245F. Multiple mutations A245F_L98K_S600Q )
    Mutant <pose> will not be packed and minimized as it cannot handle multiple mutations. Instead, take dumped Pose and do it then
    :param mutation_string: format str( <single letter code of original amino acid><PDB sequence position><single letter code of new amino acid> ). Split multiple mutations using '_'
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param pack_around_mut: bool( do you want to pack around the mutation(s) made? ). Default = True = Yes
    :param dump_pose: bool( do you want to dump your mutated pose? ). Default = False = No
    :param dump_dir: str( /path/to/dump/directory/for/pose ). Default = None = current working directory
    :return: newly mutated Pose
    """
    # no pack min to get best structure!! do that yourself with the dumped pose

    # print out the mutation string to be created
    if verbose:
        print "Adding the following mutations to your pose:", mutation_string
    
    # move the passed pose into a separate Pose object
    pose = Pose()
    pose.assign( input_pose )

    # apply sugar branch point constraints to sf, if desired
    if apply_sf_sugar_constraints:
        apply_sugar_constraints_to_sf( sf, pose )

    # get list of mutations to make
    mutations = mutation_string.split( '_' )

    # mutate all residues in chain A and chain B
    # mutate first residue in chain A and chain B so the rest can be a loop
    for mut in mutations:
        orig_amino_acid = mut[ 0 ]
        new_amino_acid = mut[ -1 ]

        # seq pos is always from 1 to one minus however many characters are in the string
        pdb_seq_pos = mut[ 1:( len( mut ) - 1 ) ]

        # get pose positions for both chain A
        pose_A_pos = pose.pdb_info().pdb2pose( 'A', int( pdb_seq_pos ) )
        pose_B_pos = pose.pdb_info().pdb2pose( 'B', int( pdb_seq_pos ) )

        # ensure that the original amino acid the user specified is actually there
        if pose.residue( pose_A_pos ).name1() != orig_amino_acid and pose.residue( pose_B_pos ).name1() != orig_amino_acid:
            print "Hold up! What you said was the original amino acid is actually incorrect!!"
            print "You told me there was originally a", orig_amino_acid, "at position", pdb_seq_pos, "but there actually was a", pose.residue( pose_A_pos ).name1(), ". Exiting."
            sys.exit()

        # make mutation with no pack or min (doing our own packing just to get rid of clashes for each point mutation)
        # for chain A mutation
        pose.assign( mutate_residue( pose, pose_A_pos, new_amino_acid ) )

        # for chain B mutation
        pose.assign( mutate_residue( pose, pose_B_pos, new_amino_acid ) )
        
        # pack around mutations, if desired
        if pack_around_mut:
            pose.assign( do_mutation_pack( pose_A_pos, new_amino_acid, sf, pose ) )
            pose.assign( do_mutation_pack( pose_B_pos, new_amino_acid, sf, pose ) )
        
    # dump the new pose if desired
    if dump_pose:
        pose.pdb_info().name( mutation_string )
        if dump_dir is None:
            filename = mutation_string + ".pdb"
            pose.dump_pdb( filename )
            return pose
        else:
            # add a slash to the end of the directory name if there wasn't already one there
            if not dump_dir.endswith( '/' ):
                dump_dir += '/'
            filename = dump_dir + mutation_string + ".pdb"
            pose.dump_pdb( filename )
            return pose
    
    return pose



def make_my_new_asymmetric_antibody( mutation_string, sf, input_pose, apply_sf_sugar_constraints = True, pack_around_mut = True, dump_pose = False, dump_dir = None, verbose = False ):
    """
    Turn <input_pose> into an asymmetric mutant with mutations specified by <mutation_string>. If specified, dump mutant into <dump_dir>
    Mutant <pose> will not be packed and minimized as it cannot handle multiple mutations. Instead, take dumped Pose and do it then
    :param mutation_string: format str( <single letter code of original amino acid><PDB sequence position><single letter code of new amino acid>_<chain id> ). Split multiple mutations using '+'
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param apply_sf_sugar_constraints: bool( add sugar bond anlge and distance constraints to the sf? ). Default = True
    :param pack_around_mut: bool( do you want to pack around the mutation(s) made? ). Default = True = Yes
    :param dump_pose: bool( do you want to dump your mutated pose? ). Default = False = No
    :param dump_dir: str( /path/to/dump/directory/for/pose ). Default = None = current working directory
    :return: newly mutated Pose
    """
    # no pack min to get best structure!! do that yourself with the dumped pose

    # print out the mutation string to be created
    if verbose:
        print "Adding the following mutations to your pose:", mutation_string
    
    # move the passed pose into a separate Pose object
    pose = Pose()
    pose.assign( input_pose )

    # apply sugar branch point constraints to sf, if desired
    if apply_sf_sugar_constraints:
        apply_sugar_constraints_to_sf( sf, pose )

    # get list of mutations to make
    mutations = mutation_string.split( '+' )

    # mutate all residues in chain A and chain B
    # mutate first residue in chain A and chain B so the rest can be a loop
    for mut in mutations:
        chain_id = mut.split( '_' )[ -1 ]
        mut_string = mut.split( '_' )[ 0 ]
        orig_amino_acid = mut_string[ 0 ]
        new_amino_acid = mut_string[ -1 ]

        # seq pos is always from 1 to one minus however many characters are in the mutation string
        pdb_seq_pos = mut_string[ 1:( len( mut_string ) - 1 ) ]

        # get pose position of residue
        pose_pos = pose.pdb_info().pdb2pose( chain_id, int( pdb_seq_pos ) )

        # ensure that the original amino acid the user specified is actually there
        if pose.residue( pose_pos ).name1() != orig_amino_acid:
            print "Hold up! What you said was the original amino acid is actually incorrect!!"
            print "You told me there was originally a", orig_amino_acid, "at position", pdb_seq_pos, "but there actually was a", pose.residue( pose_pos ).name1(), ". Exiting."
            sys.exit()

        # make the point mutation
        pose.assign( mutate_residue( pose, pose_pos, new_amino_acid ) )
        
        # pack the mutation, if desired
        if pack_around_mut:
            pose.assign( do_mutation_pack( pose_pos, new_amino_acid, sf, pose ) )
        
    # dump the new pose if desired
    if dump_pose:
        pose.pdb_info().name( mutation_string )
        if dump_dir is None:
            filename = mutation_string + ".pdb"
            pose.dump_pdb( filename )
            return pose
        else:
            # add a slash to the end of the directory name if there wasn't already one there
            if not dump_dir.endswith( '/' ):
                dump_dir += '/'
            filename = dump_dir + mutation_string + ".pdb"
            pose.dump_pdb( filename )
            return pose
        
    return pose



########################
#### DATA FUNCTIONS ####
########################

def get_contact_map( pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Returns a dictionary of each residue in <pose> that has a contact with another residue in the <pose> less than the given <cutoff> distance
    :param pose: Pose
    :param cutoff: int( or float( distance cutoff for what defines a contact). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: dict( key = residue num in pose, value = list of residue nums within <cutoff> distance )
    """
    if verbose:
        print "Getting the residue contact map for the Pose"

    # holds the resulting contact map
    return_dict = {}

    for seq_pos in range( 1, pose.n_residue() + 1 ):
        # holder for the residue numbers of residues within the given <cutoff> distance
        contacts = []

        # get the center for residue one
        center_1 = list( pose.residue( seq_pos ).nbr_atom_xyz() )

        # for every residue in the <pose>
        for residue in pose:
            seq_pos_2 = residue.seqpos()

            # if it's not the same residue (obviously they contact)
            if seq_pos != seq_pos_2:
                # get the center of the second residue
                center_2 = list( residue.nbr_atom_xyz() )

                # calculate the distance and store the residue number if it's below the given <cutoff>
                if calc_distance( center_1, center_2 ) < cutoff:
                    contacts.append( seq_pos_2 )

        # fill the return_dict with the contacts to that residue, only if there any
        if len( contacts ) != 0:
            return_dict[ seq_pos ] = contacts

    return return_dict



def get_contact_map_between_range1_range2( range1, range2, pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Returns a dictionary of each residue in <pose> that has a contact with another residue in the <pose> less than the given <cutoff> distance
    :param range1: list( residue pose numbers in range 1 )
    :param range2: list( residue pose numbers in range 2 )
    :param pose: Pose
    :param cutoff: int( or float( distance cutoff for what defines a contact). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: dict( key = residue num in pose, value = list of residue nums within <cutoff> distance )
    """
    if verbose:
        print "Getting the residue contact map for the Pose"

    # holds the resulting contact map
    return_dict = {}

    for seqpos_1 in range1:
        # holder for the residue numbers of residues within the given <cutoff> distance
        contacts = []

        # get the center for residue one
        center_1 = list( pose.residue( seqpos_1 ).nbr_atom_xyz() )

        # for every residue in the <pose>
        for seqpos_2 in range2:
            # get the center of the second residue
            center_2 = list( pose.residue( seqpos_2 ).nbr_atom_xyz() )

            # calculate the distance and store the residue number if it's below the given <cutoff>
            if calc_distance( center_1, center_2 ) < cutoff:
                contacts.append( seqpos_2 )

        # fill the return_dict with the contacts to that residue, only if there any
        if len( contacts ) != 0:
            return_dict[ seqpos_1 ] = contacts

    return return_dict



def calc_res_Fnat_recovered_between_range1_range2( decoy, decoy_r1, decoy_r2, native, cutoff = CUTOFF_DISTANCE, return_more_data = False, verbose = False ):
    """
    Uses a residue contact map to determine if residues in the <decoy> between range <decoy_r1> and <decoy_r2> match the native contacts made in the corresponding residue range using the <cutoff> distance
    :param decoy: decoy Pose
    :param decoy_r1: list( residue numbers on side 1 )
    :param decoy_r2: list( residue numbers on side 2 )
    :param native: native Pose
    :param cutoff: int( or float( distance cutoff for what defines a contact). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param return_more_data: bool( do you want to have this function return decoy and native number of contacts as well? ) Default is False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( residue Fnat recovered )
    :return: if return_more_data == True, int( decoy contacts ), int( native contacts ), float( residue Fnat recovered )
    """
    # get the corresponding native residue numbers in the decoy ranges
    native_r1 = []
    for resnum in decoy_r1:
        native_r1.append( decoy_to_native_res_map[ resnum ] )
    native_r2 = []
    for resnum in decoy_r2:
        native_r2.append( decoy_to_native_res_map[ resnum ] )
    
    # get the contact maps for the decoy and native poses
    decoy_contact_map = get_contact_map_between_range1_range2( decoy_r1, decoy_r2, decoy, cutoff = cutoff, verbose = verbose )
    native_contact_map = get_contact_map_between_range1_range2( native_r1, native_r2, native, cutoff = cutoff, verbose = verbose )

    # get the number of contacts total found in the native
    num_native_contacts = sum( [ len( contacts ) for contacts in native_contact_map.values() ] )

    # determine the difference in contacts between the decoy and native
    num_decoy_recovered_contacts = 0
    for decoy_resnum in decoy_contact_map.keys():
        # get the corresponding native residue number
        corresponding_native_resnum = decoy_to_native_res_map[ decoy_resnum ]

        # for each contact this residue makes within the decoy
        for decoy_contact_resnum in decoy_contact_map[ decoy_resnum ]:
            # get the corresponding native residue number for the contact
            corresponding_native_contact_resnum = decoy_to_native_res_map[ decoy_contact_resnum ]

            # check to see if this is a contact made in the native pose
            try:
                native_contacts_at_this_decoy_resnum = native_contact_map[ corresponding_native_resnum ]

                # if the contact made in the decoy is the same as the contact made in the native, count it
                if corresponding_native_contact_resnum in native_contacts_at_this_decoy_resnum:
                    num_decoy_recovered_contacts += 1

            # if this residue doesn't make contacts in the native pose, skip it
            except KeyError:
                pass
            
    # calculate Fnat
    Fnat = round( ( float( num_decoy_recovered_contacts ) / float( num_native_contacts ) ) * 100, 2 )

    if return_more_data:
        num_decoy_contacts = sum( [ len( contacts ) for contacts in decoy_contact_map.values() ] )
        return num_decoy_contacts, num_native_contacts, Fnat
    else:
        return Fnat



def count_interface_residue_contacts( JUMP_NUM, pose, cutoff = CUTOFF_DISTANCE, return_more_data = False, verbose = False ):
    """
    Counts the number of residue "contacts" in a <pose> between interface residues given a <cutoff>
    :param JUMP_NUM: int( Jump number that defines the interface )
    :param pose: Pose
    :param cutoff: int( or float( cutoff distance in Angstroms). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param return_more_data: bool( do you want to have this function return fraction of protein-protein, protein-carbohydrate, and carbohydrate-carbohydrate contacts as well? ) Default is False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( number of residue contacts at interface ), list( str( res1resname_res1chain_res1pdbnum+res2resname_res2chain_res2pdbnum ) )
    :return: if return_more_data == True, float( interface res contacts ), list( unique contacts made ), float( fraction pro-pro ), float( fraction pro-carb ), float( fraction carb-carb )
    """
    # check to see that the jump number is valid
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print "You gave me an invalid jump number, try again"
        sys.exit()

    # get a list of all residue numbers from side 1 to 2
    range1 = []
    range2 = []

    # find all of the residue numbers that correspond to side 1 and side 2
    for ii in range( 1, pose.total_residue() + 1 ):
        if ii < pose.fold_tree().downstream_jump_residue( JUMP_NUM ):
            range1.append( ii )
        else:
            range2.append( ii )

    # check to see that neither of the lists are empty
    if len( range1 ) == 0 or len( range2 ) == 0:
        print
        print "It appears that the jump number", JUMP_NUM, "does not actually define an interface"
        sys.exit()

    # instantiate counter and the list to hold the unique contact names
    contacts = 0
    contact_list = []

    # instantiate contact type data lists
    pro_pro_contacts = 0
    pro_carb_contacts = 0
    carb_carb_contacts = 0
    
    # loop through each residue of each residue of side 1 specified and see if its within <cutoff> distance of any residue from side 2
    for res_num_1 in range1:
        # get the chain and seqpos of this residue
        res1 = pose.residue( res_num_1 )
        res1_name = res1.name3()
        res1_chain = pose.pdb_info().chain( res_num_1 )
        res1_pdb_num = pose.pdb_info().pose2pdb( res_num_1 ).split( ' ' )[0]
        
        # get the center of res_num_1
        res1_center = res1.nbr_atom_xyz()
        
        # loop over each residue on side 2
        for res_num_2 in range2:
            # get the chain and seqpos of this residue
            res2 = pose.residue( res_num_2 )
            res2_name = res2.name3()
            res2_chain = pose.pdb_info().chain( res_num_2 )
            res2_pdb_num = pose.pdb_info().pose2pdb( res_num_2 ).split( ' ' )[0]
            
            # get the center of res_num_2
            res2_center = res2.nbr_atom_xyz()

            # if the two residues are within the appropriate distance given by <cutoff>
            dist = res1_center.distance( res2_center )
            if dist <= cutoff:
                # makes unique names for each contact for ease of analysis
                uniq_name_res1 = res1_name + '_' + res1_chain + '_' + str( res1_pdb_num )
                uniq_name_res2 = res2_name + '_' + res2_chain + '_' + str( res2_pdb_num )

                # unique name is the combination of the information from both res1 and res2
                unique_name = uniq_name_res1 + '+' + uniq_name_res2
                
                # if this specific contact has not already been found, add it to the list and up the contacts counter
                if unique_name not in contact_list:
                    contact_list.append( unique_name )
                    contacts += 1

                ## record the type of contact made
                # protein to protein
                if res1.is_protein() and res2.is_protein():
                    pro_pro_contacts += 1
                # carbhydrate to carbohydrate
                elif res1.is_carbohydrate() and res2.is_carbohydrate():
                    carb_carb_contacts += 1
                # protein to carbohydrate
                else:
                    pro_carb_contacts += 1

    # calculate fraction of contact types
    pro_pro_fraction = round( float( pro_pro_contacts ) / float( contacts ), 2 )
    pro_carb_fraction = round( float( pro_carb_contacts ) / float( contacts ), 2 )
    carb_carb_fraction = round( float( carb_carb_contacts ) / float( contacts ), 2 )

    if return_more_data:
        return contacts, contact_list, pro_pro_fraction, pro_carb_fraction, carb_carb_fraction
    else:
        return contacts, contact_list



def count_residue_contacts_between_range1_range2( range1, range2, pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Counts the number of residue "contacts" in a <pose> between <range1> and <range2>
    :param range1: list( residue pose numbers in range 1 )
    :param range2: list( residue pose numbers in range 2 )
    :param pose: Pose
    :param cutoff: int( or float( cutoff distance in Angstroms). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( number of residue contacts at interface ), list( str( res1resname_res1chain_res1pdbnum+res2resname_res2chain_res2pdbnum ) )
    """
    # instantiate counter and the list to hold the unique contact names
    contacts = 0
    contact_list = []
    
    # loop through each residue of each residue of side 1 specified and see if its within <cutoff> distance of any residue from side 2
    for res_num_1 in range1:
        # get the chain and seqpos of this residue
        res1 = pose.residue( res_num_1 )
        res1_name = res1.name3()
        res1_chain = pose.pdb_info().chain( res_num_1 )
        res1_pdb_num = pose.pdb_info().pose2pdb( res_num_1 ).split( ' ' )[0]
        
        # get the center of res_num_1
        res1_center = res1.nbr_atom_xyz()
        
        # loop over each residue on side 2
        for res_num_2 in range2:
            # get the chain and seqpos of this residue
            res2 = pose.residue( res_num_2 )
            res2_name = res2.name3()
            res2_chain = pose.pdb_info().chain( res_num_2 )
            res2_pdb_num = pose.pdb_info().pose2pdb( res_num_2 ).split( ' ' )[0]
            
            # get the center of res_num_2
            res2_center = res2.nbr_atom_xyz()

            # if the two residues are within the appropriate distance given by <cutoff>
            dist = res1_center.distance( res2_center )
            if dist <= cutoff:
                # makes unique names for each contact for ease of analysis
                uniq_name_res1 = res1_name + '_' + res1_chain + '_' + str( res1_pdb_num )
                uniq_name_res2 = res2_name + '_' + res2_chain + '_' + str( res2_pdb_num )

                # unique name is the combination of the information from both res1 and res2
                unique_name = uniq_name_res1 + '+' + uniq_name_res2
                
                # if this specific contact has not already been found, add it to the list and up the contacts counter
                if unique_name not in contact_list:
                    contact_list.append( unique_name )
                    contacts += 1

    return contacts, contact_list



def count_atomic_contacts_between_range1_range2( range1, range2, pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Counts the number of non-hydrogen atomic contacts in a <pose> between residues of <range1> and <range2> given the <cutoff> distance
    :param range1: list of ints( Pose residue numbers of side 1 )
    :param range2: list of ints( Pose residue numbers of side 2 )
    :param pose: Pose
    :param cutoff: int( or float( cutoff distance in Angstroms). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( number of atomic contacts at interface ), list( str( atm1 res name + atm name + atm1 chain + '_' + atm2 res name + atm2 atm name + atm2 chain
    """
    if verbose:
        print "Finding atomic contacts (except hydrogens) between", len( range1 ), "residues on side 1 and", len( range2 ), "residues on side 2"

    # instantiate counter and the list to hold the unique contact names
    contacts = 0
    contact_list = []
    
    # loop through each atom of each residue of side 1 specified and see if its within <cutoff> distance of any atom from side 2
    for res_num_1 in range1:
        # get the number of atoms in residue 1
        n_atms_1 = pose.residue( res_num_1 ).natoms()

        # loop over each atom
        for atm_1 in range( 1, n_atms_1 + 1 ):
            # if it's not a hydrogen atom, get its...
            if not pose.residue( res_num_1 ).atom_is_hydrogen( atm_1 ):
                # xyz coordinates
                atm1_xyz = list( pose.residue( res_num_1 ).atom( atm_1 ).xyz() )
                # residue name
                atm1_res_name = pose.residue( res_num_1 ).name3()
                # atom name
                atm1_atm_name = pose.residue( res_num_1 ).atom_name( atm_1 ).replace( ' ', '' )
                # and which chain it's on
                atm1_chain = pose.pdb_info().chain( res_num_1 )

                # do the same for residue 2
                for res_num_2 in range2:
                    # get the number of atoms in residue 2
                    n_atms_2 = pose.residue( res_num_2 ).natoms()

                    for atm_2 in range( 1, n_atms_2 + 1 ):
                        # if it's not a hydrogen atom, get its...
                        if not pose.residue( res_num_2 ).atom_is_hydrogen( atm_2 ):
                            # xyz coordinates
                            atm2_xyz = list( pose.residue( res_num_2 ).atom( atm_2 ).xyz() )
                            # residue name
                            atm2_res_name = pose.residue( res_num_2 ).name3()
                            # atom name
                            atm2_atm_name = pose.residue( res_num_2 ).atom_name( atm_2 ).replace( ' ', '' )
                            # and which chain it's on
                            atm2_chain = pose.pdb_info().chain( res_num_2 )

                            # if the two atoms are within the appropriate distance given by <cutoff>
                            atm1_atm2_dist = calc_distance( atm1_xyz, atm2_xyz )
                            if atm1_atm2_dist <= cutoff:
                                # makes unique names for each contact for ease of analysis
                                uniq_name_atm1 = atm1_res_name + '_' + atm1_atm_name + '_' + atm1_chain + '_' + str( res_num_1 )
                                uniq_name_atm2 = atm2_res_name + '_' + atm2_atm_name + '_' + atm2_chain + '_' + str( res_num_2 )

                                # unique name is the combination of the atom information from both atom1 and atom2
                                unique_name = uniq_name_atm1 + '.' + uniq_name_atm2

                                # if this specific contact has not already been found, add it to the list and up the contacts counter
                                if unique_name not in contact_list:
                                    contact_list.append( unique_name )
                                    contacts += 1

    return contacts, contact_list



def count_interface_atomic_contacts( JUMP_NUM, pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Counts the atom-to-atom (discluding hydrogens) contacts given a cutoff of <cutoff> Angstroms at the interface given the jump number <JUMP_NUM>
    Returns the number of contacts as a float and a list of the unique contacts made between the atoms of side 1 to side 2
    :param JUMP_NUM: int( Jump number that defines the interface )
    :param pose: Pose
    :param cutoff: int( or float( cutoff distance in Angstroms). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( number of atomic contacts at interface ), list( str( atm1 res name + atm name + atm1 chain + '_' + atm2 res name + atm2 atm name + atm2 chain
    """
    # check to see that the jump number is valid
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print "You gave me an invalid jump number, try again"
        sys.exit()

    if verbose:
        print "Counting interface contacts across jump number", JUMP_NUM, "discluding hydrogen atoms..."

    # get a list of all residue numbers from side 1 to 2
    side_1_list = []
    side_2_list = []

    # find all of the residue numbers that correspond to side 1 and side 2
    for ii in range( 1, pose.total_residue() + 1 ):
        if ii < pose.fold_tree().downstream_jump_residue( JUMP_NUM ):
            side_1_list.append( ii )
        else:
            side_2_list.append( ii )

    # check to see that neither of the lists are empty
    if len( side_1_list ) == 0 or len( side_2_list ) == 0:
        print
        print "It appears that the jump number", JUMP_NUM, "does not actually define an interface"
        sys.exit()

    # find all residues from side 1 and side 2 that are within 3 * CUTOFF_DISTANCE angstroms of each other
    # speeds up the counting process a bit
    side_1_list_close = []
    side_2_list_close = []

    # triple the given <cutoff> distance
    TRIPLE_CUTOFF_DISTANCE = cutoff * 3

    if verbose:
        print "Finding the residues within 3 times the given cutoff distance to speed up counting a bit"

    # loop over each residue from side 1
    for res_num_1 in side_1_list:
        # find the center of the residue in side 1
        center = pose.residue( res_num_1 ).nbr_atom_xyz()

        # now loop over the residues in side 2 for comparison
        for res_num_2 in side_2_list:
            # compare the distance to the other residue in side 2  -  using triple the distance at first to theoretically speed up the counting process
            if center.distance( pose.residue( res_num_2 ).nbr_atom_xyz() ) <= TRIPLE_CUTOFF_DISTANCE:
                # if the residues are within triple the <cutoff>, and if they haven't already been added, add them to the close list
                if res_num_1 not in side_1_list_close:
                    side_1_list_close.append( res_num_1 )

                if res_num_2 not in side_2_list_close:
                    side_2_list_close.append( res_num_2 )

    # loops over side 1 list and see which atoms are within 5 Ang of side 2 for native pose
    num_interface_contacts, interface_contact_list = count_atomic_contacts_between_range1_range2( side_1_list_close, side_2_list_close, pose, cutoff )

    return float( num_interface_contacts ), interface_contact_list



def calc_interface_sasa( pose, JUMP_NUM ):
    """
    Use rosetta.calc_total_sasa to compute the SASA of the total pose - SASA of the split-apart pose at <JUMP_NUM>
    :param pose: Pose
    :param JUMP_NUM: int( valid Jump number defining interface )
    :return: float( interface_SASA value )
    """
    # imports
    from rosetta import calc_total_sasa
    
    # make sure a valid JUMP_NUM was passed in
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print JUMP_NUM, "is an invalid Jump number. Exiting"
        sys.exit()
    
    # make and split the temporary pose
    temp_pose = Pose()
    temp_pose.assign( pose )
    jump = temp_pose.jump( JUMP_NUM )

    #TODO-get current xyz location and multiply by 500 or something instead
    vec = xyzVector_Real( 1000, 1000, 1000 )
    jump.set_translation( vec )
    temp_pose.set_jump( JUMP_NUM, jump )
    
    # calculate the SASA values
    total_sasa = calc_total_sasa( pose, PROBE_RADIUS )
    split_sasa = calc_total_sasa( temp_pose, PROBE_RADIUS )
    delta_interface_sasa = total_sasa - split_sasa
    
    return delta_interface_sasa



def analyze_interface( pose, JUMP_NUM, pack_separated = True ):
    """
    Use rosetta.protocols.analysis.Interface Analyzer to compute various interface metrics
    :param pose: Pose
    :param JUMP_NUM: int( valid Jump number defining interface )
    :param pack_separated: bool( Do you want to pack the protein after you split them apart? ). Default = True
    :return: float( interface_SASA value )
    """
    # make sure a valid JUMP_NUM was passed in
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print JUMP_NUM, "is an invalid Jump number. Exiting"
        sys.exit()

    # instantiate an InterfaceAnalyzer mover
    IAmover = IAM()

    # set the interface jump to the jump number passed in
    IAmover.set_interface_jump( JUMP_NUM )

    # TODO-see what's worth calculating and then actually return the value
    # set what you want to calculate
    ##IAmover.set_compute_interface_delta_hbond_unsat( True )
    ##IAmover.set_compute_interface_sc( True )
    IAmover.set_compute_separated_sasa( True )

    # set the relevant options
    IAmover.set_input_pose( pose )
    IAmover.set_pack_separated( pack_separated )

    print "Analyzing interface..."
    IAmover.reset_status()
    IAmover.apply( pose )

    # retrieve data
    # TODO-retrieve new, relevant data
    interface_dSASA = IAmover.get_interface_delta_sasa()
    ##unsat_hbond =  IAmover.get_interface_delta_hbond_unsat()
    ##interface_dG = IAmover.get_interface_dG()  # I already calculate this with the fxn I wrote
    ##num_interface_residues = IAmover.get_num_interface_residues()


    return interface_dSASA



def get_interface_score( JUMP_NUM, sf, pose ):
    """
    Given a jump number that defines the interface, calculates Rosetta's ddG interface
    Splits apart the two domains defined by the <JUMP_NUM>, scores it, then subtracts that from the total score of the <pose>  -  result is the interface score
    :param JUMP_NUM: int( valid Jump number of the interface )
    :param sf: ScoreFunction
    :param pose: Pose
    :return: float( ddG interface score )
    """
    # get start score
    start_score = sf( pose )

    # make and split the temporary pose
    temp_pose = Pose()
    temp_pose.assign( pose )
    jump = temp_pose.jump( JUMP_NUM )

    #TODO-get current xyz location and multiply by 500 or something instead
    vec = xyzVector_Real( 1000, 1000, 1000 )
    jump.set_translation( vec )
    temp_pose.set_jump( JUMP_NUM, jump )

    # get and return interface score
    split_apart_score = sf( temp_pose )
    interface_score = start_score - split_apart_score

    return interface_score



def get_phi_psi_omega_of_res( pose, seqpos ):
    """
    Returns the phi, psi, and omega values of the passed residue <seqpos> in <pose>
    :param pose: Pose
    :param seqpos: int( the sequence position of the residue of interest )
    :return: float( phi ), float( psi ), float( omega )
    """
    try:
        return pose.phi( seqpos ), pose.psi( seqpos ), pose.omega( seqpos )
    except:
        print "%s does not seem to be a residue in the passed Pose. Check your input" %seqpos
        return None
        
    

def determine_amino_acid_composition( pose ):
    """
    Prints the amino acid composition of the <pose> using Pandas DataFrame
    Ignores sugars ( single letter code 'Z' ) in composition analysis
    :param pose: Pose
    :return: a DataFrame or dictionary of the data, depending on if Pandas can be imported
    """
    # initialize a dictionary for the data
    # won't be used if Pandas import worked
    data_dict = {}

    # build lists for pandas DataFrame
    res_count = []
    percentage = []

    # get number of amino acid residues
    res_total = ( pose.total_residue() - pose.sequence().count( 'Z' ) )

    for res_type in AA_list:
        count = float( pose.sequence().count( res_type ) )
        res_count.append( count )
        percentage.append( round( ( count / res_total * 100 ), 2 ) )

    # if global Pandas import was successful, create a DataFrame
    if pandas_on:
        df_AA_composition = pd.DataFrame()
        df_AA_composition["Res_Type"] = AA_list
        df_AA_composition["Count"] = res_count
        df_AA_composition["Percentage"] = percentage
        
    # else create a dictionary
    else:
        data_dict["Res_Type"] = AA_list
        data_dict["Count"] = res_count
        data_dict["Percentage"] = percentage
        
    # print out relevant data
    print
    print "This protein has..."
    if pandas_on:
        print df_AA_composition.sort( "Percentage", ascending = False )
        return df_AA_composition
    else:
        print data_dict
        return data_dict



def compare_all_rotamers( pose1, pose2, E_diff = 5 ):
    """
    Uses is_similar_rotamer on each residue, only if they are the same amino acids
    Will exit if the two Poses don't have the same number of amino acids - can't guarantee the validity of the comparison otherwise
    Skips over sugars and glycines
    :param pose1: Pose 1
    :param pose2: Pose 2
    :param E_diff: int( or float( the energy difference required to be considered different ). Default = +/- 5
    :return: a DataFrame or dictionary of the data, depending on if Pandas can be imported (does not print)
    """
    # ensure that the two poses have the same number of residues
    if not pose1.total_residue() == pose2.total_residue():
        print
        print "Check your poses, they don't have the same number of residues. Exiting"
        sys.exit()

    # get the absolute value of the <E_diff>, and also get its negation
    E_diff = abs( E_diff )
    neg_E_diff = -1 * E_diff

    # make the lists to hold data
    resi_num_pose = []
    resi_chain = []
    resi_name = []

    # compare each residue rotamer
    total_residues = pose1.total_residue()
    for ii in range( 1, total_residues + 1 ):
        res1_AA = pose1.residue( ii ).name1()
        res2_AA = pose2.residue( ii ).name1()

        # only compare if the two residues are the same amino acid, and are not sugars (I'm assuming sugars will change enough)
        if res1_AA == res2_AA:
            res1 = pose1.residue( ii )
            res2 = pose2.residue( ii )

            # if they aren't sugars
            if not res1.is_carbohydrate():
                if not res1.name1() != 'G':
                    if not res2.is_carbohydrate():
                        if not res2.name1() != 'G':
                            for chi in range( 1, len( res1.chi() ) + 1 ):
                                diff = res2.chi()[ chi ] - res1.chi()[ chi ]
                                if diff >= E_diff or diff <= neg_E_diff:
                                    resi_num_pose.append( ii )
                                    # get the chain of the residue by going from pose --> PDB numbering and getting it from there
                                    resi_chain.append( pose1.pdb_info().pose2pdb( ii )[ -2:len( pose1.pdb_info().pose2pdb( ii ) ) ] )
                                    resi_name.append( res1.name1() )

    # make and return Pandas DataFrame if Pandas was imported
    if pandas_on:
        df = pd.DataFrame()
        df["Pose Num"] = resi_num_pose
        df["Chain"] = resi_chain
        df["AA"] = resi_name
        
        return df
    
    # if Pandas wasn't imported, make and return a dictionary of the data
    else:
        # initialize a dictionary for the data
        data_dict = {}
        
        # append data
        data_dict["Pose Num"] = resi_num_pose
        data_dict["Chain"] = resi_chain
        data_dict["AA"] = resi_name
        
        return data_dict



def check_E_per_residue( sf, pose, energy_cutoff = 1.5, verbose = False ):
    """
    Returns a list of the residues with a total energy value greater than <energy_cutoff>
    :param sf: ScoreFunction
    :param pose: Pose
    :param energy_cutoff: int( or float( residue energy cutoff where residues of interest have a greater energy value ) ). Default = 1.5
    :param verbose: bool( if you want the function to print out the high-energy residues, set to True ). Default = False
    :return: list( Residue objects of residues with total energy above <energy_cutoff> )
    """
    # get the absolute value of the passed energy_cutoff
    abs_cutoff = abs( energy_cutoff )
    
    # score the pose to get access to its energy data
    sf( pose )
    
    # iterate through each residue in pose looking at its energy
    high_E_residues = []
    for residue in pose:
        energy = pose.energies().residue_total_energy( residue.seqpos() )
        if energy > abs_cutoff:
            if verbose:
                print residue.name(), '\t', pose.pdb_info().pose2pdb( residue.seqpos() ), '\t', energy
            high_E_residues.append( residue )
            
    return high_E_residues





############################
#### INITIALIZE ROSETTA ####
############################

if __name__ == '__main__':
    # initialize rosetta with sugar flags
    from rosetta import init
    
    print "Initializing Rosetta with sugar flags"
    #init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records" )
    init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -write_pdb_link_records" )

############################
#### INITIALIZE ROSETTA ####
############################
