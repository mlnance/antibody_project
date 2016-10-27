#!/usr/bin/python
__author__="morganlnance"

'''
Various functions to be used in 3ay4 glycan modeling
Most compiled using the worker functions from antibody_functions
'''


############################
#### CONSTANT VARIABLES ####
############################

## Pose numbering information ONLY relevant to native PDB 3ay4
native_Fc_chain_A_nums = range( 1, 215 + 1 )
native_Fc_glycan_A_nums = range( 216, 223 + 1 )
native_Fc_glycan_A_nums_except_core_GlcNAc = range( 217, 223 + 1 )
native_Fc_chain_B_nums = range( 224, 439 + 1 )
native_Fc_glycan_B_nums = range( 440, 447 + 1 )
native_Fc_glycan_B_nums_except_core_GlcNAc = range( 441, 447 + 1 )
native_FcR_protein_nums = range( 448, 607 + 1 )
native_FcR_main_glycan_nums = range( 608, 615 + 1 )
native_FcR_three_mer_nums = range( 616, 618 + 1 )
native_FcR_glycan_nums = range( 608, 618 + 1 )
native_branch_points = [ 69, 218, 292, 442, 478, 595, 608, 610 ]
native_Fc_glycan_branch_point_nums = [ 218, 442 ]
native_Fc_glycan_branch_point_nums_with_ASN = [ 69, 218, 292, 442 ]
native_Fc_protein_nums = []
native_Fc_protein_nums.extend( native_Fc_chain_A_nums )
native_Fc_protein_nums.extend( native_Fc_chain_B_nums )
native_Fc_glycan_nums = []
native_Fc_glycan_nums.extend( native_Fc_glycan_A_nums )
native_Fc_glycan_nums.extend( native_Fc_glycan_B_nums )
native_Fc_glycan_nums_except_core_GlcNAc = []
native_Fc_glycan_nums_except_core_GlcNAc.extend( native_Fc_glycan_A_nums_except_core_GlcNAc )
native_Fc_glycan_nums_except_core_GlcNAc.extend( native_Fc_glycan_B_nums_except_core_GlcNAc )
native_order_nums = range( 1, 618 + 1 )
native_Fc_protein_chains = [ 'A', 'B' ]
native_FcR_protein_chains = [ 'C' ]
native_Fc_glycan_chains = [ 'D', 'E', 'F', 'G' ]
native_Fc_glycan_A_chains = [ 'D', 'E' ]
native_Fc_glycan_B_chains = [ 'F', 'G' ]
native_FcR_glycan_chains = [ 'H', 'I', 'J', 'K' ]



###################
#### FUNCTIONS ####
###################

def initialize_rosetta( constant_seed = False, debug = False ):
    """
    Initialize Rosetta and mute basic, core, and protocols.
    If constant_seed == True, use default constant seed 1111111
    If debug == True, use default constant seed and do not mute Rosetta
    """
    # imports
    from rosetta import init


    print "Initializing Rosetta with sugar flags"

    # makes Rosetta quiet and sugar I/O ready
    #init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records" )
    if constant_seed:
        init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -write_pdb_link_records -constant_seed" )
    elif debug:
        init( extra_options="-include_sugars -override_rsd_type_limit -write_pdb_link_records -constant_seed -out:level 400" )
    else:
        init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -write_pdb_link_records" )



def load_pose( pose_filename ):
    """
    Load pose from a filename
    :param pose_filename: str( /path/to/pose/filename )
    :return: a Rosetta Pose
    """
    # imports
    from rosetta import Pose, pose_from_file, FoldTree

    
    # create Pose object from filename
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
    # imports
    import sys
    from rosetta import get_fa_scorefxn, score_type_from_name
    from rosetta.core.scoring import fa_atr


    # argument check - check the passed argument is a dict
    if not isinstance( weights_dict, dict ):
        print "You didn't give me a dictionary for your input. I need a dict of ScoreType (or name) : weight. Exiting."
        sys.exit()
        
    # get a standard fa_scorefxn to start with
    sf = get_fa_scorefxn()
    
    # for each entry of the dictionary, change the weight
    for scoretype_name, scoretype_weight in weights_dict.items():
        # if the key is a string
        if isinstance( scoretype_name, str ):
            try:
                scoretype = score_type_from_name( scoretype_name )
            except:
                print "\nThe string name: '%s' does not appear to be a valid ScoreType. Exiting" %scoretype_name
                sys.exit()
            
            # set the weight
            sf.set_weight( scoretype, scoretype_weight )

        # if the argument is a ScoreType object
        elif isinstance( scoretype_name, type( fa_atr ) ):
            # adjust the weight in the scorefxn using the corresponding weight given
            sf.set_weight( scoretype_name, scoretype_weight )

        # else, I don't know what they gave me as a scoretype
        else: 
            print "I'm not sure what '%s' is from your ScoreType key in your <weights_dict> argument. Exiting" %scoretype_name
            sys.exit()

    # if verbose, print out the weights of the new ScoreFunction
    if verbose:
        print "\nNew score weights sf:\n%s\n" %( "\n".join( [ "%s: %s" %( str( name ), sf.get_weight( name ) ) for name in sf.get_nonzero_weighted_scoretypes() ] ) )

    # return the newly weighted fa_scorefxn
    return sf



def native_3ay4_Fc_glycan_LCM_reset( mm, input_pose, use_population_ideal_LCM_reset = False ):
    '''
    Reset the 3ay4 Fc glycan as determined by carbohydrate residues with BB set to True in MoveMap <mm>
    Default action is an LCM_reset with stdev +/- 1
    Uses data from the LCM default.table
    <use_population_ideal_LCM_reset> makes the reset only use the ideal values (with appropriate population weights). Essentially, stdev = 0
    :param mm: MoveMap ( residues with BB set to True and that are carbohydrates are available to be reset )
    :param input_pose: Pose
    :param use_population_ideal_LCM_reset: bool( use population ideals only? ) Default = False
    :return: Pose
    '''
    # TODO add this option back?
    #:param use_ideal_LCM_reset: bool( use only the ideal value from the highest population? ) Default = False

    # imports
    import sys
    from rosetta import MoveMap
    from rosetta.protocols.carbohydrates import LinkageConformerMover


    # check input arguments
    #if use_population_ideal_LCM_reset and use_ideal_LCM_reset:
    #    print "\nYou want me to use both the population ideal and the main ideal. This does not make sense to me. Ensure you are only setting one of these arguments to True. Exiting\n"
    #    sys.exit()

    # copy in the input_pose
    pose = input_pose.clone()

    # residues that are allowed to move are based on residues with BB set to True in the MoveMap and the residue is a carbohydrate
    # I'm copying this idea from GlycanRelaxMover, but overall it makes sense anyway. Only things with BB freedom should be sampled
    carbohydrate_res_nums = [ res_num for res_num in range( 1, pose.n_residue() + 1 ) if mm.get_bb( res_num ) and pose.residue( res_num ).is_carbohydrate() ]

    # for each residue in <carbohydrate_res_nums>
    for res_num in carbohydrate_res_nums:
        # make a BB MoveMap for this single residue
        lcm_mm = MoveMap()
        lcm_mm.set_bb( res_num, True )

        # create a LinkageConformerMover with the residue MoveMap
        lcm = LinkageConformerMover()
        lcm.set_movemap( lcm_mm )

        # if the user only wants to use ideals, but within the different population clusters
        if use_population_ideal_LCM_reset:
            # set_idealize_torsions uses ideal values instead of sampling from stdev
            lcm.set_idealize_torsions( True )
            # use_conformer_population_stats is True by default, but setting for clarity
            lcm.set_use_conformer_population_stats( True )
            # apply the LCM
            lcm.apply( pose )

        # else, standard LCM reset using (by default) population data and a stdev of 1
        else:
            # setting these options for clarity
            lcm.set_use_conformer_population_stats( True )
            lcm.set_x_standard_deviations( 1 )
            # apply the LCM
            lcm.apply( pose )

    return pose



def native_3ay4_Fc_glycan_random_reset( mm, input_pose ):
    '''
    Reset the 3ay4 Fc glycan randomly. Each torsion is chosen from -180 to 180
    :param mm: MoveMap ( residues with BB set to True and that are carbohydrates are available to be reset )
    :param input_pose: Pose
    :return: Pose
    '''
    # imports
    import sys
    from random import uniform
    from rosetta import MoveMap


    # copy in the input_pose
    pose = input_pose.clone()

    # residues that are allowed to move are based on residues with BB set to True in the MoveMap and the residue is a carbohydrate
    # I'm copying this idea from GlycanRelaxMover, but overall it makes sense anyway. Only things with BB freedom should be sampled
    carbohydrate_res_nums = [ res_num for res_num in range( 1, pose.n_residue() + 1 ) if mm.get_bb( res_num ) and pose.residue( res_num ).is_carbohydrate() ]

    # for each residue in <carbohydrate_res_nums>
    for res_num in carbohydrate_res_nums:
        # reset the phi, psi, and omega to some random number between -180 and 180
        pose.set_phi( res_num, uniform( -180, 180 ) )
        pose.set_psi( res_num, uniform( -180, 180 ) )
        pose.set_omega( res_num, uniform( -180, 180 ) )

    return pose



def add_constraints_to_pose( constraint_file, input_pose ):
    '''
    Add the constraints found in the <constraint_file> to the <input_pose>
    Don't forget to turn on the appropriate ScoreFunction weights yourself! (like atom_pair_constraint or dihedral_constraint)
    :param constraint_file: str( /path/to/cst file )
    :param input_pose: Pose
    :return: Pose
    '''
    # imports
    from rosetta.protocols.simple_moves import ConstraintSetMover


    # copy the input_pose
    pose = input_pose.clone()

    # set the constraint from the file
    constraint_setter = ConstraintSetMover()
    constraint_setter.constraint_file( constraint_file )
    constraint_setter.apply( pose )

    return pose



def SugarSmallMover( mm, nmoves, angle_max, input_pose, move_all_torsions = True, use_sugar_bb_angle_max = False ):
    """
    Randomly resets the phi, psi, and omega values of <nmoves> residues found in the <mm> MoveMap that are carbohydrates and have their BB freedom turned on
    Math works out to where the max motion to either direction or starting position is angle_max/2 ( old_value +/- angle_max/2 )
    Emulates the SmallMover but with the additional omega mover
    Can either move all torsions found on the chosen residue (Default action), or can set move_all_torsions to False and will instead pick one of the available torsions of that residue to perturb
    use_sugar_bb_angle_max means to use an angle max for phi and psi (and omega that Jason told me) that comes from the sugar_bb graph (just eyeballed)
    :param mm: MoveMap ( residues with BB set to True and that are carbohydrates are available to be moved )
    :param nmoves: int( how many moves should be allowed in one call to the SugarSmallMover? )
    :param angle_max: int( or float( the max angle around which the phi/psi/omega could move ) )
    :param input_pose: Pose
    :param move_all_torsions: bool( do you want to move all BackBone torsions of the residues at the same time? If not, a single torsion will be randomly chosen and perturbed ) Default = True
    :param use_sugar_bb_angle_max: bool( use sugar_bb phi/psi/omega data to get angle_max ) Default = False
    :return: Pose
    """
    # imports
    from random import choice
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg
    from rosetta.core.id import MainchainTorsionType, phi_dihedral, psi_dihedral, omega_dihedral
    from rosetta.core.pose.carbohydrates import get_glycosidic_torsion, set_glycosidic_torsion, \
        get_reference_atoms


    # copy the input pose
    pose = input_pose.clone()

    # residues that are allowed to move are based on residues with BB set to True in the MoveMap and the residue is a carbohydrate
    # I'm copying this idea from GlycanRelaxMover, but overall it makes sense anyway. Only things with BB freedom should be sampled
    moveable_res_nums = [ res_num for res_num in range( 1, pose.n_residue() + 1 ) if mm.get_bb( res_num ) and pose.residue( res_num ).is_carbohydrate() ]

    # from rosetta.protocols.simple_moves:BackboneMover.cc file for SmallMover
    big_angle = angle_max
    small_angle = big_angle / 2.0

    # for as many moves as specified
    for ii in range( nmoves ):
        # pick a residue to sample
        res_num = choice( moveable_res_nums )

        # check which glycosidic torsions it has (except omega3_dihedral, that doesn't work at the moment)
        moveable_torsions = []
        for torsion_name in MainchainTorsionType.names:
            # have to skip omega3, hence doing by name
            if torsion_name != "omega3_dihedral":
                # if this torsion has reference atoms, it exists
                # bool( len([]) ) = False, so it will return False if there are no reference atoms (ie. it doesn't exist)
                if bool( len( get_reference_atoms( MainchainTorsionType.names[ torsion_name ], pose, res_num ) ) ):
                    # this torsion type has reference atoms, thus it exists. Keep it
                    moveable_torsions.append( MainchainTorsionType.names[ torsion_name ] )
        # keep one torsion if the user doesn't want all torsions sampled
        if move_all_torsions is False:
            moveable_torsions = [ choice( moveable_torsions ) ]

        # get the current torsions for the moveable_torsions and perturb them according to angle_max
        for moveable_torsion in moveable_torsions:
            # get the current torsion value
            old_torsion_value = get_glycosidic_torsion( moveable_torsion, pose, res_num )

            # if use_sugar_bb_angle_max is set to True, change the angle_max, small_angle, and big_angle
            # depending on if the torsion being sampled is a phi, psi, or omega
            if use_sugar_bb_angle_max:
                # phi and omega
                if moveable_torsion == phi_dihedral or moveable_torsion == omega_dihedral:
                    big_angle = 30
                    small_angle = big_angle / 2.0
                # psi
                elif moveable_torsion == psi_dihedral:
                    big_angle = 100
                    small_angle = big_angle / 2.0

            # perturb this torsion randomly
            # this specific format is pulled from rosetta.protocols.simple_moves:ShearMover::make_move
            new_torsion_value = periodic_range( old_torsion_value - small_angle + rg().uniform() * big_angle, 360.0 )

            # set the new torsion
            set_glycosidic_torsion( moveable_torsion, pose, res_num, new_torsion_value )

    return pose



def get_ramp_score_weight( current_weight, target_weight, current_step, total_steps ):
    """
    Given the <current_weight> and the <target_weight>, use the <current_step> and the <total_steps> to determine how much the weight should be increased or decreased for this particular round
    Current and Total steps -- Say you're doing 100 (1-100) rounds, if you're on round 57, current_step = 57, total_steps = 100
    :param current_weight: int( or float( value of your current ScoreType weight ) )
    :param target_weight: int( or float( value of your target ScoreType weight ) )
    :param current_step: int( current step of how many rounds you're doing )
    :param total_steps: int( number of rounds you're doing )
    :return: float( the new weight to set for your ScoreType for this particular round )
    """
    # imports
    import sys


    # need to be able to do this a minimum of 10 times, otherwise the math will break
    if total_steps < 10:
        print
        print "You need to run a loop at least more than 10 times - otherwise the math in here won't work"
        sys.exit()

    # adjust the total_steps so that it's actually 10% less than what it actually is
    # this is so that the final 10% of the simulation will run with the target score weights
    total_steps = int( total_steps * 0.9 )

    # make appropriate weight adjustments
    if current_weight > target_weight:  # we're decreasing the weight
        amount_left = target_weight - current_weight
        moves_left = float( total_steps - current_step )
        add_weight = amount_left / moves_left  # should be negative
        new_weight = current_weight + add_weight

    elif target_weight > current_weight:  # we're increasing the weight
        amount_left = target_weight - current_weight
        moves_left = float( total_steps - current_step )
        add_weight = amount_left / moves_left  # should be positive
        new_weight = current_weight + add_weight

    else:
        new_weight = current_weight

    return new_weight



def get_ramp_angle_max( current_angle_max, target_angle_max, current_step, total_steps ):
    """
    Given the <current_angle_max> and the <target_angle_max>, use the <current_step> and the <total_steps> to determine how much the angle_max should be increased or decreased for this particular round
    Current and Total steps -- Say you're doing 100 (1-100) rounds, if you're on round 57, current_step = 57, total_steps = 100
    :param current_angle_max: int( or float( value of your current angle_max ) )
    :param target_angle_max: int( or float( value of your target angle_max ) )
    :param current_step: int( current step of how many rounds you're doing )
    :param total_steps: int( number of rounds you're doing )
    :return: float( the new angle angle_max to set for your SugarSmall/ShearMover for this particular round )
    """
    # imports
    import sys


    # need to be able to do this a minimum of 10 times, otherwise the math will break
    if total_steps < 10:
        print
        print "You need to run a loop at least more than 10 times - otherwise the math in here won't work"
        sys.exit()

    # adjust the total_steps so that it's actually 10% less than what it actually is
    # this is so that the final 10% of the simulation will run with the target score weights
    total_steps = int( total_steps * 0.9 )

    # make appropriate angle_max adjustments
    if current_angle_max > target_angle_max:  # we're decreasing the angle_max
        amount_left = target_angle_max - current_angle_max
        moves_left = float( total_steps - current_step )
        add_angle_max = amount_left / moves_left  # should be negative
        new_angle_max = current_angle_max + add_angle_max

    elif target_angle_max > current_angle_max:  # we're increasing the angle_max
        amount_left = target_angle_max - current_angle_max
        moves_left = float( total_steps - current_step )
        add_angle_max = amount_left / moves_left  # should be positive
        new_angle_max = current_angle_max + add_angle_max

    else:
        new_angle_max = current_angle_max

    return new_angle_max



def get_res_nums_within_radius( res_num_in, input_pose, radius, include_res_num = False ):
    """
    Use the nbr_atom_xyz to find residue numbers within <radius> of <pose_num> in <pose>
    :param res_num_in: int( Pose residue number )
    :param input_pose: Pose
    :param radius: int or float( radius around <pose_num> to use to select resiudes )
    :param include_res_num: bool( do you want to include <res_num> in the return list? ) Default = False
    :return: list( Pose residue numbers within <radius> of <pose_num>
    """
    # clone the <input_pose>
    pose = input_pose.clone()

    # container for the centers of each residue in pose
    centers_of_res = []

    # fill up the centers container
    for res_num in range( 1, pose.n_residue() + 1 ):
        center = pose.residue( res_num ).nbr_atom_xyz()
        centers_of_res.append( center )

    # container for residues inside the <radius>
    res_nums_in_radius = []

    # nbr_xyz of the residue of interest
    res_num_xyz = pose.residue( res_num_in ).nbr_atom_xyz()

    for res_num in range( 1, pose.n_residue() + 1 ):
        # this will get the xyz of the residue of interest, but it will be removed from the final list if desired
        # (since it will be added as 0 will always be less than <radius>)
        # get the center of the residue
        center = pose.residue( res_num ).nbr_atom_xyz()

        # keep the residue number if the nbr_atom_xyz is less than <radius>
        if center.distance( res_num_xyz ) <= radius:
            res_nums_in_radius.append( res_num )

    # if the user didn't want the residue of interest in the return list, remove it
    if not include_res_num:
        res_nums_in_radius.remove( res_num_in )

    return res_nums_in_radius



def get_res_nums_within_radius_of_residue_list( residues, input_pose, radius, include_res_nums = False ):
    """
    Find all residue numbers around the list of <residues> given in <input_pose> within <radius> Angstroms.
    Set <include_residues> if you want to include the list of passed <residues> in the return list of residue numbers.
    :param residues: list( Pose residue numbers )
    :param input_pose: Pose
    :param radius: int() or float( radius in Angstroms )
    :param include_res_nums: bool( do you want to include the passed <residues> in the return list of resiude numbers? ) Default = False
    :return: list( residues around passed <residues> list within <radius> Angstroms
    """
    # argument check: ensure passed <residues> argument is a list
    if type( residues ) != list:
        print "\nArgument error. You're supposed to past me a list of residue numbers for the <residues> argument. Returning None."
        return None

    # use get_res_nums_within_radius to get all residue numbers
    residues_within_radius = []
    for res_num in residues:
        residues_within_radius.extend( get_res_nums_within_radius( res_num, input_pose, radius, include_res_num = include_res_nums ) )


    # get the set of the list and sort the residue numbers
    set_of_residues_within_radius = [ res for res in set( residues_within_radius ) ]

    # it is possible that there are still residues from <residues> in the list, so remove them one by one if desired
    if not include_res_nums:
        for res in residues:
            try:
                set_of_residues_within_radius.remove( res )
            except ValueError:
                pass

    # sort
    set_of_residues_within_radius.sort()

    return set_of_residues_within_radius



def spin_carbs_connected_to_prot( mm, input_pose, spin_using_ideal_omegas = True ):
    """
    The intent of this spin is to set the carbohydrate involved in the protein-carbohydrate connection into a reasonable starting position after the reset. This is because current LCM data found in the default.table was collected for surface glycans. These data are not reflective of the glycans found in the Ig system
    :param mm: MoveMap ( residues with BB set to True and that are carbohydrates are available to be reset )
    :param input_pose: Pose
    :param spin_using_ideal_omegas: bool( set omega1 (and omega2) to either 180, 60, or -60? If not, it would be those values +/- 0-20 ) Default = True
    :return: Pose
    """
    # imports
    from random import choice
    from rosetta.core.pose.carbohydrates import find_seqpos_of_saccharides_parent_residue, \
        get_reference_atoms, set_glycosidic_torsion
    from rosetta.core.id import omega_dihedral, omega2_dihedral
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg


    # copy the input pose
    pose = input_pose.clone()

    # use the MoveMap to see which residues have a parent connection to a protein
    # residues that are allowed to move are based on residues with BB set to True in the MoveMap and the residue is a carbohydrate
    carbohydrate_res_nums = [ res_num for res_num in range( 1, pose.n_residue() + 1 ) if mm.get_bb( res_num ) and pose.residue( res_num ).is_carbohydrate() ]

    # find the seqpos of the parent residue of each moveable residue and determine if it is a protein
    residues_to_spin = []
    for glyc_res_num in carbohydrate_res_nums:
        # skip residues that are of VariantType LOWER_TERMINUS as they, like Batman, do not have a parent
        if not pose.residue( glyc_res_num ).is_lower_terminus():
            # get residue number of this carbohydrate's parent and check if it is a protein residue
            if pose.residue( find_seqpos_of_saccharides_parent_residue( pose.residue( glyc_res_num ) ) ).is_protein():
                residues_to_spin.append( glyc_res_num )

    # reset values for omega 1 and omega 2
    # omega has stronger steric constraints
    possible_omega_values = [ 180, 60, -60 ]

    # for each carbohydrate residue connected to the protein
    for glyc_res_num in residues_to_spin:
        # bool( len([]) ) = False, so it will return False if there are no reference atoms (ie. it doesn't exist)
        # if an omega 1 exists (which it should...)
        if bool( len( get_reference_atoms( omega_dihedral, pose, glyc_res_num ) ) ):
            # get a base omega value to start at
            new_omega = choice( possible_omega_values )
            if not spin_using_ideal_omegas:
                # this specific format is pulled from rosetta.protocols.simple_moves:ShearMover::make_move
                # small_angle = 15 meaning the max you can move your angle in one direction ( + or - )
                # big_angle = 30 meaning the entire range available from your current ( + and - small_angle )
                new_omega = periodic_range( new_omega - 15.0 + rg().uniform() * 30.0, 360.0 )
            set_glycosidic_torsion( omega_dihedral, pose, glyc_res_num, new_omega )
        # if an omega2 exists
        if bool( len( get_reference_atoms( omega2_dihedral, pose, glyc_res_num ) ) ):
            # get a base omega value to start at
            new_omega2 = choice( possible_omega_values )
            if not spin_using_ideal_omegas:
                # this specific format is pulled from rosetta.protocols.simple_moves:ShearMover::make_move
                # small_angle = 15 meaning the max you can move your angle in one direction ( + or - )
                # big_angle = 30 meaning the entire range available from your current ( + and - small_angle )
                new_omega2 = periodic_range( new_omega2 - 15.0 + rg().uniform() * 30.0, 360.0 )
            set_glycosidic_torsion( omega2_dihedral, pose, glyc_res_num, new_omega2 )

    return pose



def glycosylate_working_pose( input_pose, glyco_file, glyco_sites ):
    """
    Glycosylate the <input_pose> using the <glyco_file> at the <glyco_sites> on the ND2 atom (assumes an ASN N-linked glycosylation, for now)
    Can work with PDB numbers (must have chain!) or Pose numbers for the glycosylation sites
    PDB number example: 123B. Or can give a list [ 123A, 123B ]
    Pose number example: 98. Or can give a list [ 98, 209 ]
    Assumes your glyco_file actually exists (doesn't use os.path.isfile)
    :param input_pose: Pose
    :param glyco_file: str( /path/to/glyco.iupac file )
    :param glyco_sites: PDB numbers list( [ 123A, 123B ] ) or Pose numbers list( [ 35, 78 ] )
    :return: glycosylated Pose
    """
    # imports
    import re, sys
    from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file


    # copy over the input_pose
    pose = input_pose.clone()

    # convert all the glyco_sites to Pose numbers first
    # this is because the PDB numbering changes as soon as you glycosylate the Pose (not the Pose numbering)
    # so you need the Pose numbering which remains absolute in order to properly glycosylate the Pose
    all_glyco_sites = {}
    for glyco_site in glyco_sites:
        # remove any whitespace
        glyco_site.replace( ' ', '' )

        # if the numbers are integers, then Pose numbers were passed
        try:
            #### POSE NUM ####
            all_glyco_sites[ glyco_site ] = int( glyco_site )

        # if the numbers didn't turn into integers, then a PDB number with a chain was passed
        except ValueError:
            # use a regular expression to match numbers and upper and lower case letters
            # should be something like 123A
            glyco_site_split = re.match(r"([0-9]+)([A-Za-z]+)", glyco_site )

            # if this resulted in a splitting, turn this PDB number and chain into a pose number
            #### PDB NUM ####
            if glyco_site_split is not None:
                # get the Pose number conversion from the PDB number <residue_num><residue_chain> format
                glyco_site_split = glyco_site_split.groups()
                glyco_site_pose_num = pose.pdb_info().pdb2pose( glyco_site_split[ 1 ], int( glyco_site_split[ 0 ] ) )
                # if glyco_site_pose_num came back as 0, then the chain and/or number doesn't exist in the Pose
                if glyco_site_pose_num == 0:
                    print "\nYou gave me an invalid PDB number as there is no to-Pose conversion. Does PDB number %s look right to you? Does this chain and number exist?\n" %glyco_site
                    sys.exit()
                # otherwise, this PDB number is a valid Pose number
                all_glyco_sites[ glyco_site ] = glyco_site_pose_num

            # if the regular expression returned None, then the PDB code wasn't in the right format that I asked for (number)(string)
            else:
                print "\nThere is something wrong with the glycosylation site you gave me. Does PDB number %s look right to you? Make sure it's in the order of <residue number><residue chain> with no spaces. Should look like 123A\n" %glyco_site
                sys.exit()

    # glycosylate all of the Pose-numbered glycosylation sites
    # it's stored as a dictionary as to relay to the user which glyco_site they originally gave messed up (if it does)
    for orig_glyco_site, glyco_site in all_glyco_sites.items():
        try:
            # since this should be a Pose number now, try it
            glycosylate_pose_by_file( pose, glyco_site, "ND2", glyco_file )
        except:
            # maybe the glycosylation site is in the right format but what they gave me is wrong?
            # maybe this is a Pose number but not a good glycosylation site?
            print "\nThere is something wrong with the glycosylation site you gave me. Does Pose number %s look right to you? Are you sure this residue can be glycosylated? Or maybe your glyco_file doesn't exist or work.\n" %orig_glyco_site
            sys.exit()

    return pose



def get_chains( input_pose, residue_range = None ):
    """
    Get a list of the chain ID's in the Pose. Can get the chain ID's associated with specific residues in <residue_range>
    :param input_pose: Pose
    :return: list( chains ) such as [ 'A', 'B', 'C' ]
    """
    chains = []

    # if no residue_range was given, replace it with all the residues in the input_pose
    if residue_range is None:
        residue_range = range( 1, input_pose.n_residue() + 1 )

    for res_num in residue_range:
        # example output of input_pose.pdb_info().pose2pdb( 45 ) is "245 B ", so strip whitespace and get the second item in the output
        chain_id = input_pose.pdb_info().pose2pdb( res_num ).strip().split( ' ' )[1]
        if chain_id not in chains:
            chains.append( chain_id )

    chains.sort()
    return chains
