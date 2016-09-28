#!/usr/bin/python
__author__="morganlnance"

'''
Various functions to be used in 3ay4 glycan modeling
Most compiled using the worker functions from antibody_functions
'''


def native_3ay4_Fc_glycan_LCM_reset( input_pose, residue_numbers, use_population_ideal_LCM_reset = False ):
    '''
    Reset the 3ay4 Fc glycan with freedom determined by <movemap_in>
    Default action is an LCM_reset with stdev +/1
    Uses data from the LCM default.table
    <use_population_ideal_LCM_reset> makes the reset only use the ideal values (with appropriate population weights). Essentially, stdev = 0
    :param input_pose: Pose
    :param residue_numbers: list( the Pose numbers for the residues of interest )
    :param use_population_ideal_LCM_reset: bool( use population ideals only? ) Default = False
    :param use_ideal_LCM_reset: bool( use only the ideal value from the highest population? ) Default = False
    '''
    # imports
    import sys
    from rosetta import MoveMap
    from rosetta.protocols.carbohydrates import LinkageConformerMover

    # check input arguments
    #if use_population_ideal_LCM_reset and use_ideal_LCM_reset:
    #    print "\nYou want me to use both the population ideal and the main ideal. This does not make sense to me. Ensure you are only setting one of these arguments to True. Exiting\n"
    #    sys.exit()

    # copy in the input_pose
    testing_pose = input_pose.clone()

    # for each residue in <residue_numbers>
    for res_num in residue_numbers:
        # make a BB MoveMap for this single residue
        res_mm = MoveMap()
        res_mm.set_bb( res_num, True )

        # create a LinkageConformerMover with the residue MoveMap
        lcm = LinkageConformerMover()
        lcm.set_movemap( res_mm )

        # if the user only wants to use ideals, but within the different population clusters
        if use_population_ideal_LCM_reset:
            # set_idealize_torsions uses ideal values instead of sampling from stdev
            lcm.set_idealize_torsions( True )
            # use_conformer_population_stats is True by default, but setting for clarity
            lcm.set_use_conformer_population_stats( True )
            # apply the LCM
            lcm.apply( testing_pose )

        # else, standard LCM reset using (by default) population data and a stdev of 1
        else:
            # setting these options for clarity
            lcm.set_use_conformer_population_stats( True )
            lcm.set_x_standard_deviations( 1 )
            # apply the LCM
            lcm.apply( testing_pose )

    return testing_pose



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

    # set the constraint from the file
    constraint_setter = ConstraintSetMover()
    constraint_setter.constraint_file( constraint_file )
    constraint_setter.apply( input_pose )

    return input_pose



def SugarSmallMover( seqpos, input_pose, angle_max, set_phi = True, set_psi = True, set_omega = True ):
    """
    Randomly resets the phi, psi, and omega values of the sugar residue <seqpos> in <input_pose> to old_value +/- angle_max/2
    Emulates the SmallMover but with the additional omega mover
    :param seqpos: int( the pose number for the residue )
    :param input_pose: Pose
    :param angle_max: int( or float( the max angle around which the phi/psi/omega could move ) )
    :param set_phi: bool( do you want to change the phi angle? ) Default = True
    :param set_psi: bool( do you want to change the psi angle? ) Default = True
    :param set_omega: bool( do you want to change the omega angle? ) Default = True
    :return: Pose
    """
    # imports
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg
    from rosetta.core.id import phi_dihedral, psi_dihedral, omega_dihedral
    from rosetta.core.pose.carbohydrates import set_glycosidic_torsion


    # copy the input pose
    pose = input_pose.clone()

    # from rosetta.protocols.simple_moves:BackboneMover.cc file for SmallMover
    big_angle = angle_max
    small_angle = big_angle / 2.0

    # get current phi, psi, and omega
    old_phi = pose.phi( seqpos )
    old_psi = pose.psi( seqpos )
    old_omega = pose.omega( seqpos )

    # get random values for phi, psi, and omega
    # this specific format is pulled from rosetta.protocols.simple_moves:ShearMover::make_move
    new_phi = periodic_range( old_phi - small_angle + rg().uniform() * big_angle, 360.0 )
    new_psi = periodic_range( old_psi - small_angle + rg().uniform() * big_angle, 360.0 )
    new_omega = periodic_range( old_omega - small_angle + rg().uniform() * big_angle, 360.0 )

    # set the new values
    if set_phi:
        set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
    if set_psi:
        set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
    if set_omega:
        set_glycosidic_torsion( omega_dihedral, pose, seqpos, new_omega )

    return pose



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
    # imports
    import sys
    from rosetta.core.scoring import score_type_from_name


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
