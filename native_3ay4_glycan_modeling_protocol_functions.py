#!/usr/bin/python
__author__="morganlnance"

'''
Various functions to be used in 3ay4 glycan modeling
Most compiled using the worker functions from antibody_functions
'''


def native_3ay4_Fc_glycan_LCM_reset( input_pose, residue_numbers, bb_freedom = True, chi_freedom = False, use_population_ideal_LCM_reset = False ):
    '''
    Reset the 3ay4 Fc glycan with freedom determined by <movemap_in>
    Default action is an LCM_reset with stdev +/1
    Uses data from the LCM default.table
    <use_population_ideal_LCM_reset> makes the reset only use the ideal values (with appropriate population weights). Essentially, stdev = 0
    :param input_pose: Pose
    :param residue_numbers: list( the Pose numbers for the residues of interest )
    :param bb_freedom: bool( allow bb freedom of each residue in <residue_numbers> in the MoveMap? ) Default = True
    :param chi_freedom: bool( allow chi freedom of each residue in <residue_numbers> in the MoveMap? ) Default = False
    :param use_population_ideal_LCM_reset: bool( use population ideals only? ) Default = False
    :param use_ideal_LCM_reset: bool( use only the ideal value from the highest population? ) Default = False
    '''
    # imports
    import sys
    from rosetta import MoveMap
    from rosetta.protocols.carbohydrates import LinkageConformerMover

    # check input arguments
    if use_population_ideal_LCM_reset and use_ideal_LCM_reset:
        print "\nYou want me to use both the population ideal and the main ideal. This does not make sense to me. Ensure you are only setting one of these arguments to True. Exiting\n"
        sys.exit()

    # copy in the input_pose
    testing_pose = input_pose.clone()

    # for each residue in <residue_numbers>
    for res_num in residue_numbers:
        # make a MoveMap for this single residue
        res_mm = MoveMap()

        # set bb and chi to True/False, as specified by the user
        # backbone
        if bb_freedom:
            res_mm.set_bb( res_num, True )
        # chi
        if chi_freedom:
            res_mm.set_chi( res_num, True )

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
