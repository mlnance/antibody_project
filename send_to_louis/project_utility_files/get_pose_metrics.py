#!/usr/bin/python
__author__ = "morganlnance"



#####################################
#### DATA FROM LOW E NATIVE 3AY4 ####
#####################################

native_nhbonds = 445
working_nhbonds = 423  # low E native without Fc glycan




def main( in_working, working_info, in_native, native_info, in_sf, JUMP_NUM, decoy_num, dump_dir, util_dir, MC_acceptance_rate = None, native_constraint_file = None ):
    """
    Return a space-delimited string containing various pose metrics.
    :param in_working: decoy Pose()
    :param working_info: obj( class object that holds relevant chain and residue designations for the working pose )
    :param in_native: native Pose() 
    :param native_info: obj( class object that holds relevant chain and residue designations for the native pose )
    :param in_sf: ScoreFunction()
    :param JUMP_NUM: int( JUMP_NUM that defines the interface )
    :param decoy_num: int( the number of the decoy for use when dumping its Fc glycan )
    :param dump_dir: str( /path/to/dump_dir for the temp pdb files made. Files will be deleted )
    :param util_dir: str( /path/to/utilities directory
    :param MC_acceptance_rate: float( the MonteCarlo acceptance rate of your protocol, if relevant ). Default = None
    :param native_constraint_file: str( /path/to/corresponding native constraint file used to be used on native to get accurate delta_total_score )
    :return: str( pose metrics )
    """
    #################
    #### IMPORTS ####
    #################

    import sys
    try:
        sys.path.append( util_dir )
    except:
        pass

    # Rosetta functions
    from rosetta import Vector1, calc_interaction_energy
    from rosetta.protocols.simple_moves import ConstraintSetMover
    from toolbox import get_hbonds
    
    # Rosetta functions I wrote out
    from antibody_functions import calc_interface_sasa, decoy_to_native_res_map, \
        get_contact_map_between_range1_range2, get_contact_map_with_JUMP_NUM, \
        analyze_contact_map, calc_Fnat_with_contact_maps
    from pose_metrics_util import Fc_glycan_rmsd, Fc_glycan_hbonds, \
        pseudo_interface_energy_3ay4
    
    
    
    #############################
    #### METRIC CALCULATIONS ####
    #############################
    
    # holds all relevant metric data and corresponding label
    metric_data = []

    # copy the ScoreFunction and Poses passed in
    sf = in_sf.clone()
    working = in_working.clone()
    native = in_native.clone()

    ## delta total score calculation
    # give the native pose the constraint file passed, if any
    if native_constraint_file is not None:
        constraint_setter = ConstraintSetMover()
        constraint_setter.constraint_file( native_constraint_file )
        constraint_setter.apply( native )
    working_total_score = sf( working )
    native_total_score = sf( native )
    delta_total_score = working_total_score - native_total_score
    metric_data.append( "delta_total_score:" )
    metric_data.append( str( delta_total_score ) )
    

    # glycan RMSD calculation
    glycan_rmsd = Fc_glycan_rmsd( working, native, working_info.Fc_glycan_chains, native_info.native_Fc_glycan_chains, decoy_num, dump_dir )
    metric_data.append( "glycan_rmsd:" )
    metric_data.append( str( glycan_rmsd ) )


    ## pseudo-inferface energy 
    # ( full protein score - Fc-FcR glycan score [ except the short glycan away from interface ] )
    working_pseudo_interface_energy = pseudo_interface_energy_3ay4( working, sf, 
                                                                    native = False, 
                                                                    pmm = None )
    native_pseudo_interface_energy = pseudo_interface_energy_3ay4( native, sf, 
                                                                   native = True, 
                                                                   pmm = None )
    delta_pseudo_interface_energy = working_pseudo_interface_energy - native_pseudo_interface_energy
    metric_data.append( "pseudo_interface_energy:" )
    metric_data.append( str( working_pseudo_interface_energy ) )
    metric_data.append( "delta_pseudo_interface_energy:" )
    metric_data.append( str( delta_pseudo_interface_energy ) )


    # delta standard interaction energy ( across an interface defined by a JUMP number )
    working_interaction_energy = calc_interaction_energy( working, sf, Vector1( [ JUMP_NUM ] ) )
    native_interaction_energy = calc_interaction_energy( native, sf, Vector1( [ JUMP_NUM ] ) )
    delta_interaction_energy = working_interaction_energy - native_interaction_energy
    metric_data.append( "std_interface_interaction_energy:" )
    metric_data.append( str( working_interaction_energy ) )
    metric_data.append( "delta_std_interface_interaction_energy:" )
    metric_data.append( str( delta_interaction_energy ) )


    # delta hbonds in total complex
    working_hbonds = get_hbonds( working )
    native_hbonds = get_hbonds( native )
    delta_hbonds = working_hbonds.nhbonds() - native_hbonds.nhbonds()
    metric_data.append( "hbonds:" )
    metric_data.append( str( working_hbonds.nhbonds() ) )
    metric_data.append( "delta_hbonds:" )
    metric_data.append( str( delta_hbonds ) )


    #################
    # get the contact maps for the working and the native to use for other metric calculations
    Fc_glycan_to_Fc_protein_CUTOFF = 10
    working_Fc_glycan_to_Fc_protein_contact_map, working_Fc_glycan_to_Fc_protein_tot_contacts = get_contact_map_between_range1_range2( working_info.Fc_glycan_nums, 
                                                                                                                                       working_info.Fc_protein_nums, 
                                                                                                                                       working, 
                                                                                                                                       cutoff = Fc_glycan_to_Fc_protein_CUTOFF,
                                                                                                                                       return_more_info = True )
    native_Fc_glycan_to_Fc_protein_contact_map, native_Fc_glycan_to_Fc_protein_tot_contacts = get_contact_map_between_range1_range2( native_info.native_Fc_glycan_nums, 
                                                                                                                                     native_info.native_Fc_protein_nums, 
                                                                                                                                     native, 
                                                                                                                                     cutoff = Fc_glycan_to_Fc_protein_CUTOFF, 
                                                                                                                                     return_more_info = True )

    Fc_glycan_to_FcR_glycan_CUTOFF = 10
    working_Fc_glycan_to_FcR_glycan_contact_map, working_Fc_glycan_to_FcR_glycan_tot_contacts = get_contact_map_between_range1_range2( working_info.Fc_glycan_nums, 
                                                                                                                                       working_info.FcR_glycan_nums, 
                                                                                                                                       working, 
                                                                                                                                       cutoff = Fc_glycan_to_FcR_glycan_CUTOFF, 
                                                                                                                                       return_more_info = True )
    native_Fc_glycan_to_FcR_glycan_contact_map, native_Fc_glycan_to_FcR_glycan_tot_contacts = get_contact_map_between_range1_range2( native_info.native_Fc_glycan_nums, 
                                                                                                                                     native_info.native_FcR_glycan_nums, 
                                                                                                                                     native, 
                                                                                                                                     cutoff = Fc_glycan_to_FcR_glycan_CUTOFF, 
                                                                                                                                     return_more_info = True )
    #################


    # delta Fc-glycan to protein contacts
    Fc_glycan_to_Fc_protein_Fnat_contacts_recovered = calc_Fnat_with_contact_maps( working_Fc_glycan_to_Fc_protein_contact_map, 
                                                                                   working, 
                                                                                   native_Fc_glycan_to_Fc_protein_contact_map, 
                                                                                   native, 
                                                                                   decoy_to_native_res_map = decoy_to_native_res_map )
    delta_Fc_glycan_to_Fc_protein_tot_contacts = working_Fc_glycan_to_Fc_protein_tot_contacts - native_Fc_glycan_to_Fc_protein_tot_contacts
    metric_data.append( "Fc_glycan_to_Fc_protein_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( working_Fc_glycan_to_Fc_protein_tot_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_Fc_protein_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( delta_Fc_glycan_to_Fc_protein_tot_contacts ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_Fnat_contacts_recovered_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( Fc_glycan_to_Fc_protein_Fnat_contacts_recovered ) )


    # delta Fc-glycan to FcR glycan contacts
    Fc_glycan_to_FcR_glycan_Fnat_contacts_recovered = calc_Fnat_with_contact_maps( working_Fc_glycan_to_FcR_glycan_contact_map, 
                                                                                   working, 
                                                                                   native_Fc_glycan_to_FcR_glycan_contact_map, 
                                                                                   native, 
                                                                                   decoy_to_native_res_map = decoy_to_native_res_map )
    delta_Fc_glycan_to_FcR_glycan_tot_contacts = working_Fc_glycan_to_FcR_glycan_tot_contacts - native_Fc_glycan_to_FcR_glycan_tot_contacts
    metric_data.append( "Fc_glycan_to_FcR_glycan_tot_contacts_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    metric_data.append( str( working_Fc_glycan_to_FcR_glycan_tot_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_FcR_glycan_tot_contacts_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    metric_data.append( str( delta_Fc_glycan_to_FcR_glycan_tot_contacts ) )
    metric_data.append( "Fc_glycan_to_FcR_glycan_Fnat_contacts_recovered_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    metric_data.append( str( Fc_glycan_to_FcR_glycan_Fnat_contacts_recovered ) )

    
    # delta interface residue contacts - interface between Fc and FcR defined by JUMP_NUM 2
    intf_CUTOFF = 8
    working_intf_contact_map, working_intf_tot_contacts = get_contact_map_with_JUMP_NUM( JUMP_NUM, 
                                                                                         working, 
                                                                                         cutoff = intf_CUTOFF, 
                                                                                         return_more_info = True )
    native_intf_contact_map, native_intf_tot_contacts = get_contact_map_with_JUMP_NUM( JUMP_NUM, 
                                                                                       native, 
                                                                                       cutoff = intf_CUTOFF, 
                                                                                       return_more_info = True )
    intf_Fnat = calc_Fnat_with_contact_maps( working_intf_contact_map, 
                                             working, 
                                             native_intf_contact_map, 
                                             native, 
                                             decoy_to_native_res_map = decoy_to_native_res_map )
    delta_intf_tot_contacts = working_intf_tot_contacts - native_intf_tot_contacts
    metric_data.append( "interface_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_tot_contacts ) )
    metric_data.append( "delta_interface_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_tot_contacts ) )
    metric_data.append( "interface_Fnat_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnat ) )
    '''
    metric_data.append( "interface_pro_pro_fraction:" )
    metric_data.append( str( working_pro_pro_fraction ) )
    metric_data.append( "delta_interface_pro_pro_fraction:" )
    metric_data.append( str( delta_pro_pro_fraction ) )
    metric_data.append( "interface_pro_carb_fraction:" )
    metric_data.append( str( working_pro_carb_fraction ) )
    metric_data.append( "delta_interface_pro_carb_fraction:" )
    metric_data.append( str( delta_pro_carb_fraction ) )
    metric_data.append( "interface_carb_carb_fraction:" )
    metric_data.append( str( working_carb_carb_fraction ) )
    metric_data.append( "delta_interface_carb_carb_fraction:" )
    metric_data.append( str( delta_carb_carb_fraction ) )
    '''


    # delta interface sasa
    working_interface_sasa = calc_interface_sasa( working, JUMP_NUM )
    native_interface_sasa = calc_interface_sasa( native, JUMP_NUM )
    delta_interface_sasa = working_interface_sasa - native_interface_sasa
    metric_data.append( "interface_sasa:" )
    metric_data.append( str( working_interface_sasa ) )
    metric_data.append( "delta_interface_sasa:" )
    metric_data.append( str( delta_interface_sasa ) )

    
    # MonteCarlo acceptance rate - if relevant
    if MC_acceptance_rate is not None:
        metric_data.append( "MonteCarlo_acceptance_rate:" )
        metric_data.append( str( MC_acceptance_rate ) )

    # create metrics string
    metrics = ' '.join( metric_data )    
    
    return metrics
