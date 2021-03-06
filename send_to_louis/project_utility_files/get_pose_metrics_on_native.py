#!/usr/bin/python
__author__ = "morganlnance"



def main( in_working, working_info, in_native, native_info, in_sf, JUMP_NUM, decoy_num, dump_dir, util_dir, MC_acceptance_rate = None, native_constraint_file = None, res_of_torsional_interest1 = None, res_of_torsional_interest2 = None ):
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
    :param native_constraint_file: str( /path/to/constraint file used to be used on native to get accurate delta_total_score )
    :param res_of_torsional_interest1: int( which residue do you want to get the phi, psi, and omega of? ) Default = None
    :param res_of_torsional_interest2: int( which residue do you want to get the phi, psi, and omega of? ) Default = None
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
    from antibody_functions import calc_interface_sasa, \
        get_contact_map_between_range1_range2, get_contact_map_with_JUMP_NUM, \
        analyze_contact_map, calc_Fnats_with_contact_maps, \
        get_scoretype_with_biggest_score_diff, get_res_with_biggest_score_diff
    from pose_metrics_util import Fc_glycan_metrics, \
        pseudo_interface_energy_3ay4, check_GlcNAc_to_Phe_contacts, \
        avg_dist_from_native_Fc_glycan_ring_atoms
    
    
    
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
    

    # delta hbonds in total complex
    working_hbonds = get_hbonds( working )
    native_hbonds = get_hbonds( native )
    delta_hbonds = working_hbonds.nhbonds() - native_hbonds.nhbonds()
    metric_data.append( "hbonds:" )
    metric_data.append( str( working_hbonds.nhbonds() ) )
    metric_data.append( "delta_hbonds:" )
    metric_data.append( str( delta_hbonds ) )
    

    # glycan RMSD, tot score, hbond, and sasa calculation
    working_Fc_glycan_data = Fc_glycan_metrics( working, native, working_info.native_Fc_glycan_chains, native_info.native_Fc_glycan_chains, sf, decoy_num, dump_dir )
    delta_Fc_glycan_tot_score = working_Fc_glycan_data.Fc_glycan_tot_score - working_Fc_glycan_data.native_Fc_glycan_tot_score
    delta_Fc_glycan_internal_hbonds = working_Fc_glycan_data.Fc_glycan_internal_hbonds - working_Fc_glycan_data.native_Fc_glycan_internal_hbonds
    delta_Fc_glycan_hbonds_contributed = working_Fc_glycan_data.Fc_glycan_hbonds_contributed - working_Fc_glycan_data.native_Fc_glycan_hbonds_contributed
    delta_Fc_glycan_sasa_contributed = working_Fc_glycan_data.Fc_glycan_sasa_contributed - working_Fc_glycan_data.native_Fc_glycan_sasa_contributed
    metric_data.append( "Fc_glycan_rmsd:" )
    metric_data.append( str( working_Fc_glycan_data.Fc_glycan_rmsd ) )
    metric_data.append( "Fc_glycan_tot_score:" )
    metric_data.append( str( working_Fc_glycan_data.Fc_glycan_tot_score ) )
    metric_data.append( "delta_Fc_glycan_tot_score:" )
    metric_data.append( str( delta_Fc_glycan_tot_score ) )
    metric_data.append( "Fc_glycan_internal_hbonds:" )
    metric_data.append( str( working_Fc_glycan_data.Fc_glycan_internal_hbonds ) )
    metric_data.append( "delta_Fc_glycan_internal_hbonds:" )
    metric_data.append( str( delta_Fc_glycan_internal_hbonds ) )
    metric_data.append( "Fc_glycan_hbonds_contributed:" )
    metric_data.append( str( working_Fc_glycan_data.Fc_glycan_hbonds_contributed ) )
    metric_data.append( "delta_Fc_glycan_hbonds_contributed:" )
    metric_data.append( str( delta_Fc_glycan_hbonds_contributed ) )
    metric_data.append( "Fc_glycan_sasa_contributed:" )
    metric_data.append( str( working_Fc_glycan_data.Fc_glycan_sasa_contributed ) )
    metric_data.append( "delta_Fc_glycan_sasa_contributed:" )
    metric_data.append( str( delta_Fc_glycan_sasa_contributed ) )


    # average ring atom distances (not rmsd as there is no alignment) of decoy to native
    for ring_residue in working_info.native_Fc_glycan_nums_except_core_GlcNAc:
        # decoy numbering is the same as native, so just use the same number
        ring_atoms_avg_dist = avg_dist_from_native_Fc_glycan_ring_atoms( working, ring_residue, native, ring_residue )
        metric_data.append( "ring_atoms_avg_dist_res_%s:" %str( ''.join( working.pdb_info().pose2pdb( ring_residue ).strip().split( ' ' ) ) ) )
        metric_data.append( str( round( ring_atoms_avg_dist, 2 ) ) )
    

    ## pseudo-inferface energy 
    # ( full protein score - Fc-FcR main glycan score
    working_pseudo_interface_energy = pseudo_interface_energy_3ay4( working, sf, 
                                                                    native = True, 
                                                                    pmm = None )
    native_pseudo_interface_energy = pseudo_interface_energy_3ay4( native, sf, 
                                                                   native = True, 
                                                                   pmm = None )
    delta_pseudo_interface_energy = working_pseudo_interface_energy - native_pseudo_interface_energy
    metric_data.append( "pseudo_interface_energy:" )
    metric_data.append( str( working_pseudo_interface_energy ) )
    metric_data.append( "delta_pseudo_interface_energy:" )
    metric_data.append( str( delta_pseudo_interface_energy ) )


    # determine which ScoreType contributes the biggest difference in total score compared to native
    working_total_score_data = get_scoretype_with_biggest_score_diff( working, native, sf )
    metric_data.append( "delta_biggest_pos_score_diff_by_scoretype:" )
    metric_data.append( str( working_total_score_data.biggest_pos_delta_score_by_scoretype ) )
    metric_data.append( "biggest_pos_score_diff_scoretype:" )
    metric_data.append( str( working_total_score_data.str_score_type_pos ) )
    metric_data.append( "delta_biggest_neg_score_diff_by_scoretype:" )
    metric_data.append( str( working_total_score_data.biggest_neg_delta_score_by_scoretype ) )
    metric_data.append( "biggest_neg_score_diff_scoretype:" )
    metric_data.append( str( working_total_score_data.str_score_type_neg ) )


    # determine the residue with the biggest score difference compared to native
    working_residue_score_data = get_res_with_biggest_score_diff( working, native, sf )
    metric_data.append( "delta_res_biggest_score_diff_tot_score:" )
    metric_data.append( str( working_residue_score_data.biggest_delta_score ) )
    metric_data.append( "decoy_res_num:" )
    metric_data.append( str( ''.join( working.pdb_info().pose2pdb( working_residue_score_data.decoy_num ).strip().split( ' ' ) ) ) )
    

    # delta standard interaction energy ( across an interface defined by a JUMP number )
    working_interaction_energy = calc_interaction_energy( working, sf, Vector1( [ JUMP_NUM ] ) )
    native_interaction_energy = calc_interaction_energy( native, sf, Vector1( [ JUMP_NUM ] ) )
    delta_interaction_energy = working_interaction_energy - native_interaction_energy
    metric_data.append( "std_interface_interaction_energy:" )
    metric_data.append( str( working_interaction_energy ) )
    metric_data.append( "delta_std_interface_interaction_energy:" )
    metric_data.append( str( delta_interaction_energy ) )
    

    # check if the GlcNAc above the Gal residue contacts the Phe residue within 5 Angstroms ( 4.68 contact distance in native )
    GlcNAc_to_Phe_cutoff = 5
    working_GlcNAc_to_Phe_contacts = check_GlcNAc_to_Phe_contacts( working, GlcNAc_to_Phe_cutoff, native = True )
    metric_data.append( "GlcNAc_to_its_Phe_contacts_%sA:" %( str( GlcNAc_to_Phe_cutoff ) ) )
    metric_data.append( str( working_GlcNAc_to_Phe_contacts ) )
    

    #################
    # get the contact maps for the working and the native to use for other metric calculations
    Fc_glycan_to_Fc_protein_CUTOFF = 8
    working_Fc_glycan_to_Fc_protein_contact_map, working_Fc_glycan_to_Fc_protein_tot_contacts = get_contact_map_between_range1_range2( working_info.native_Fc_glycan_nums_except_core_GlcNAc, 
                                                                                                                                       working_info.native_Fc_protein_nums, 
                                                                                                                                       working, 
                                                                                                                                       cutoff = Fc_glycan_to_Fc_protein_CUTOFF,
                                                                                                                                       return_more_info = True )
    native_Fc_glycan_to_Fc_protein_contact_map, native_Fc_glycan_to_Fc_protein_tot_contacts = get_contact_map_between_range1_range2( native_info.native_Fc_glycan_nums_except_core_GlcNAc, 
                                                                                                                                     native_info.native_Fc_protein_nums, 
                                                                                                                                     native, 
                                                                                                                                     cutoff = Fc_glycan_to_Fc_protein_CUTOFF, 
                                                                                                                                     return_more_info = True )

    Fc_glycan_to_FcR_glycan_CUTOFF = 8
    working_Fc_glycan_to_FcR_glycan_contact_map, working_Fc_glycan_to_FcR_glycan_tot_contacts = get_contact_map_between_range1_range2( working_info.native_Fc_glycan_nums_except_core_GlcNAc, 
                                                                                                                                       working_info.native_FcR_glycan_nums, 
                                                                                                                                       working, 
                                                                                                                                       cutoff = Fc_glycan_to_FcR_glycan_CUTOFF, 
                                                                                                                                       return_more_info = True )
    native_Fc_glycan_to_FcR_glycan_contact_map, native_Fc_glycan_to_FcR_glycan_tot_contacts = get_contact_map_between_range1_range2( native_info.native_Fc_glycan_nums_except_core_GlcNAc, 
                                                                                                                                     native_info.native_FcR_glycan_nums, 
                                                                                                                                     native, 
                                                                                                                                     cutoff = Fc_glycan_to_FcR_glycan_CUTOFF, 
                                                                                                                                     return_more_info = True )

    '''
    intf_CUTOFF = 8
    working_intf_contact_map, working_intf_tot_contacts = get_contact_map_with_JUMP_NUM( JUMP_NUM, 
                                                                                         working, 
                                                                                         cutoff = intf_CUTOFF, 
                                                                                         return_more_info = True )
    native_intf_contact_map, native_intf_tot_contacts = get_contact_map_with_JUMP_NUM( JUMP_NUM, 
                                                                                       native, 
                                                                                       cutoff = intf_CUTOFF, 
                                                                                       return_more_info = True )
    '''

    #################

    ## analyze the contact maps
    # Fc glycan to Fc protein contact map analysis
    working_Fc_glycan_to_Fc_protein_data_holder = analyze_contact_map( working_Fc_glycan_to_Fc_protein_contact_map, working )
    native_Fc_glycan_to_Fc_protein_data_holder = analyze_contact_map( native_Fc_glycan_to_Fc_protein_contact_map, native )

    # Fc glycan to protein Fnat calculations
    Fc_glycan_to_Fc_protein_Fnats_data_holder = calc_Fnats_with_contact_maps( working_Fc_glycan_to_Fc_protein_contact_map, 
                                                                              working, 
                                                                              native_Fc_glycan_to_Fc_protein_contact_map, 
                                                                              native )

    delta_Fc_glycan_to_Fc_protein_tot_contacts = working_Fc_glycan_to_Fc_protein_tot_contacts - native_Fc_glycan_to_Fc_protein_tot_contacts
    delta_Fc_glycan_to_Fc_protein_carb_to_polar_contacts = working_Fc_glycan_to_Fc_protein_data_holder.carb_to_polar_contacts - native_Fc_glycan_to_Fc_protein_data_holder.carb_to_polar_contacts
    delta_Fc_glycan_to_Fc_protein_carb_to_nonpolar_contacts = working_Fc_glycan_to_Fc_protein_data_holder.carb_to_nonpolar_contacts - native_Fc_glycan_to_Fc_protein_data_holder.carb_to_nonpolar_contacts
    delta_Fc_glycan_to_Fc_protein_carb_to_aromatic_contacts = working_Fc_glycan_to_Fc_protein_data_holder.carb_to_aromatic_contacts - native_Fc_glycan_to_Fc_protein_data_holder.carb_to_aromatic_contacts
    #delta_Fc_glycan_to_Fc_protein_contact_distance_avg = working_Fc_glycan_to_Fc_protein_data_holder.contact_distance_avg - native_Fc_glycan_to_Fc_protein_data_holder.contact_distance_avg
    #delta_Fc_glycan_to_Fc_protein_contact_distance_max = working_Fc_glycan_to_Fc_protein_data_holder.contact_distance_max - native_Fc_glycan_to_Fc_protein_data_holder.contact_distance_max
    #delta_Fc_glycan_to_Fc_protein_contact_distance_min = working_Fc_glycan_to_Fc_protein_data_holder.contact_distance_min - native_Fc_glycan_to_Fc_protein_data_holder.contact_distance_min
    metric_data.append( "Fc_glycan_to_Fc_protein_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( working_Fc_glycan_to_Fc_protein_tot_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_Fc_protein_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( delta_Fc_glycan_to_Fc_protein_tot_contacts ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_Fnat_tot_contacts_recovered_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( Fc_glycan_to_Fc_protein_Fnats_data_holder.Fnat_tot_contacts_recovered ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_carb_to_polar_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( working_Fc_glycan_to_Fc_protein_data_holder.carb_to_polar_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_Fc_protein_carb_to_polar_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( delta_Fc_glycan_to_Fc_protein_carb_to_polar_contacts ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_Fnat_carb_to_polar_contacts_recovered_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( Fc_glycan_to_Fc_protein_Fnats_data_holder.Fnat_carb_to_polar_contacts_recovered ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_carb_to_nonpolar_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( working_Fc_glycan_to_Fc_protein_data_holder.carb_to_nonpolar_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_Fc_protein_carb_to_nonpolar_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( delta_Fc_glycan_to_Fc_protein_carb_to_nonpolar_contacts ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_Fnat_carb_to_nonpolar_contacts_recovered_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( Fc_glycan_to_Fc_protein_Fnats_data_holder.Fnat_carb_to_nonpolar_contacts_recovered ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_carb_to_aromatic_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( working_Fc_glycan_to_Fc_protein_data_holder.carb_to_aromatic_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_Fc_protein_carb_to_aromatic_tot_contacts_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( delta_Fc_glycan_to_Fc_protein_carb_to_aromatic_contacts ) )
    metric_data.append( "Fc_glycan_to_Fc_protein_Fnat_carb_to_aromatic_contacts_recovered_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    metric_data.append( str( Fc_glycan_to_Fc_protein_Fnats_data_holder.Fnat_carb_to_aromatic_contacts_recovered ) )
    #metric_data.append( "Fc_glycan_to_Fc_protein_contact_distance_avg_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    #metric_data.append( str( working_Fc_glycan_to_Fc_protein_data_holder.contact_distance_avg ) )
    #metric_data.append( "delta_Fc_glycan_to_Fc_protein_contact_distance_avg_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    #metric_data.append( str( delta_Fc_glycan_to_Fc_protein_contact_distance_avg ) )
    #metric_data.append( "Fc_glycan_to_Fc_protein_contact_distance_max_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    #metric_data.append( str( working_Fc_glycan_to_Fc_protein_data_holder.contact_distance_max ) )
    #metric_data.append( "delta_Fc_glycan_to_Fc_protein_contact_distance_max_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    #metric_data.append( str( delta_Fc_glycan_to_Fc_protein_contact_distance_max ) )
    #metric_data.append( "Fc_glycan_to_Fc_protein_contact_distance_min_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    #metric_data.append( str( working_Fc_glycan_to_Fc_protein_data_holder.contact_distance_min ) )
    #metric_data.append( "delta_Fc_glycan_to_Fc_protein_contact_distance_min_%sA:" %( str( Fc_glycan_to_Fc_protein_CUTOFF ) ) )
    #metric_data.append( str( delta_Fc_glycan_to_Fc_protein_contact_distance_min ) )
    
    # Fc glycan to FcR glycan contact map analysis
    working_Fc_glycan_to_FcR_glycan_data_holder = analyze_contact_map( working_Fc_glycan_to_FcR_glycan_contact_map, working )
    native_Fc_glycan_to_FcR_glycan_data_holder = analyze_contact_map( native_Fc_glycan_to_FcR_glycan_contact_map, native )

    # Fc glycan to FcR glycan Fnat calculations
    Fc_glycan_to_FcR_glycan_Fnats_data_holder = calc_Fnats_with_contact_maps( working_Fc_glycan_to_FcR_glycan_contact_map, 
                                                                              working, 
                                                                              native_Fc_glycan_to_FcR_glycan_contact_map, 
                                                                              native )

    delta_Fc_glycan_to_FcR_glycan_tot_contacts = working_Fc_glycan_to_FcR_glycan_tot_contacts - native_Fc_glycan_to_FcR_glycan_tot_contacts
    #delta_Fc_glycan_to_FcR_glycan_contact_distance_avg = working_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_avg - native_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_avg
    #delta_Fc_glycan_to_FcR_glycan_contact_distance_max = working_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_max - native_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_max
    #delta_Fc_glycan_to_FcR_glycan_contact_distance_min = working_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_min - native_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_min
    metric_data.append( "Fc_glycan_to_FcR_glycan_tot_contacts_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    metric_data.append( str( working_Fc_glycan_to_FcR_glycan_tot_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_FcR_glycan_tot_contacts_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    metric_data.append( str( delta_Fc_glycan_to_FcR_glycan_tot_contacts ) )
    metric_data.append( "Fc_glycan_to_FcR_glycan_Fnat_tot_contacts_recovered_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    metric_data.append( str( Fc_glycan_to_FcR_glycan_Fnats_data_holder.Fnat_tot_contacts_recovered ) )
    #metric_data.append( "Fc_glycan_to_FcR_glycan_contact_distance_avg_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    #metric_data.append( str( working_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_avg ) )
    #metric_data.append( "delta_Fc_glycan_to_FcR_glycan_contact_distance_avg_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    #metric_data.append( str( delta_Fc_glycan_to_FcR_glycan_contact_distance_avg ) )
    #metric_data.append( "Fc_glycan_to_FcR_glycan_contact_distance_max_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    #metric_data.append( str( working_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_max ) )
    #metric_data.append( "delta_Fc_glycan_to_FcR_glycan_contact_distance_max_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    #metric_data.append( str( delta_Fc_glycan_to_FcR_glycan_contact_distance_max ) )
    #metric_data.append( "Fc_glycan_to_FcR_glycan_contact_distance_min_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    #metric_data.append( str( working_Fc_glycan_to_FcR_glycan_data_holder.contact_distance_min ) )
    #metric_data.append( "delta_Fc_glycan_to_FcR_glycan_contact_distance_min_%sA:" %( str( Fc_glycan_to_FcR_glycan_CUTOFF ) ) )
    #metric_data.append( str( delta_Fc_glycan_to_FcR_glycan_contact_distance_min ) )


    '''
    # interface residue contact map analysis
    working_intf_data_holder = analyze_contact_map( working_intf_contact_map, working )
    native_intf_data_holder = analyze_contact_map( native_intf_contact_map, native )

    # interface residue Fnat calculations - interface between Fc and FcR defined by JUMP_NUM 2
    intf_Fnats_data_holder = calc_Fnats_with_contact_maps( working_intf_contact_map, 
                                                           working, 
                                                           native_intf_contact_map, 
                                                           native )

    delta_intf_tot_contacts = working_intf_tot_contacts - native_intf_tot_contacts
    delta_intf_pro_pro_contacts = working_intf_data_holder.pro_pro_contacts - native_intf_data_holder.pro_pro_contacts
    delta_intf_pro_carb_contacts = working_intf_data_holder.pro_carb_contacts - native_intf_data_holder.pro_carb_contacts
    delta_intf_carb_carb_contacts = working_intf_data_holder.carb_carb_contacts - native_intf_data_holder.carb_carb_contacts
    delta_intf_carb_to_polar_contacts = working_intf_data_holder.carb_to_polar_contacts - native_intf_data_holder.carb_to_polar_contacts
    delta_intf_carb_to_nonpolar_contacts = working_intf_data_holder.carb_to_nonpolar_contacts - native_intf_data_holder.carb_to_nonpolar_contacts
    delta_intf_carb_to_aromatic_contacts = working_intf_data_holder.carb_to_aromatic_contacts - native_intf_data_holder.carb_to_aromatic_contacts
    #delta_intf_contact_distance_avg = working_intf_data_holder.contact_distance_avg - native_intf_data_holder.contact_distance_avg
    #delta_intf_contact_distance_max = working_intf_data_holder.contact_distance_max - native_intf_data_holder.contact_distance_max
    #delta_intf_contact_distance_min = working_intf_data_holder.contact_distance_min - native_intf_data_holder.contact_distance_min
    metric_data.append( "interface_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_tot_contacts ) )
    metric_data.append( "delta_interface_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_tot_contacts ) )
    metric_data.append( "interface_Fnat_tot_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnats_data_holder.Fnat_tot_contacts_recovered ) )
    metric_data.append( "intf_pro_pro_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_data_holder.pro_pro_contacts ) )
    metric_data.append( "delta_intf_pro_pro_tot_contacts%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_pro_pro_contacts ) )
    metric_data.append( "interface_Fnat_pro_pro_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnats_data_holder.Fnat_pro_pro_contacts_recovered ) )
    metric_data.append( "intf_pro_carb_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_data_holder.pro_carb_contacts ) )
    metric_data.append( "delta_intf_pro_carb_tot_contacts%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_pro_carb_contacts ) )
    metric_data.append( "interface_Fnat_pro_carb_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnats_data_holder.Fnat_pro_carb_contacts_recovered ) )
    metric_data.append( "intf_carb_carb_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_data_holder.carb_carb_contacts ) )
    metric_data.append( "delta_intf_carb_carb_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_carb_carb_contacts ) )
    metric_data.append( "interface_Fnat_carb_carb_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnats_data_holder.Fnat_carb_carb_contacts_recovered ) )
    metric_data.append( "intf_carb_to_polar_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_data_holder.carb_to_polar_contacts ) )
    metric_data.append( "delta_intf_carb_to_polar_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_carb_to_polar_contacts ) )
    metric_data.append( "interface_Fnat_carb_to_polar_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnats_data_holder.Fnat_carb_to_polar_contacts_recovered ) )
    metric_data.append( "intf_carb_to_nonpolar_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_data_holder.carb_to_nonpolar_contacts ) )
    metric_data.append( "delta_intf_carb_to_nonpolar_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_carb_to_nonpolar_contacts ) )
    metric_data.append( "interface_Fnat_carb_to_nonpolar_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnats_data_holder.Fnat_carb_to_nonpolar_contacts_recovered ) )
    metric_data.append( "intf_carb_to_aromatic_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( working_intf_data_holder.carb_to_aromatic_contacts ) )
    metric_data.append( "delta_intf_carb_to_aromatic_tot_contacts_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( delta_intf_carb_to_aromatic_contacts ) )
    metric_data.append( "interface_Fnat_carb_to_aromatic_contacts_recovered_%sA:" %( str( intf_CUTOFF ) ) )
    metric_data.append( str( intf_Fnats_data_holder.Fnat_carb_to_aromatic_contacts_recovered ) )
    #metric_data.append( "intf_contact_distance_avg_%sA:" %( str( intf_CUTOFF ) ) )
    #metric_data.append( str( working_intf_data_holder.contact_distance_avg ) )
    #metric_data.append( "delta_intf_contact_distance_avg_%sA:" %( str( intf_CUTOFF ) ) )
    #metric_data.append( str( delta_intf_contact_distance_avg ) )
    #metric_data.append( "intf_contact_distance_max_%sA:" %( str( intf_CUTOFF ) ) )
    #metric_data.append( str( working_intf_data_holder.contact_distance_max ) )
    #metric_data.append( "delta_intf_contact_distance_max_%sA:" %( str( intf_CUTOFF ) ) )
    #metric_data.append( str( delta_intf_contact_distance_max ) )
    #metric_data.append( "intf_contact_distance_min_%sA:" %( str( intf_CUTOFF ) ) )
    #metric_data.append( str( working_intf_data_holder.contact_distance_min ) )
    #metric_data.append( "delta_intf_contact_distance_min_%sA:" %( str( intf_CUTOFF ) ) )
    #metric_data.append( str( delta_intf_contact_distance_min ) )
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

    if res_of_torsional_interest1 is not None:
        metric_data.append( "APhi" )
        metric_data.append( str( working.phi( res_of_torsional_interest1 ) ) )
        metric_data.append( "APsi" )
        metric_data.append( str( working.psi( res_of_torsional_interest1 ) ) )
        metric_data.append( "AOmega" )
        metric_data.append( str( working.omega( res_of_torsional_interest1 ) ) )
    if res_of_torsional_interest2 is not None:
        metric_data.append( "BPhi" )
        metric_data.append( str( working.phi( res_of_torsional_interest2 ) ) )
        metric_data.append( "BPsi" )
        metric_data.append( str( working.psi( res_of_torsional_interest2 ) ) )
        metric_data.append( "BOmega" )
        metric_data.append( str( working.omega( res_of_torsional_interest2 ) ) )

    # create metrics string
    metrics = ' '.join( metric_data )    
    
    return metrics
