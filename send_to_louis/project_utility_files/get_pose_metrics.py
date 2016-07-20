#!/usr/bin/python
__author__ = "morganlnance"



def main( in_working, in_native, in_sf, JUMP_NUM, working_Fc_glycan_chains, working_Fc_glycan_resnums, working_FcR_glycan_resnums, native_Fc_glycan_chains, native_Fc_glycan_resnums, native_FcR_glycan_resnums, decoy_num, dump_dir, MC_acceptance_rate = None, constraint_file_used = None ):
    """
    Return a space-delimited string containing various pose metrics.
    :param in_working: decoy Pose()
    :param in_native: native Pose()
    :param in_sf: ScoreFunction()
    :param JUMP_NUM: int( JUMP_NUM that defines the interface )
    :param working_Fc_glycan_chains: list( the chain id's for the working Fc glycan ). Ex = [ 'H', 'I' ]
    :param working_Fc_glycan_resnums: list( pose residue numbers for the working Fc glycans ).
    :param working_FcR_glycan_resnums: list( pose residue numbers for the working FcR glycans ).
    :param native_Fc_glycan_chains: list( the chain id's for the native Fc glycan ). Ex = [ 'D', 'E' ]
    :param native_Fc_glycan_resnums: list( pose residue numbers for the native Fc glycans ).
    :param native_FcR_glycan_resnums: list( pose residue numbers for the native FcR glycans ).
    :param decoy_num: int( the number of the decoy for use when dumping its Fc glycan )
    :param dump_dir: str( /path/to/dump_dir for the temp pdb files made. Files will be deleted )
    :param MC_acceptance_rate: float( the MonteCarlo acceptance rate of your protocol, if relevant ). Default = None
    :param constraint_file_used: str( /path/to/constraint file used to be used on native to get accurate delta_total_score )
    :return: str( pose metrics )
    """
    #################
    #### IMPORTS ####
    #################

    import sys
    try:
        sys.path.append( "/home/mlnance/project_utility_files" )
        sys.path.append( "/Users/Research/pyrosetta_dir/project_utility_files" )
    except:
        pass

    # Rosetta functions
    from rosetta import Vector1, calc_interaction_energy
    from rosetta.protocols.simple_moves import ConstraintSetMover
    from toolbox import get_hbonds
    
    # Rosetta functions I wrote out
    from antibody_functions import count_interface_residue_contacts, \
        calc_interface_sasa, count_residue_contacts_between_range1_range2, \
        calc_res_Fnat_recovered_between_range1_range2, decoy_Fc_glycan, \
        decoy_Fc_protein, decoy_FcR_main_glycan
    from pose_metrics_util import Fc_glycan_rmsd, pseudo_interface_energy_3ay4
    
    
    
    #############################
    #### METRIC CALCULATIONS ####
    #############################
    
    # holds all relevant metric data and corresponding label
    metric_data = []

    # copy the scorefunction and Poses passed in
    sf = in_sf.clone()
    working = in_working.clone()
    native = in_native.clone()

    ## delta total score calculation
    # give the native pose the constraint file passed, if any
    if constraint_file is not None:
        constraint_setter = ConstraintSetMover()
        constraint_setter.constraint_file( input_args.constraint_file )
        constraint_setter.apply( working )
    working_total_score = sf( working )
    native_total_score = sf( native )
    

    # glycan RMSD calculation
    glycan_rmsd = Fc_glycan_rmsd( working, native, working_Fc_glycan_chains, native_Fc_glycan_chains, decoy_num, dump_dir )
    metric_data.append( "glycan_rmsd:" )
    metric_data.append( str( glycan_rmsd ) )


    # pseudo-inferface energy 
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


    # delta standard interaction energy
    working_interaction_energy = calc_interaction_energy( working, sf, Vector1( [ JUMP_NUM ] ) )
    native_interaction_energy = calc_interaction_energy( native, sf, Vector1( [ JUMP_NUM ] ) )
    delta_interaction_energy = working_interaction_energy - native_interaction_energy
    metric_data.append( "std_interface_interaction_energy:" )
    metric_data.append( str( working_interaction_energy ) )
    metric_data.append( "delta_std_interface_interaction_energy:" )
    metric_data.append( str( delta_interaction_energy ) )


    # delta hbonds
    working_hbonds = get_hbonds( working )
    native_hbonds = get_hbonds( native )
    delta_hbonds = working_hbonds.nhbonds() - native_hbonds.nhbonds()
    metric_data.append( "hbonds:" )
    metric_data.append( str( working_hbonds.nhbonds() ) )
    metric_data.append( "delta_hbonds:" )
    metric_data.append( str( delta_hbonds ) )


    # delta Fc-glycan to protein contacts
    working_Fc_glycan_to_protein_contacts, native_Fc_glycan_to_protein_contacts, Fc_glycan_to_protein_Fnat_recovered_contacts = calc_res_Fnat_recovered_between_range1_range2( working, decoy_Fc_glycan, decoy_Fc_protein, native, cutoff = 10, return_more_data = True )
    delta_Fc_glycan_to_protein_contacts = working_Fc_glycan_to_protein_contacts - native_Fc_glycan_to_protein_contacts
    metric_data.append( "glycan_to_protein_contacts:" )
    metric_data.append( str( working_Fc_glycan_to_protein_contacts ) )
    metric_data.append( "delta_glycan_to_protein_contacts:" )
    metric_data.append( str( delta_Fc_glycan_to_protein_contacts ) )
    metric_data.append( "Fc_glycan_to_protein_Fnat_recovered_contacts:" )
    metric_data.append( str( Fc_glycan_to_protein_Fnat_recovered_contacts ) )


    # delta Fc-glycan to FcR glycan contacts
    working_Fc_glycan_to_FcR_glycan_contacts, native_Fc_glycan_to_FcR_glycan_contacts, Fc_glycan_to_FcR_glycan_Fnat_recovered_contacts = calc_res_Fnat_recovered_between_range1_range2( working, decoy_Fc_glycan, decoy_FcR_main_glycan, native, cutoff = 10, return_more_data = True )
    delta_Fc_glycan_to_FcR_glycan_contacts = working_Fc_glycan_to_FcR_glycan_contacts - native_Fc_glycan_to_FcR_glycan_contacts
    metric_data.append( "Fc_glycan_to_FcR_glycan_contacts:" )
    metric_data.append( str( working_Fc_glycan_to_FcR_glycan_contacts ) )
    metric_data.append( "delta_Fc_glycan_to_FcR_glycan_contacts:" )
    metric_data.append( str( delta_Fc_glycan_to_FcR_glycan_contacts ) )
    metric_data.append( "Fc_glycan_to_FcR_glycan_Fnat_recovered_contacts:" )
    metric_data.append( str( Fc_glycan_to_FcR_glycan_Fnat_recovered_contacts ) )

    
    # delta interface residue contacts
    cutoff = 8
    working_intf_contacts, working_contact_list = count_interface_residue_contacts( JUMP_NUM, 
                                                                                    working, 
                                                                                    cutoff = cutoff )
    native_intf_contacts, native_contact_list = count_interface_residue_contacts( JUMP_NUM, 
                                                                                  native, 
                                                                                  cutoff = cutoff )
    delta_interface_res_contacts = working_intf_contacts - native_intf_contacts
    metric_data.append( "interface_res_contacts_%s_A:" %( str( cutoff ) ) )
    metric_data.append( str( working_intf_contacts ) )
    metric_data.append( "delta_interface_res_contacts_%s_A:" %( str( cutoff ) ) )
    metric_data.append( str( delta_interface_res_contacts ) )


    # delta Fc-glycan to protein contacts
    

    
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
