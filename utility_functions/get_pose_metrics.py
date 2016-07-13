#!/usr/bin/python
__author__ = "morganlnance"



def main( working, native, sf, JUMP_NUM, working_Fc_glycan_chains, working_Fc_glycan_resnums, working_FcR_glycan_resnums, native_Fc_glycan_chains, native_Fc_glycan_resnums, native_FcR_glycan_resnums, decoy_num, dump_dir, MC_acceptance_rate = None ):
    """
    Return a space-delimited string containing various pose metrics.
    :param working: decoy Pose()
    :param native: native Pose()
    :param sf: ScoreFunction()
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
    from toolbox import get_hbonds
    
    # Rosetta functions I wrote out
    from antibody_functions import count_interface_residue_contacts, \
        calc_interface_sasa, count_residue_contacts_between_range1_range2
    from pose_metrics_util import Fc_glycan_rmsd, pseudo_interface_energy_3ay4
    
    
    
    #############################
    #### METRIC CALCULATIONS ####
    #############################
    
    # holds all relevant metric data and corresponding label
    metric_data = []
    

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


    ## delta Fc-glycan to protein contacts
    # get the protein residue numbers for the working and native pose
    working_protein_resnums = []
    for res in working:
        if not res.is_carbohydrate():
            working_protein_resnums.append( res.seqpos() )
    native_protein_resnums = []
    for res in native:
        if not res.is_carbohydrate():
            native_protein_resnums.append( res.seqpos() )

    working_Fc_glycan_to_protein_contacts, working_Fc_contact_map = count_residue_contacts_between_range1_range2( working_Fc_glycan_resnums, working_protein_resnums, working, cutoff = 10 )
    native_Fc_glycan_to_protein_contacts, native_Fc_contact_map = count_residue_contacts_between_range1_range2( native_Fc_glycan_resnums, native_protein_resnums, native, cutoff = 10 )
    delta_Fc_glycan_to_protein_contacts = working_Fc_glycan_to_protein_contacts - native_Fc_glycan_to_protein_contacts
    glycan_to_protein_contacts_Frecovered = round( ( float( working_Fc_glycan_to_protein_contacts ) / float( native_Fc_glycan_to_protein_contacts ) ) * 100, 2 )
    metric_data.append( "glycan_to_protein_contacts:" )
    metric_data.append( str( working_Fc_glycan_to_protein_contacts ) )
    metric_data.append( "delta_glycan_to_protein_contacts:" )
    metric_data.append( str( delta_Fc_glycan_to_protein_contacts ) )
    metric_data.append( "glycan_to_protein_contacts_Frecovered" )
    metric_data.append( str( glycan_to_protein_contacts_Frecovered ) )


    # delta Fc-glycan to FcR glycan contacts
    working_Fc_glycan_to_FcR_glycan_contacts, working_FcR_contact_map = count_residue_contacts_between_range1_range2( working_Fc_glycan_resnums, working_FcR_glycan_resnums, working, cutoff = 10 )
    native_Fc_glycan_to_FcR_glycan_contacts, native_FcR_contact_map = count_residue_contacts_between_range1_range2( native_Fc_glycan_resnums, native_FcR_glycan_resnums, native, cutoff = 10 )
    delta_Fc_glycan_to_FcR_glycan_contacts = working_Fc_glycan_to_FcR_glycan_contacts - native_Fc_glycan_to_FcR_glycan_contacts
    glycan_to_FcR_glycan_contacts_Frecovered = round( ( float( working_Fc_glycan_to_FcR_glycan_contacts ) / float( native_Fc_glycan_to_FcR_glycan_contacts ) ) * 100, 3 )
    metric_data.append( "glycan_to_FcR_glycan_contacts:" )
    metric_data.append( str( working_Fc_glycan_to_FcR_glycan_contacts ) )
    metric_data.append( "delta_glycan_to_FcR_glycan_contacts:" )
    metric_data.append( str( delta_Fc_glycan_to_FcR_glycan_contacts ) )
    metric_data.append( "glycan_to_FcR_glycan_contacts_Frecovered" )
    metric_data.append( str( glycan_to_FcR_glycan_contacts_Frecovered ) )

    
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
