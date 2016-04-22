#!/usr/bin/python
__author__ = "morganlnance"



def get_pose_metrics( working, native, sf, JUMP_NUM ):
    """
    Return a space-delimited string containing various pose metrics.
    :param working: decoy Pose()
    :param native: native Pose()
    :param sf: ScoreFunction()
    :param JUMP_NUM: int( interface JUMP_NUM )
    :return: str( pose metrics )
    """
    #################
    #### IMPORTS ####
    #################

    from rosetta import Vector1, calc_interaction_energy, \
        calc_Fnat, calc_total_sasa
    from rosetta.core.scoring import non_peptide_heavy_atom_RMSD
    from toolbox import get_hbonds
    
    
    
    #############################
    #### METRIC CALCULATIONS ####
    #############################
    
    # holds all relevant metric data and corresponding label
    metric_data = []
    
    # ligand RMSD calculation
    ligand_rmsd = non_peptide_heavy_atom_RMSD( working, native )
    metric_data.append( "ligand_rmsd:" )
    metric_data.append( str( ligand_rmsd ) )
    
    # delta interaction energy
    working_interaction_energy = calc_interaction_energy( working, sf, Vector1( [ JUMP_NUM ] ) )
    native_interaction_energy = calc_interaction_energy( native, sf, Vector1( [ JUMP_NUM ] ) )
    delta_interaction_energy = working_interaction_energy - native_interaction_energy
    metric_data.append( "delta_interface_interaction_energy:" )
    metric_data.append( str( delta_interaction_energy ) )
    
    # fraction of native contacts
    #Fnat = calc_Fnat( working, native, sf, Vector1( [ JUMP_NUM ] ) )
    #metric_data.append( "Fnat:" )
    #metric_data.append( str( Fnat ) )
        
    # delta hbonds
    delta_hbonds = get_hbonds( working ).nhbonds() - get_hbonds( native ).nhbonds()
    metric_data.append( "delta_hbonds:" )
    metric_data.append( str( delta_hbonds ) )
    
    # create metrics string
    metrics = ' '.join( metric_data )
    
    
    return metrics




def read_fasc_file( fasc_filename, terms_to_collect ):
    """
    Read the <fasc_filename> and return a dictionary of the scoring data
    :param fasc_filename: str( /path/to/.fasc file
    :param terms_to_collect: list( of str( valid terms to collect from the .fasc file ) )
    :return: dict( .fasc scoring data )
    """
    #################
    #### IMPORTS ####
    #################

    import csv
    
    
    ##########################
    #### FASC FILE READER ####
    ##########################
    
    with open( fasc_filename, "rb" ) as f:
        # holds fasc data
        fasc_data = []
        
        for line in f.readlines():
            # remove carriage returns
            line = line.rstrip()
            
            # split the line on spaces
            space_split = line.split( ' ' )
            
            # ignore empty spaces
            for data in space_split:
                if data != ' ' and data != '':
                    # remove the trailing ':' character
                    fasc_data.append( data.replace( ':', '' ) )
    
    
    return fasc_data
