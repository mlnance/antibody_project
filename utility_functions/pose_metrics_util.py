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



class fasc_dict( dict ):
    """
    Allows you to create your own version of a dictionary so you can add attributes to it
    """
    pass



def read_fasc_file( fasc_filename ):
    """
    Read the <fasc_filename> and return a dictionary of the scoring data
    :param fasc_filename: str( /path/to/.fasc file
    :return: dict( .fasc scoring data )
    """
    #################
    #### IMPORTS ####
    #################

    import re
    
    try:
        import pandas
    except:
        import csv
    
    
    ##########################
    #### FASC FILE READER ####
    ##########################
     
    # create the dictionary that will hold the fasc_data_dict and its corresponding attributes
    fasc_data_dict = fasc_dict()
    fasc_data_dict.nstruct = None
    fasc_data_dict.pdb_name = None
    
    # open up the fasc_file
    with open( fasc_filename, "rb" ) as fh:
        for line in fh:
            # remove newline characters
            line = line.rstrip()
            
            # this should be the very first line of the .fasc file only
            if line.startswith( "pdb" ):
                # need to replace this specific space for clarity
                line = line.replace( "pdb name", "pdb_name" )
                
                # replace X amount of spaces that follow ':'
                arr = re.split( "[:\s]+", line)
            
    # turn the list of line data into an iterator
    iterator = iter( arr )
    
    # cycle through the first iterator to get the pdb_name and the nstruct
    # first iterator should always have this data
    for ii in iterator:
        key = str( ii )
        try:
            value = str( next( iterator ) )
            if key == "pdb_name":
                fasc_data_dict.pdb_name = value
            elif key == "nstruct":
                fasc_data_dict.nstruct = value
        except:
            pass
        
    # will hold valid decoy numbers in the fasc_data_dict object for easy access
    fasc_data_dict.decoy_nums = []
        
    # cycle through the .fasc file again to get the scoring terms and their values
    with open( fasc_filename, 'r' ) as fh:
        for line in fh:
            line = line.rstrip()

            # if the line is NOT the first line of the .fasc file
            if not line.startswith( "pdb" ):
                # split on multiple spaces after the ':' again
                arr = re.split( "[:\s]+", line )
                
                iterator = iter( arr )
                data_holder = {}
                for ii in iterator:
                    key = str( ii )
                    try:
                        # put the decoy data in the temporary holder
                        value = str( next( iterator ) )
                        data_holder[ key ] = value
                        
                        # pull out the filename decoy number
                        if key == "filename":
                            decoy_num = int( value.replace( ".pdb", '' ).split( '_' )[-1] )
                            fasc_data_dict.decoy_nums.append( decoy_num )
                    except:
                        pass
                    
                # add the decoy data to the fasc_data_dict object given the decoy number
                fasc_data_dict[ decoy_num ] = data_holder
                
    # just because I'm picky - sort the decoy_nums
    fasc_data_dict.decoy_nums.sort()
    
    # return the fasc data
    return fasc_data_dict



def get_score_term_from_fasc_data_dict( fasc_data_dict, score_term ):
    """
    After getting a <fasc_data_dict> from running read_fasc_file, pull out the <score_term> for each decoy in the .fasc file
    :param fasc_data_dict: dict( data dictionary from .fasc file )
    :param score_term: str( the score term you want back for each decoy
    :return: dict( decoy : score_term value )
    """
    # iterate through the dictionary skipping decoys that don't have that term
    score_term_dict = {}
    for decoy_num in fasc_data_dict.decoy_nums:
        try:
            score_term_dict[ decoy_num ] = float( fasc_data_dict[ decoy_num ][ score_term ] )
        # if it can only be a string
        except ValueError:
            score_term_dict[ decoy_num ] = fasc_data_dict[ decoy_num ][ score_term ]
        except:
            pass
        
    return score_term_dict