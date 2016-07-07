#!/usr/bin/python
__author__ = "morganlnance"



def pseudo_interface_score_3ay4( pose, sf, native = False, pmm = None ):
    """
    Attempts to get pseudo-interface score of a glycosylated 3ay4 decoy
    Lots of hard coding here - works on a decoy pose as Rosetta renumbers the Pose a bit
    Makes the two ASN connections to the Fc A and B glycans JUMPs instead of chemical EDGEs
    :return: 
    """
    from rosetta import FoldTree, Pose
    from rosetta.numeric import xyzVector_Real

    # save the original FoldTree
    orig_ft = pose.fold_tree()

    # get the score of the whole complex
    start_score = sf( pose )

    # hard code the new FoldTree specific to a glycosylated decoy of 3ay4
    if not native:
        ft = FoldTree()
        ft.add_edge( 1, 215, -1 )
        ft.add_edge( 1, 216, 1 )
        ft.add_edge( 216, 431, -1 )
        ft.add_edge( 1, 432, 2 )
        ft.add_edge( 432, 591, -1 )
        ft.add_edge( 579, 592, "ND2", "C1" )
        ft.add_edge( 592, 596, -1 )
        ft.add_edge( 594, 597, "O6", "C1" )
        ft.add_edge( 597, 598, -1 )
        ft.add_edge( 592, 599, "O6", "C1" )
        ft.add_edge( 462, 600, "ND2", "C1" )
        ft.add_edge( 600, 602, -1 )
        ft.add_edge( 69, 603, 3 )
        ft.add_edge( 603, 607, -1 )
        ft.add_edge( 605, 608, "O6", "C1" )
        ft.add_edge( 608, 610, -1 )
        ft.add_edge( 281, 611, 4 )
        ft.add_edge( 611, 615, -1 )
        ft.add_edge( 613, 616, "O6", "C1" )
        ft.add_edge( 616, 618, -1 )

    # hard code the new FoldTree specific to the native 3ay4
    else:
        ft = FoldTree()
        ft.add_edge( 1, 215, -1 )
        ft.add_edge( 69, 216, 1 )
        ft.add_edge( 216, 220, -1 )
        ft.add_edge( 218, 221, "O6", "C1" )
        ft.add_edge( 221, 223, -1 )
        ft.add_edge( 1, 224, 2 )
        ft.add_edge( 224, 439, -1 )
        ft.add_edge( 292, 440, 3 )
        ft.add_edge( 440, 444, -1 )
        ft.add_edge( 442, 445, "O6", "C1" )
        ft.add_edge( 445, 447, -1 )
        ft.add_edge( 1, 448, 4 )
        ft.add_edge( 448, 607, -1 )
        ft.add_edge( 595, 608, "ND2", "C1" )
        ft.add_edge( 608, 612, -1 )
        ft.add_edge( 610, 613, "O6", "C1" )
        ft.add_edge( 613, 614, -1 )
        ft.add_edge( 608, 615, "O6", "C1" )
        ft.add_edge( 478, 616, "ND2", "C1" )
        ft.add_edge( 616, 618, -1 )
        
    # make a temporary Pose and give the new FoldTree to it
    #try:
    #    pmm.keep_history( True )
    #except:
    #    pass
    temp_pose = Pose()
    temp_pose.assign( pose )
    #try:
    #    pmm.apply( temp_pose )
    #except:
    #    pass
    temp_pose.fold_tree( ft )

    # split apart the two Fc sugars one-by-one
    # if decoy structure
    if not native:
        jump = temp_pose.jump( 3 ) # sugar A
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 3, jump )
        #try:
        #    pmm.apply( temp_pose )
        #except:
        #    pass
        jump = temp_pose.jump( 4 ) # sugar B
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 4, jump )
        #try:
        #    pmm.apply( temp_pose )
        #except:
        #    pass

    # else native structure
    else:
        jump = temp_pose.jump( 1 ) # sugar A
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 1, jump )
        #try:
        #    pmm.apply( temp_pose )
        #except:
        #    pass
        jump = temp_pose.jump( 3 ) # sugar B
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 3, jump )
        #try:
        #    pmm.apply( temp_pose )
        #except:
        #    pass        

    # score the split-apart Pose
    split_score = sf( temp_pose )

    # get the pseduo-interface score
    pseudo_interface_score = start_score - split_score

    return pseudo_interface_score



def get_pose_metrics( working, native, sf, JUMP_NUM, working_Fc_glycan_chains, native_Fc_glycan_chains, decoy_num, dump_dir, MC_acceptance_rate = None ):
    """
    Return a space-delimited string containing various pose metrics.
    :param working: decoy Pose()
    :param native: native Pose()
    :param sf: ScoreFunction()
    :param JUMP_NUM: int( JUMP_NUM that defines the interface )
    :param working_Fc_glycan_chains: list( the chain id's for the working Fc glycan ). Ex = [ 'H', 'I' ]
    :param native_Fc_glycan_chains: list( the chain id's for the native Fc glycan ). Ex = [ 'D', 'E' ]
    :param decoy_num: int( the number of the decoy for use when dumping its Fc glycan )
    :param dump_dir: str( /path/to/dump_dir for the temp pdb files made. Files will be deleted )
    :param MC_acceptance_rate: float( the MonteCarlo acceptance rate of your protocol, if relevant ). Default = None
    :return: str( pose metrics )
    """
    #################
    #### IMPORTS ####
    #################
    
    # Rosetta functions
    from rosetta import Pose, Vector1, calc_interaction_energy
    from rosetta.core.scoring import non_peptide_heavy_atom_RMSD
    from toolbox import get_hbonds
    
    # Rosetta functions I wrote out
    from antibody_functions import count_interface_residue_contacts, \
        calc_interface_sasa, load_pose
    
    # utility functions
    import os
    from util import dump_pdb_by_chain
    
    
    
    #############################
    #### METRIC CALCULATIONS ####
    #############################
    
    # holds all relevant metric data and corresponding label
    metric_data = []
    
    # glycan RMSD calculation
    if dump_dir.endswith( '/' ):
        working_filename = "%stemp_working_%s.pdb" %( dump_dir, str( decoy_num ) )
        native_filename = "%stemp_native_%s.pdb" %( dump_dir, str( decoy_num ) )
    else:
        working_filename = "%s/temp_working_%s.pdb" %( dump_dir, str( decoy_num ) )
        native_filename = "%s/temp_native_%s.pdb" %( dump_dir, str( decoy_num ) )
    dump_pdb_by_chain( working_filename, working, working_Fc_glycan_chains, decoy_num, dump_dir = dump_dir )
    dump_pdb_by_chain( native_filename, native, native_Fc_glycan_chains, decoy_num, dump_dir = dump_dir )

    temp_working = Pose()
    try:
        temp_working.assign( load_pose( working_filename ) )
    except:
        pass

    temp_native = Pose()
    try:
        temp_native.assign( load_pose( native_filename ) )
    except:
        pass

    try:
        glycan_rmsd = non_peptide_heavy_atom_RMSD( temp_working, temp_native )
    except:
        glycan_rmsd = "nan"
        pass
    metric_data.append( "glycan_rmsd:" )
    metric_data.append( str( glycan_rmsd ) )
    
    # delete the files
    try:
        os.popen( "rm %s" %working_filename )
        os.popen( "rm %s" %native_filename )
    except:
        pass

    # pseduo-inferface energy 
    # ( full protein score - Fc-FcR glycan score [ except the short glycan away from interface ] )
    pseudo_interface_score = pseudo_interface_score_3ay4( working, sf, native = False, pmm = None ):    
    metric_data.append( "pseudo_interface_energy:" )
    metric_data.append( str( pseudo_interface_score ) )


    # delta interaction energy
    working_interaction_energy = calc_interaction_energy( working, sf, Vector1( [ JUMP_NUM ] ) )
    native_interaction_energy = calc_interaction_energy( native, sf, Vector1( [ JUMP_NUM ] ) )
    delta_interaction_energy = working_interaction_energy - native_interaction_energy
    metric_data.append( "delta_interface_interaction_energy:" )
    metric_data.append( str( delta_interaction_energy ) )


    # delta hbonds
    working_hbonds = get_hbonds( working )
    native_hbonds = get_hbonds( native )
    delta_hbonds = working_hbonds.nhbonds() - native_hbonds.nhbonds()
    metric_data.append( "delta_hbonds:" )
    metric_data.append( str( delta_hbonds ) )

    
    # delta interface residue contacts
    cutoff = 8
    working_intf_contacts, working_contact_list = count_interface_residue_contacts( JUMP_NUM, 
                                                                                    working, 
                                                                                    cutoff = cutoff )
    native_intf_contacts, native_contact_list = count_interface_residue_contacts( JUMP_NUM, 
                                                                                  native, 
                                                                                  cutoff = cutoff )
    delta_interface_res_contacts = working_intf_contacts - native_intf_contacts
    metric_data.append( "delta_interface_res_contacts_%s_A:" %( str( cutoff ) ) )
    metric_data.append( str( delta_interface_res_contacts ) )

    
    # delta interface sasa
    working_interface_sasa = calc_interface_sasa( working, JUMP_NUM )
    native_interface_sasa = calc_interface_sasa( native, JUMP_NUM )
    delta_interface_sasa = working_interface_sasa - native_interface_sasa
    metric_data.append( "delta_interface_sasa:" )
    metric_data.append( str( delta_interface_sasa ) )

    
    # MonteCarlo acceptance rate - if relevant
    if MC_acceptance_rate is not None:
        metric_data.append( "MonteCarlo_acceptance_rate:" )
        metric_data.append( str( MC_acceptance_rate ) )
    
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
                arr = re.split( "[:\s]+", line )
            
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
