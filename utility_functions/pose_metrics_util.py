#!/usr/bin/python
__author__ = "morganlnance"



def pseudo_interface_energy_3ay4( pose, sf, native = False, pmm = None ):
    """
    Attempts to get pseudo-interface energy of a glycosylated 3ay4 decoy
    Lots of hard coding here - works on a decoy pose as Rosetta renumbers the Pose a bit
    Makes the two ASN connections to the Fc A and B glycans JUMPs instead of chemical EDGEs
    :param pose: Pose
    :param sf: ScoreFunction
    :param native: bool( is this the native 3ay4 or a decoy? Answer determines how FoldTree gets coded )
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: float( pseudo interface energy )
    """
    from rosetta import FoldTree, Pose
    from rosetta.numeric import xyzVector_Real
    from rosetta.core.scoring import score_type_from_name


    # set atom_pair_constraint weight to 0
    sf.set_weight( score_type_from_name( "atom_pair_constraint" ), 0.0 )

    # get the score of the whole complex
    start_score = sf( pose )

    # hard code the new FoldTree specific to a glycosylated decoy of 3ay4
    if not native:
        ft = FoldTree()
        ft.add_edge( 1, 215, -1 )
        ft.add_edge( 1, 216, 1 )    # beginning of chain A to beginning of chain B
        ft.add_edge( 216, 431, -1 )
        ft.add_edge( 1, 432, 2 )    # beginning of chain A to beginning of chain C
        ft.add_edge( 432, 591, -1 )
        ft.add_edge( 579, 592, "ND2", "C1" )
        ft.add_edge( 592, 596, -1 )
        ft.add_edge( 594, 597, "O6", "C1" )
        ft.add_edge( 597, 598, -1 )
        ft.add_edge( 592, 599, "O6", "C1" )
        ft.add_edge( 462, 600, "ND2", "C1" )
        ft.add_edge( 600, 602, -1 )
        ft.add_edge( 69, 603, 3 )   # ASN 297 A to core GlcNAc H
        ft.add_edge( 603, 607, -1 )
        ft.add_edge( 605, 608, "O6", "C1" )
        ft.add_edge( 608, 610, -1 )
        ft.add_edge( 284, 611, 4 )  # ASN 297 B to core GlcNAc J
        ft.add_edge( 611, 615, -1 )
        ft.add_edge( 613, 616, "O6", "C1" )
        ft.add_edge( 616, 618, -1 )

    # hard code the new FoldTree specific to the native 3ay4
    else:
        ft = FoldTree()
        ft.add_edge( 1, 215, -1 )
        ft.add_edge( 69, 216, 1 )   # ASN 297 A to core GlcNAc D
        ft.add_edge( 216, 220, -1 )
        ft.add_edge( 218, 221, "O6", "C1" )
        ft.add_edge( 221, 223, -1 )
        ft.add_edge( 1, 224, 2 )    # beginning of chain A to beginning of chain B
        ft.add_edge( 224, 439, -1 )
        ft.add_edge( 292, 440, 3 )  # ASN 297 B to core GlcNAc E
        ft.add_edge( 440, 444, -1 )
        ft.add_edge( 442, 445, "O6", "C1" )
        ft.add_edge( 445, 447, -1 )
        ft.add_edge( 1, 448, 4 )    # beginning of chain A to beginning of chain C
        ft.add_edge( 448, 607, -1 )
        ft.add_edge( 595, 608, "ND2", "C1" )
        ft.add_edge( 608, 612, -1 )
        ft.add_edge( 610, 613, "O6", "C1" )
        ft.add_edge( 613, 614, -1 )
        ft.add_edge( 608, 615, "O6", "C1" )
        ft.add_edge( 478, 616, "ND2", "C1" )
        ft.add_edge( 616, 618, -1 )
        
    # make a temporary Pose and give the new FoldTree to it
    try:
        pmm.keep_history( True )
    except:
        pass
    temp_pose = Pose()
    temp_pose.assign( pose )
    temp_pose.fold_tree( ft )
    try:
        pmm.apply( temp_pose )
    except:
        pass

    # split apart the two Fc sugars one-by-one
    # if decoy structure -- the two glycans are now the last two new jumps
    if not native:
        jump = temp_pose.jump( 3 ) # sugar A
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 3, jump )
        try:
            pmm.apply( temp_pose )
        except:
            pass
        jump = temp_pose.jump( 4 ) # sugar B
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 4, jump )
        try:
            pmm.apply( temp_pose )
        except:
            pass

    # else native structure -- the two glycans are the first and third new jumps
    else:
        jump = temp_pose.jump( 1 ) # sugar A
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 1, jump )
        try:
            pmm.apply( temp_pose )
        except:
            pass
        jump = temp_pose.jump( 3 ) # sugar B
        vec = xyzVector_Real( 1000, 1000, 1000 )
        jump.set_translation( vec )
        temp_pose.set_jump( 3, jump )
        try:
            pmm.apply( temp_pose )
        except:
            pass        

    # score the split-apart Pose
    split_score = sf( temp_pose )

    # get the pseudo-interface score
    # total - split = interface ( ie. interface + split = total )
    pseudo_interface_energy = start_score - split_score

    return pseudo_interface_energy



def Fc_glycan_rmsd( working, native, working_Fc_glycan_chains, native_Fc_glycan_chains, decoy_num, dump_dir ):
    '''
    Return the glycan RMSD of the two Fc glycans in 3ay4 (may work for other PDBs, but I don't know yet)
    :param working: decoy Pose()
    :param native: native Pose()
    :param working_Fc_glycan_chains: list( the chain id's for the working Fc glycan ). Ex = [ 'H', 'I' ]
    :param native_Fc_glycan_chains: list( the chain id's for the native Fc glycan ). Ex = [ 'D', 'E' ]
    :param decoy_num: int( the number of the decoy for use when dumping its Fc glycan )
    :param dump_dir: str( /path/to/dump_dir for the temp pdb files made. Files will be deleted )
    :return: float( Fc glycan rmsd )
    '''
    #################
    #### IMPORTS ####
    #################

    # Rosetta functions
    from rosetta import Pose
    from rosetta.core.scoring import non_peptide_heavy_atom_RMSD

    # Rosetta functions I wrote out
    from antibody_functions import load_pose

    # utility functions
    import os
    from util import dump_pdb_by_chain, id_generator



    # get temporary files to work with
    id = id_generator()
    if dump_dir.endswith( '/' ):
        working_filename = "%s%s_temp_working_%s.pdb" %( dump_dir, id, str( decoy_num ) )
        native_filename = "%s%s_temp_native_%s.pdb" %( dump_dir, id, str( decoy_num ) )
    else:
        working_filename = "%s/%s_temp_working_%s.pdb" %( dump_dir, id, str( decoy_num ) )
        native_filename = "%s/%s_temp_native_%s.pdb" %( dump_dir, id, str( decoy_num ) )

    # dump out the Fc glycans by their chain id's
    dump_pdb_by_chain( working_filename, working, working_Fc_glycan_chains, decoy_num, dump_dir = dump_dir )
    dump_pdb_by_chain( native_filename, native, native_Fc_glycan_chains, decoy_num, dump_dir = dump_dir )

    # load in the Fc glycans
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

    # calculate the glycan rmsd
    try:
        glycan_rmsd = non_peptide_heavy_atom_RMSD( temp_working, temp_native )
    except:
        glycan_rmsd = "nan"
        pass
    
    # delete the files
    try:
        os.popen( "rm %s" %working_filename )
        os.popen( "rm %s" %native_filename )
    except:
        pass

    return glycan_rmsd
