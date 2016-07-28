#!/usr/bin/python
__author__ = "morganlnance"



# TODO: need to figure out the ordering of the FoldTree or something because the jumps get added successfully, but then crashes when you try to access them. kinematics::Atom bad method call to Atom tree
def pseudo_glycan_interface_energy( sugar_nums, in_sf, in_pose, pmm = None ):
    """
    Gets chemical edges from <in_pose> and turns each one into a jump if it contains a protein residue on one side and a sugar residue found in <sugar_nums> on the other.
    <sugar_nums> contains the glycan(s) that are connected to <in_pose> by chemical edge(s).
    If all the sugar numbers of interest in <sugar_nums> are not connected to each other and then the protein by a chemical edge, they will not be moved from the complex and the score will be inaccurate
    :param sugar_nums: list( Pose numbers for sugar residues of interest )
    :param in_sf: ScoreFunction
    :param in_pose: Pose
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: float( pseudo interface energy )
    """
    import sys
    from rosetta.core.scoring import score_type_from_name


    # check that all residues in <sugar_nums> are actually sugars
    for resnum in sugar_nums:
        if not in_pose.residue( resnum ).is_carbohydrate():
            print "\nYou gave me a non-sugar residue number in your <sugar_nums> argument. I don't like the looks of this. Exiting."
            sys.exit()

    # clone the input pose
    pose = in_pose.clone()

    # set atom_pair_constraint weight to 0 in the sf clone
    sf = in_sf.clone()
    sf.set_weight( score_type_from_name( "atom_pair_constraint" ), 0.0 )

    # get the score of the whole complex
    start_score = sf( pose )

    # get the chemical edges found in the Pose
    chem_edges = pose.fold_tree().get_chemical_edges()

    # collect the JUMP nums that are added
    jump_nums_added = []

    # for each edge, see which one has a protein residue and a sugar residue found in <sugar_nums>
    for chem_edge in chem_edges:
        # turn the chemical edge into a string
        chem_edge_str = str( chem_edge ).strip()

        # split the string on whitespace as to access residue numbers
        chem_edge_split_with_empty_string = chem_edge_str.split( ' ' )

        # remove the empty string occurances ( '' ) as they may be of variable occurance
        chem_edge_split = []
        for char in chem_edge_split_with_empty_string:
            if char != '':
                chem_edge_split.append( char )

        # collect the relevant information from the split string
        upstream_connect = int( chem_edge_split[1] )
        downstream_connect = int( chem_edge_split[2] )

        ## check for the two possible cases of interest
        # if the upstream is a protein residue and the downstream is a sugar in <sugar_nums>
        if pose.residue( upstream_connect ).is_protein() and downstream_connect in sugar_nums:
            # remove the chemical edge using the chem_edge object
            pose.fold_tree().delete_edge( chem_edge )

            # add a new jump edge using the information from the chemical edge
            add_jump_num = pose.fold_tree().num_jump() + 1
            jump_nums_added.append( add_jump_num )
            pose.fold_tree().add_edge( upstream_connect, downstream_connect, add_jump_num )

            # reorder the new FoldTree on the first residue
            pose.fold_tree().reorder( 1 )
        
        # elif the upstream is a sugar in <sugar_nums> and the downstream is a protein residue
        elif pose.residue( downstream_connect ).is_protein() and upstream_connect in sugar_nums:
            # remove the chemical edge using the chem_edge object
            pose.fold_tree().delete_edge( chem_edge )

            # add a new jump edge using the information from the chemical edge
            add_jump_num = pose.fold_tree().num_jump() + 1
            jump_nums_added.append( add_jump_num )
            pose.fold_tree().add_edge( upstream_connect, downstream_connect, add_jump_num )

            # reorder the new FoldTree on the first residue
            pose.fold_tree().reorder( 1 )

        # else, skip this chemical edge
        else:
            pass

    # send the Pose to PyMOL
    try:
        pmm.keep_history( True )
        pmm.apply( pose )
    except:
        pass

    # split apart the <sugar_nums>
    for jump_num in jump_nums_added:
        try:
            jump = pose.jump( jump_num )
        except:
            return pose
        vec = jump.get_translation() * 1000
        jump.set_translation( vec )
        pose.set_jump( jump_num, jump )
        try:
            pmm.apply( pose )
        except:
            pass

    return pose



def pseudo_interface_energy_3ay4( pose, in_sf, native = False, pmm = None ):
    """
    Attempts to get pseudo-interface energy of a glycosylated 3ay4 decoy
    Lots of hard coding here - works on a decoy pose as Rosetta renumbers the Pose a bit
    Makes the two ASN connections to the Fc A and B glycans JUMPs instead of chemical EDGEs
    :param pose: Pose
    :param in_sf: ScoreFunction
    :param native: bool( is this the native 3ay4 or a decoy? Answer determines how FoldTree gets coded )
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: float( pseudo interface energy )
    """
    from rosetta import FoldTree, Pose
    from rosetta.numeric import xyzVector_Real
    from rosetta.core.scoring import score_type_from_name


    # set atom_pair_constraint weight to 0
    sf = in_sf.clone()
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



def Fc_glycan_rmsd( working, native, working_Fc_glycan_chains, native_Fc_glycan_chains, decoy_num, dump_dir, return_hbonds_too = False ):
    '''
    Return the glycan RMSD contribution of the two Fc glycans in 3ay4 (may work for other PDBs, but I don't know yet)
    :param working: decoy Pose()
    :param native: native Pose()
    :param working_Fc_glycan_chains: list( the chain id's for the working Fc glycan ). Ex = [ 'H', 'I' ]
    :param native_Fc_glycan_chains: list( the chain id's for the native Fc glycan ). Ex = [ 'D', 'E' ]
    :param decoy_num: int( the number of the decoy for use when dumping its Fc glycan )
    :param dump_dir: str( /path/to/dump_dir for the temp pdb files made. Files will be deleted )
    :param return_hbonds_too: bool( do you also want to return the number of hbonds in the working Fc glycan? ) Default = False
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
        working_filename = "%s%s_temp_working_glyc_%s.pdb" %( dump_dir, id, str( decoy_num ) )
        native_filename = "%s%s_temp_native_glyc_%s.pdb" %( dump_dir, id, str( decoy_num ) )
    else:
        working_filename = "%s/%s_temp_working_glyc_%s.pdb" %( dump_dir, id, str( decoy_num ) )
        native_filename = "%s/%s_temp_native_glyc_%s.pdb" %( dump_dir, id, str( decoy_num ) )

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

    # get the hbonds of the working Fc glycans, if desired
    if return_hbonds_too:
        from rosetta import get_fa_scorefxn
        from toolbox import get_hbonds

        # score the Pose so that its hbond energies gets updated
        sf = get_fa_scorefxn()
        sf( temp_working )

        working_Fc_glycan_hbonds = get_hbonds( temp_working ).nhbonds()
        return glycan_rmsd, working_Fc_glycan_hbonds

    # or just return the glycan rmsd
    else:
        return glycan_rmsd



def Fc_glycan_hbonds( working, working_Fc_glycan_chains, decoy_num, dump_dir ):
    '''
    Return the number of hbonds the Fc glycan contributes to the Pose. Total hbonds in Pose - Total hbonds in Pose without Fc glycans = hbonds contributed by Fc glycans
    :param working: decoy Pose()
    :param native: native Pose()
    :param working_Fc_glycan_chains: list( the chain id's for the working Fc glycan ). Ex = [ 'H', 'I' ]
    :param decoy_num: int( the number of the decoy for use when dumping its Fc glycan )
    :param dump_dir: str( /path/to/dump_dir for the temp pdb files made. Files will be deleted )
    :return: float( Fc glycan rmsd )
    '''
    #################
    #### IMPORTS ####
    #################

    # Rosetta functions
    from rosetta import Pose, get_fa_scorefxn
    from toolbox import get_hbonds

    # Rosetta functions I wrote out
    from antibody_functions import load_pose

    # utility functions
    import os
    from util import dump_pdb_by_chain, id_generator



    # get temporary files to work with
    id = id_generator()
    if dump_dir.endswith( '/' ):
        working_filename = "%s%s_hbtemp_working_no_glyc_%s.pdb" %( dump_dir, id, str( decoy_num ) )
    else:
        working_filename = "%s/%s_hbtemp_working_no_glyc_%s.pdb" %( dump_dir, id, str( decoy_num ) )

    # get the chain id's of everything discluding the passed Fc glycan chain id's
    working_pose_chains = []
    for res in working:
        chain_id = working.pdb_info().chain( res.seqpos() )
        if ( chain_id not in working_pose_chains ) and ( chain_id not in working_Fc_glycan_chains ):
            working_pose_chains.append( chain_id )

    # dump out the pose without its Fc glycans by the chain id's
    dump_pdb_by_chain( working_filename, working, working_pose_chains, decoy_num, dump_dir = dump_dir )

    # load in the working Pose without the Fc glycans
    temp_working = Pose()
    try:
        temp_working.assign( load_pose( working_filename ) )
    except:
        pass

    # score the Poses so their hbond energies get updated
    sf = get_fa_scorefxn()
    sf( working )
    sf( temp_working )

    # get the number of hbonds in the Pose without the Fc glycans
    with_Fc_glycan_hbonds = get_hbonds( working )
    no_Fc_glycan_hbonds = get_hbonds( temp_working )
    Fc_glycan_hbonds = with_Fc_glycan_hbonds.nhbonds() - no_Fc_glycan_hbonds.nhbonds()

    # delete the files
    try:
        os.popen( "rm %s" %working_filename )
    except:
        pass

    return Fc_glycan_hbonds



def check_GlcNAc_to_Phe_contacts( working, cutoff, native = False ):
    """
    Returns 0 if GlcNAc-2 (on Gal) contacts Phe-243, 1 if at least one does, and 2 if they both do within the given <cutoff>
    :param working: Pose
    :param cutoff: float( the distance cutoff to count this as a contact )
    :param native: bool( is this the native Pose? Or a glycosylated decoy of 3ay4? ) Default = False
    return int( 0, 1, or 2 )
    """
    from antibody_functions import calc_distance

    # get the residue objects by using pre-determined Pose numbers
    if not native:
        GlcNAc_side_A = working.residue( 609 )
        Phe_side_A = working.residue( 15 )
        GlcNAc_side_B = working.residue( 617 )
        Phe_side_B = working.residue( 230 )
    else:
        GlcNAc_side_A = working.residue( 222 )
        Phe_side_A = working.residue( 15 )
        GlcNAc_side_B = working.residue( 446 )
        Phe_side_B = working.residue( 238 )

    # get the distances in question
    dist_A = calc_distance( list( GlcNAc_side_A.nbr_atom_xyz() ), list( Phe_side_A.nbr_atom_xyz() ) )
    dist_B = calc_distance( list( GlcNAc_side_B.nbr_atom_xyz() ), list( Phe_side_B.nbr_atom_xyz() ) )

    # check to see if they are within the cutoff distance
    num_Phe_GlcNAc_contacts = 0
    if dist_A <= cutoff:
        num_Phe_GlcNAc_contacts += 1
    if dist_B <= cutoff:
        num_Phe_GlcNAc_contacts += 1

    return num_Phe_GlcNAc_contacts
