#!/usr/bin/python

"""
Notes 9/15/16: find the glycans that are together and need to be sampled (that basically is working). Now need to figure out how to delete them in an approprite manner. Currently assuming that the glycans are found in order. Perhaps should check which glycan numbering bunch is lower (order each glycan from 1, 2, 3 etc. 1 deleted first, 2 deleted second and numbers adjusted by how many were delete in 1, 3 adjusted by num deleted in 1 and 2, etc). Keep relevant lines, delete HETNAM, ATOM, and LINK lines when needed
"""

from antibody_functions import initialize_rosetta, load_pose
from rosetta import get_fa_scorefxn, PyMOL_Mover
from rosetta.core.chemical.carbohydrates import CarbohydrateInfo
from rosetta.core.pose.carbohydrates import find_seqpos_of_saccharides_parent_residue
import sys
util_path = "project_utility_files/"
sys.path.append( util_path )
from line_definitions import PDB_line, HETNAM_line, SSBOND_line, LINK_line



initialize_rosetta()
pose_path = "/Users/mlnance/pyrosetta_dir/project_structs/fa_intra_rep_lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb"
native = load_pose( pose_path )
pmm = PyMOL_Mover()
pmm.keep_history( True )


native_Fc_glycan_A_nums = range( 216, 223 + 1 )
native_Fc_glycan_A_nums_except_core_GlcNAc = range( 217, 223 + 1 )
native_Fc_glycan_A_nums_except_core_GlcNAc_pdb = { 'D': range( 2, 5 + 1 ), 'E': range( 1, 3 + 1 ) }
native_Fc_glycan_B_nums = range( 440, 447 + 1 )
native_Fc_glycan_B_nums_except_core_GlcNAc = range( 441, 447 + 1 )
native_Fc_glycan_B_nums_except_core_GlcNAc_pdb = { 'F': range( 2, 5 + 1 ), 'G': range( 1, 3 + 1 ) }
A_to_B_glycan = { 'D': 'F', 'E': 'G' }
# since this is symmetrical, we should map the side1 to side2 numbers
side1_to_side2 = {}
for index_num in range( len( native_Fc_glycan_A_nums_except_core_GlcNAc ) ):
    side1_to_side2[ native_Fc_glycan_A_nums_except_core_GlcNAc[ index_num ] ] = native_Fc_glycan_B_nums_except_core_GlcNAc[ index_num ]


A_keys = native_Fc_glycan_A_nums_except_core_GlcNAc_pdb.keys()
for A_key in A_keys:
    B_key = A_to_B_glycan[ A_key ]


with open( pose_path, "rb" ) as fh:
    lines = fh.readlines()
pdb_lines = []
for line in lines:
    keep = True
    if line.startswith( "HETNAM" ):
        line = HETNAM_line( line )
    elif line.startswith( "SSBOND" ):
        line = SSBOND_line( line )
    elif line.startswith( "LINK" ):
        line = LINK_line( line )
    elif line.startswith( "ATOM" ) or line.startswith( "TER" ):
        line = PDB_line( line )
    else:
        keep = False

    if keep:
        pdb_lines.append( line )
    


keep_these_chains = [ 'A', 'B' ]
glycans_on_chain = {}
chem_edges = native.fold_tree().get_chemical_edges()
for edge in chem_edges:
    if native.residue( edge.start() ).is_protein() and native.pdb_info().chain( edge.start() ) in keep_these_chains:
        if native.residue( edge.stop() ).is_carbohydrate():
            glycan_residues = [ edge.stop() ]
            check_parent = edge.stop() + 1
            while True:
                try:
                    parent = find_seqpos_of_saccharides_parent_residue( native.residue( check_parent ) )
                except:
                    break
                if parent in glycan_residues:
                    glycan_residues.append( check_parent )
                    check_parent += 1
                else:
                    break
            # remove the core residue ( the edge.stop() connected to edge.start() )
            glycan_residues.remove( edge.stop() )
            glycans_on_chain[ native.pdb_info().chain( edge.start() ) ] = glycan_residues


# dict{ chain: lines keep } then do the set of the lines to keep and keep the ones that are common between both chains
keys = glycans_on_chain.keys()
keys.sort()
keep_lines = {}
second_chain = False
#for chain in keys:
'''
This doesn't work because it's not keeping the residues upper from the test residue. Have to address that somehow
'''
for chain in [ 'A' ]:
    # glycan on this chain
    glycan_residues = glycans_on_chain[ chain ]

    for index in range( len( glycan_residues ) ):
        # residues in this glycan
        test_res = glycan_residues[ index ]
        # uniq name = resnum_reschain
        test_res_uniq_name = '_'.join( native.pdb_info().pose2pdb( test_res ).split() )
        del_res = glycan_residues[ index + 1 : ]
        nres_del = len( del_res )
        keep_these_lines = []

        for line in pdb_lines:
            keep_line = True
            if line.line.startswith( "HETNAM" ):
                if '_'.join( [ str( line.res_num ), line.res_chain ] ) != test_res_uniq_name:
                    keep_line = False
            elif line.line.startswith( "LINK" ):
                if '_'.join( [ str( line.res1_num ), line.res1_chain ] ) != test_res_uniq_name:
                    if '_'.join( [ str( line.res2_num ), line.res2_chain ] ) != test_res_uniq_name:
                        # this is in a LINK record, so check if the residue is the upper or lower branch (keep if res is a BRANCH_LOWER_TERMINUS_VARIANT)
                        if not "BRANCH_LOWER_TERMINUS_VARIANT" in native.residue_type( test_res ).variant_types():
                            keep_line = False
            elif line.line.startswith( "ATOM" ):
                if '_'.join( [ str( line.res_num ), line.res_chain ] ) != test_res_uniq_name:
                    keep_line = False
            else:
                # TER line
                pass

            if keep_line:
                keep_these_lines.append( line.line )

        keep_lines[ nres_del ] = keep_these_lines
