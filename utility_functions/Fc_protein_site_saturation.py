#!/usr/bin/python
__author__="morgalnance"


# import Rosetta functions
from antibody_functions import initialize_rosetta, load_pose, \
    make_pack_rotamers_mover, get_res_nums_within_radius, \
    native_Fc_glycan_nums, native_Fc_protein_nums, \
    get_contact_map_between_range1_range2
from rosetta import Pose, get_fa_scorefxn, PyMOL_Mover, \
    MoveMap, MinMover
from rosetta.core.chemical import VariantType
from rosetta.core.pose import remove_variant_type_from_pose_residue
from rosetta.core.pose import add_variant_type_to_pose_residue

# mutation-related imports
from rosetta import pose_from_sequence, ResidueFactory
#from rosetta.protocols.simple_moves import MutateResidue
#from toolbox import mutate_residue

# general imports
import os, sys
import pandas as pd



# create global pymol object for viewing
pmm = PyMOL_Mover()
pmm.keep_history()

# path for mutational analysis data
data_dir = os.getcwd() + "/mutational_data/"
if not os.path.isdir( data_dir ):
    os.mkdir( data_dir )
data_filename = "Fc_protein_site_saturation_fresh_start_low_E_native.csv"

# initialize rosetta with sugar flags
initialize_rosetta()

# create pose object, assign info object, and send pose to pymol
native_pose = load_pose( "/Users/Research/pyrosetta_dir/structures_from_jazz/3ay4_Fc_FcgRIIIa/fresh_start/project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb" )
native_pose.pdb_info().name( "native_pose" )
info = native_pose.pdb_info()
#pmm.apply( native_pose )

# instantiate full atom score function and score pose
sf = get_fa_scorefxn()
nat_E = sf( native_pose )


# create some constant data
#AA_options = [ "ALA", "CYS", "ASP", "GLN", "PHE", "GLY", "HIS", "ILE", "LYS", "MET", "ASN", "PRO", "GLU", "ARG", "SER", "THR", "THR", "VAL", "TRP", "TYR" ]
#AA_name1_list = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
AA_options = [ "ALA", "CYS", "ASP", "GLN", "PHE", "GLY", "HIS", "ILE", "LYS", "MET", "ASN", "GLU", "ARG", "SER", "THR", "THR", "VAL", "TRP", "TYR" ]
AA_name1_list = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
PACK_RADIUS = 5.0

# make a Pose out of the residues you want to use as mutants
site_saturation_pose = pose_from_sequence( ''.join( AA_name1_list ) )

# remove the LOWER_TERMINUS_VARIANT and UPPER_TERMINUS_VARIANT VariantTypes from
# the first and last residue of the site_saturation_pose
remove_variant_type_from_pose_residue( site_saturation_pose, VariantType.LOWER_TERMINUS_VARIANT, 1 )
remove_variant_type_from_pose_residue( site_saturation_pose, VariantType.UPPER_TERMINUS_VARIANT, site_saturation_pose.n_residue() )


# insantiate lists for data (will be put into pandas dataframe)
orig_AA_list = []
new_AA_list = []
pose_position = []
pdb_position = []
pdb_chain = []
within_10A_of_Fc_glycan_res = []
native_E = []
native_E_res = []
mut_E = []
mut_E_res = []
ddG = []
ddG_res = []

# used to build a new Ala residue
res_factory = ResidueFactory()

# get the 10A contact map between Fc protein and Fc glycan for use in determining if the Fc protein mutation
# is within 10A of an Fc glycan residue
Fc_protein_to_Fc_glycan_cmap = get_contact_map_between_range1_range2( native_Fc_protein_nums, native_Fc_glycan_nums, native_pose, cutoff = 10 )


# mutate all residues to Ala
# make sure residue is not a sugar or a branch point
# also skipping Cys residues because they're just being sassy
residues_that_didnt_work = []
#for seq_pos in range(1, native_pose.total_residue() + 1):
#for seq_pos in native_Fc_protein_nums:
for seq_pos in Fc_protein_to_Fc_glycan_cmap.keys():
    res = native_pose.residue( seq_pos )

    # replacing alanines too for the heck of it. The ddG should be very close to 0 if this works
    if not res.name1() == 'C':
        if not res.is_carbohydrate():
            if not res.is_branch_point():
                # build the mutation residue (include the to-normal mutation too)
                for new_AA in AA_name1_list:
                    # store all data for this residue's mutation for the df
                    # append original residue, position, and score to list
                    orig_AA_name1 = str( native_pose.residue( seq_pos ).name1() )
                    orig_AA_list.append( orig_AA_name1 )
                    pose_position.append( seq_pos )
                    pdb_position.append( native_pose.pdb_info().pose2pdb( seq_pos ).strip().split( ' ' )[0] )
                    pdb_chain.append( native_pose.pdb_info().pose2pdb( seq_pos ).strip().split( ' ' )[1] )
                    native_E.append( nat_E )
                    new_AA_list.append( new_AA )

                    # if this is an Fc protein residue, determine if it is within 10A of the Fc glycan
                    # use the contact map created earlier to determine if this residue is near the glycan
                    # since the contact map was made between the protein and glycan, just check to see if this res is a key
                    # keeping this here for now in case I decide to do the alanine scan on the whole protein
                    if seq_pos in native_Fc_protein_nums:                        
                        within_10A_of_Fc_glycan_res.append( seq_pos in Fc_protein_to_Fc_glycan_cmap.keys() )
                        
                    # get a copy of the Pose
                    mutant = Pose()
                    mutant.assign( native_pose )

                    # get the res type from the site_saturation_pose using the list index+1 (since Pose numbering starts at 1)
                    new_AA_restype = site_saturation_pose.conformation().residue_type( AA_name1_list.index( new_AA ) + 1 )

                    # build the mutation residue
                    # this function uses the backbone information from the current residue at this position to build the new Ala
                    # if this isn't a Gly residue, use the CB information as well
                    if mutant.residue( seq_pos ).name1() != "G":
                        new_AA_res = res_factory.create_residue( new_AA_restype, mutant.residue( seq_pos ), mutant.conformation(), preserve_c_beta=True )
                    # otherwise, don't preserve CB because Gly doesn't have a CB
                    else:
                        new_AA_res = res_factory.create_residue( new_AA_restype, mutant.residue( seq_pos ), mutant.conformation() )


                    # mutate!
                    try:
                        mutant.replace_residue( seq_pos, new_AA_res, orient_backbone=True )

                        # change pose name to (orig AA)(pdb number)_to_newAA
                        mutant.pdb_info().name( "%s%s_to_%s" %( orig_AA_name1, str( mutant.pdb_info().pose2pdb( seq_pos ).strip().replace( ' ', '' ) ), new_AA ))


                        # get residue numbers (including mutation site) to be packed and minimized
                        res_nums_around_mutation_site = get_res_nums_within_radius( seq_pos, mutant, PACK_RADIUS, include_seq_pos = True )


                        # pack around mutation
                        pack_rotamers_mover = make_pack_rotamers_mover( sf, mutant,
                                                                        apply_sf_sugar_constraints = False,
                                                                        pack_branch_points = True,
                                                                        residue_range = res_nums_around_mutation_site )
                        pack_rotamers_mover.apply( mutant )

                        # minimize around mutation
                        min_mm = MoveMap()
                        for res_num in res_nums_around_mutation_site:
                            min_mm.set_bb( res_num, True )
                            min_mm.set_chi( res_num, True )
                        min_mover = MinMover( movemap_in = min_mm,
                                              scorefxn_in = sf,
                                              min_type_in = "dfpmin_strong_wolfe",
                                              tolerance_in = 0.01,
                                              use_nb_list_in = True )
                        min_mover.apply(mutant)

                        # visualize mutation
                        #pmm.apply(mutant)


                        # score and add to list for dataframe
                        new_E = sf( mutant )
                        new_E_res = mutant.energies().residue_total_energy( seq_pos )
                        nat_E_res = native_pose.energies().residue_total_energy( seq_pos )

                        mut_E.append( new_E )
                        mut_E_res.append( new_E_res )
                        native_E_res.append( nat_E_res )
                        ddG.append( new_E - nat_E )
                        ddG_res.append( new_E_res - nat_E_res )
                    except:
                        residues_that_didnt_work.append( seq_pos )
                        pass


# print data and add data to dataframe
df = pd.DataFrame()
df["orig_AA"] = orig_AA_list
df["new_AA"] = new_AA_list
df["pose_num"] = pose_position
df["pdb_num"] = pdb_position
df["pdb_chain"] = pdb_chain
df["within_10A_of_Fc_glycan_res"] = within_10A_of_Fc_glycan_res
df["native_E"] = native_E
df["mut_E"] = mut_E
df["native_E_res"] = native_E_res
df["mut_E_res"] = mut_E_res
df["ddG"] = ddG
df["ddG_res"] = ddG_res
print df

# output results to file
with open(data_dir + data_filename, 'a') as f:
    df.to_csv(data_dir + data_filename)
print 'Data written to:', data_dir + data_filename
