#!/usr/bin/python
__author__="morgalnance"


# import Rosetta functions
from antibody_functions import initialize_rosetta, load_pose, \
    native_Fc_glycan_nums, native_Fc_protein_nums, \
    native_FcR_protein_nums, get_contact_map_between_range1_range2, \
    AA_name1_list, mutate_residue, get_interface_score
from rosetta import Pose, PyMOL_Mover, get_fa_scorefxn
from toolbox import get_hbonds

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
data_filename = "Fc_protein_near_Fc_glycan_10A_site_saturation_fresh_start_low_E_native.csv"

# initialize rosetta with sugar flags
initialize_rosetta()

# create pose object, assign info object, and send pose to pymol
native_pose = load_pose( "/Users/mlnance/pyrosetta_dir/project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb" )
native_pose.pdb_info().name( "native_pose" )
info = native_pose.pdb_info()
#pmm.apply( native_pose )

# instantiate full atom score function and score native pose
sf = get_fa_scorefxn()
native_E = sf( native_pose )
native_E_interface = get_interface_score( 2, sf, native_pose )
native_nhbonds = get_hbonds( native_pose ).nhbonds()


# insantiate lists for data (will be put into pandas dataframe)
orig_AA_list = []
pose_position = []
pdb_position = []
pdb_chain = []
new_AA_list = []
was_best_mutation = []
within_10A_of_FcR_interface = []
native_E_list = []
mut_E_list = []
ddG_list = []
native_E_res_list = []
mut_E_res_list = []
ddG_res_list = []
native_E_interface_list = []
mut_E_interface_list = []
ddG_interface_list = []
dhbonds_list = []


# get the 10A contact map between Fc protein and Fc glycan
# Fc protein residues within 10A of the Fc glycan will be the ones getting mutated
# my code will skip over branch points so ASN 297 is okay to have in this contact map
Fc_protein_to_Fc_glycan_cmap = get_contact_map_between_range1_range2( native_Fc_protein_nums, 
                                                                      native_Fc_glycan_nums, 
                                                                      native_pose, 
                                                                      cutoff = 10 )

# get the 10A contact map between Fc protein residues that are within 10A of the Fc glycan to the FcR protein
Fc_protein_to_FcR_protein_cmap = get_contact_map_between_range1_range2( Fc_protein_to_Fc_glycan_cmap.keys(), 
                                                                        native_FcR_protein_nums, 
                                                                        native_pose, 
                                                                        cutoff = 10 )


# for each Fc protein residue near the Fc glycan, mutate each residue to all 20 amino acids
# and pack within 10A around the mutation site
for seq_pos in Fc_protein_to_Fc_glycan_cmap.keys():
    res = native_pose.residue( seq_pos )

    # used to store scores for each mutation to then determine which was best
    AA_to_score_dict = {}

    # skipping Cys residues because they're just being sassy
    if not res.name1() == 'C':
        # don't do sugars
        if not res.is_carbohydrate():
            # skip branch points like ASN 297 since it connects to GlcNAc 1
            if not res.is_branch_point():
                # for each amino acid
                for new_AA in AA_name1_list:
                    # get a copy of the Pose
                    mutant = Pose()
                    mutant.assign( native_pose )

                    # mutate!
                    mutant.assign( mutate_residue( seq_pos, new_AA, mutant, sf, pack_radius = 10 ) )

                    # change pose name to (orig AA)(pdb number and chain)_to_newAA
                    orig_AA_name1 = str( native_pose.residue( seq_pos ).name1() )
                    mutant.pdb_info().name( "%s%s_to_%s" %( orig_AA_name1, str( mutant.pdb_info().pose2pdb( seq_pos ).strip().replace( ' ', '' ) ), new_AA ))

                    # visualize mutation
                    #pmm.apply(mutant)

                    # store all information
                    # store all data for this residue's mutation for the df
                    # append original residue, position, and score to list
                    orig_AA_list.append( orig_AA_name1 )
                    pose_position.append( seq_pos )
                    pdb_position.append( native_pose.pdb_info().pose2pdb( seq_pos ).strip().split( ' ' )[0] )
                    pdb_chain.append( native_pose.pdb_info().pose2pdb( seq_pos ).strip().split( ' ' )[1] )
                    new_AA_list.append( new_AA )
                    native_E_list.append( native_E )
                    native_E_interface_list.append( native_E_interface )

                    # check if this Fc protein residue is within 10A of the FcR protein
                    # only Fc protein residues within 10A of the Fc glycan are used in this count, so this is
                    # checking if they're also within 10A of the FcR protein
                    if seq_pos in Fc_protein_to_FcR_protein_cmap.keys():
                        within_10A_of_FcR_interface.append( True )
                    else:
                        within_10A_of_FcR_interface.append( False )

                    # score and store the native residue
                    native_E_res = native_pose.energies().residue_total_energy( seq_pos )
                    native_E_res_list.append( native_E_res )

                    # score and store the mutation
                    # total_score
                    mut_E = sf( mutant )
                    mut_E_list.append( mut_E )
                    ddG_list.append( mut_E - native_E )
                    # residue score
                    mut_E_res = mutant.energies().residue_total_energy( seq_pos )
                    mut_E_res_list.append( mut_E_res )
                    ddG_res_list.append( mut_E_res - native_E_res )
                    # interface score
                    mut_E_interface = get_interface_score( 2, sf, mutant )
                    mut_E_interface_list.append( mut_E_interface )
                    ddG_interface_list.append( mut_E_interface - native_E_interface )
                    # collect hbond information
                    mut_nhbonds = get_hbonds( native_pose ).nhbonds()
                    dhbonds_list.append( mut_nhbonds - native_nhbonds )

                    # collect information for this particular information to score to determine which mutation was best
                    # this will be based on the ddG of the total_score
                    AA_to_score_dict[ new_AA ] = mut_E - native_E

                # determine which mutation was the best for this particular residue
                lowest_score = None
                best_mutation = None
                for amino_acid in AA_to_score_dict.keys():
                    # if this is the first mutation checked, update lowest_score and amino_acid with this mutation's information
                    if lowest_score is None:
                        lowest_score = AA_to_score_dict[ amino_acid ]
                        best_mutation = amino_acid
                    # otherwise it is not empty, so compare the other mutations to this one. Keep the lowest score
                    else:
                        if AA_to_score_dict[ amino_acid ] < lowest_score:
                            lowest_score = AA_to_score_dict[ amino_acid ]
                            best_mutation = amino_acid

                # store the final information as to which amino acid mutation was best for this particular residue
                # lists go in order, so itering through AA_name1_list like this will result in the same order as they were created
                for amino_acid in AA_name1_list:
                    if amino_acid == best_mutation:
                        was_best_mutation.append( True )
                    else:
                        was_best_mutation.append( False )



# print data and add data to dataframe
df = pd.DataFrame()
df["orig_AA"] = orig_AA_list
df["new_AA"] = new_AA_list
df["pose_num"] = pose_position
df["pdb_num"] = pdb_position
df["pdb_chain"] = pdb_chain
df["within_10A_of_FcR_interface"] = within_10A_of_FcR_interface
df["was_best_mutation"] = was_best_mutation
df["native_E"] = native_E_list
df["mut_E"] = mut_E_list
df["ddG"] = ddG_list
df["native_E_res"] = native_E_res_list
df["mut_E_res"] = mut_E_res_list
df["ddG_res"] = ddG_res_list
df["native_E_interface"] = native_E_interface_list
df["mut_E_interface"] = mut_E_interface_list
df["ddG_interface"] = ddG_interface_list
df["dhbonds"] = dhbonds_list
print df

# output results to file
with open(data_dir + data_filename, 'a') as f:
    df.to_csv(data_dir + data_filename)
print 'Data written to:', data_dir + data_filename
