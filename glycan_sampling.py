#!/usr/bin/python

'''
NOTES
STARTING POSE (3ay4 without Fc) is a base structure that was acquired from making 1000 pack/min decoys of native 3ay4 (with Fc sugar on) and then manually removing the Fc sugars of the lowest energy decoy
'''

'''
SAMPLE INPUT
run glycan_sampling.py pdb_copies_dont_touch/lowest_E_single_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb pdb_copies_dont_touch/lowest_E_single_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII_removed_Fc_sugar.pdb database/chemical/carbohydrates/common_glycans/3ay4_Fc_Glycan.iupac
'''



#########################
#### PARSE ARGUMENTS ####
#########################

import argparse

# parse and store input arguments
parser = argparse.ArgumentParser(description="Use PyRosetta to glycosylate a pose and find a low E structure")
# NATIVE SHOULD JUST BE 3AY4 STRAIGHT FROM PDB
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure.")
parser.add_argument("working_pdb_file", type=str, help="the filename of the PDB structure to be glycosylated.")
parser.add_argument("glyco_file", type=str, help="/path/to/the .iupac glycan file to be used.")
input_args = parser.parse_args()


#################
#### IMPORTS ####
#################

from antibody_functions import *
#from antibody_protocols import *
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file
#GlycanRelaxMover, LinkageConformerMover


##############################
#### PREPARE FOR PROTOCOL ####
##############################

## load up the poses given from the arguments passed
# native pose ( for comparison, really )
native_pose = load_pose( input_args.native_pdb_file )
native_pose_name = "native"
native_pose.pdb_info().name( native_pose_name )

# working pose
working_pose = load_pose( input_args.working_pdb_file )
working_pose_name = "glycosylated_pose"
working_pose.pdb_info().name( working_pose_name )

## get some numbers that will be used in pieces of this protocol
# this number is used later for resetting the core glycan
n_res_no_Fc_glycan = working_pose.n_residue()

# these numbers are of just the receptor glycan - they should be ignored sometimes
# and get the residue numbers of branch points - they should also be ignored sometimes
FcR_seqpos_nums = []
FcR_branch_point_nums = []
for res in working_pose:
    if res.is_carbohydrate():
        FcR_seqpos_nums.append( res.seqpos() )
    if res.is_branch_point():
        FcR_branch_point_nums.append( res.seqpos() )

# get a standard fa_scorefxn for protein stuff
sf = get_fa_scorefxn()

# adjust the standard fa_scorefxn for sugar stuff
sugar_sf = get_fa_scorefxn_with_given_weights( "fa_intra_rep", 0.440 )

# score the unmodified pose using the standard sf
print
print "Unmodified pose:\t\t", sf( working_pose )
print

# pymol stuff
pmm.keep_history(True)
pmm.apply( native_pose )
pmm.apply( working_pose )



# glycosylate the given working pose
# 69 and 284 are the two ASN297 residues from 3ay4
glycosylate_pose_by_file( working_pose, 69, "ND2", input_args.glyco_file )
glycosylate_pose_by_file( working_pose, 284, "ND2", input_args.glyco_file )

#working_pose.pdb_info().name( "glycosylated" )
working_pose.pdb_info().name( working_pose_name )
pmm.apply( working_pose )
print "After glycosylation:\t\t", sf( working_pose )
print


# reset the chibose core ( GlcNAc1 and 2 and Man 3, 4, and 5 )
# don't do reset for G9 and bisecting 2,6-NSCT because wouldn't be as relevant
n_res_Fc_glycan = working_pose.n_residue()
num_sugars_added = n_res_Fc_glycan - n_res_no_Fc_glycan
size_of_one_glycan = num_sugars_added / 2
A_core_GlcNAc = n_res_no_Fc_glycan + 1
B_core_GlcNAc = n_res_no_Fc_glycan + size_of_one_glycan + 1

# numbers collected from NATIVE 3ay4 PDB
# off a bit from the low energy structure given, but that's because I'm sure I'll make more of those
A_phi = -102.659710797
A_psi = 178.686597952
A_omega = -154.56647992437044
B_phi = -84.7881455098
B_psi = 177.132547367
B_omega = -162.5038699839906

# reset the core sugar residues of the glycosylated working_pose
working_pose.set_phi( A_core_GlcNAc, A_phi )
working_pose.set_psi( A_core_GlcNAc, A_psi )
working_pose.set_omega( A_core_GlcNAc, A_omega )

working_pose.set_phi( B_core_GlcNAc, B_phi )
working_pose.set_psi( B_core_GlcNAc, B_psi )
working_pose.set_omega( B_core_GlcNAc, B_omega )

#working_pose.pdb_info().name( "core_sugar_reset" )
working_pose.pdb_info().name( working_pose_name )
pmm.apply( working_pose )
print "After 3-core sugar reset:\t", sf( working_pose )
print



'''
# minimize using a sf with only the sugar_bb term on
sugar_bb_sf = get_sugar_bb_only_sf()
mm = MoveMap()
for residue in working_pose:
    if residue.is_carbohydrate():
        mm.set_bb( residue.seqpos(), True )
min_mover = MinMover( mm, sugar_bb_sf, "dfpmin_strong_wolfe", 0.01, True )
min_mover.apply( working_pose )
pmm.apply( working_pose )
print "After min:\t\t", sf( working_pose )
print
'''



## use the LinkageConformerMover to find a local sugar minima
# first get the sequence positions of the newly glycosylated Fc sugars and their branch points
Fc_sugar_nums = []
Fc_branch_point_nums = []
for res in working_pose:
    if res.is_carbohydrate():
        if res.seqpos() not in FcR_seqpos_nums:
            Fc_sugar_nums.append( res.seqpos() )
    if res.is_branch_point():
        if res.seqpos() not in FcR_branch_point_nums:
            Fc_branch_point_nums.append( res.seqpos() )
        
# make a MoveMap for these Fc sugars allowing only bb movement
mm = make_movemap_for_range( Fc_sugar_nums, allow_bb_movement = True, allow_chi_movement = False )

# add in the branch points myself
for branch_point in Fc_branch_point_nums:
    mm.set_branches( branch_point, True )

            
# TODO: uhhh this LCM thing isn't working in general
lcm = LinkageConformerMover()
lcm.set_movemap( mm )
# TODO: give it an sf with the fa_intra_rep set to .440 (not .004)
# lcm.set_scorefunction( sugar_sf ) ???????
lcm.apply( working_pose )
pmm.apply( working_pose )
print "After LCM:\t\t\t", sf( working_pose )
print


'''
# pack the Fc sugars and around them within 10 Angstroms
pack_rotamers_mover = make_pack_rotamers_mover( sf, working_pose, 
                                                apply_sf_sugar_constraints = False, 
                                                pack_branch_points = True, 
                                                residue_range = Fc_sugar_nums, 
                                                use_pack_radius = True, 
                                                pack_radius = PACK_RADIUS )
pack_rotamers_mover.apply( working_pose )
pmm.apply( working_pose )
print "After Fc sugar and 10"
print "Angstrom sphere pack:\t\t", sf( working_pose )
print
'''



# do a regular pack and minimization round
working_pose = do_pack_min( sf, working_pose, 
                            apply_sf_sugar_constraints = False, 
                            pack_branch_points = True )
#working_pose = do_pack_min( sf, working_pose, 
#                            apply_sf_sugar_constraints = False, 
#                            pack_branch_points = False )
pmm.apply( working_pose )
print "After total pack/min\t\t", sf( working_pose )
print
