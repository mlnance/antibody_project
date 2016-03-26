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
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure.")
parser.add_argument("working_pdb_file", type=str, help="the filename of the PDB structure to be glycosylated.")
parser.add_argument("glyco_file", type=str, help="/path/to/the .iupac glycan file to be used.")
input_args = parser.parse_args()


#################
#### IMPORTS ####
#################

from antibody_functions import *
from rosetta import MonteCarlo
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file
from rosetta.protocols.carbohydrates import GlycanRelaxMover, LinkageConformerMover
import os, sys
sys.path.append( "utility_functions" )
from nearby_residues_to_pickle_file import main as get_nearby_residues


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



# reset the 1st GlcNAc on the ASN in the chibose core
n_res_Fc_glycan = working_pose.n_residue()
num_sugars_added = n_res_Fc_glycan - n_res_no_Fc_glycan
size_of_one_glycan = num_sugars_added / 2
A_core_GlcNAc = n_res_no_Fc_glycan + 1
B_core_GlcNAc = n_res_no_Fc_glycan + size_of_one_glycan + 1

# numbers collected from the lowest_E decoy of 3ay4 PDB after just a total pack/min
A_phi = -102.58846630984607
A_psi = 178.68359502878405
A_omega = -153.94141254167278
B_phi = -84.79653728190505
B_psi = 177.12713080144076
B_omega = -156.4554337951647

# reset both of the core GlcNAc residue of the glycosylated working_pose
working_pose.set_phi( A_core_GlcNAc, A_phi )
working_pose.set_psi( A_core_GlcNAc, A_psi )
working_pose.set_omega( A_core_GlcNAc, A_omega )

working_pose.set_phi( B_core_GlcNAc, B_phi )
working_pose.set_psi( B_core_GlcNAc, B_psi )
working_pose.set_omega( B_core_GlcNAc, B_omega )

#working_pose.pdb_info().name( "core_sugar_reset" )
working_pose.pdb_info().name( working_pose_name )
pmm.apply( working_pose )
print "After reseting the"
print "core GlcNAc:\t\t\t", sf( working_pose )
print



# get the res nums of the Fc sugars added
Fc_sugar_nums = []
Fc_branch_point_nums = []
for res in working_pose:
    if res.is_carbohydrate():
        if res.seqpos() not in FcR_seqpos_nums:
            Fc_sugar_nums.append( res.seqpos() )
    if res.is_branch_point():
        if res.seqpos() not in FcR_branch_point_nums:
            Fc_branch_point_nums.append( res.seqpos() )
# get a list of the Fc sugars discluding the core GlcNAc residues
Fc_sugar_nums_except_core_GlcNAc = []
for res in Fc_sugar_nums:
    if res != A_core_GlcNAc and res != B_core_GlcNAc:
        Fc_sugar_nums_except_core_GlcNAc.append( res )




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



## use the LinkageConformerMover to find a local sugar minima        
# make a MoveMap for these Fc sugars allowing only bb movement
mm = make_movemap_for_range( Fc_sugar_nums_except_core_GlcNAc, allow_bb_movement = True, allow_chi_movement = False )

# TODO: Ask Jason if I should allow branch point movement
'''
# add in the branch points myself
for branch_point in Fc_branch_point_nums:
    mm.set_branches( branch_point, True )
'''

# make an appropriate MonteCarlo object
# kT is 0.7 - from antibody_functions.py
mc = MonteCarlo( working_pose, sugar_sf, kT )

# make an appropriate LinkageConformerMover
lcm = LinkageConformerMover()
lcm.set_movemap( mm )
lcm.set_x_standard_deviations( 2 )

# run the LCM 10-100 times using a MonteCarlo object to accept or reject the move
num_lcm_accept = 0
for ii in range( 100 ):
    # apply the LCM
    lcm.apply( working_pose )
    
    # accept or reject the move using the MonteCarlo object
    if mc.boltzmann( working_pose ):
        num_lcm_accept += 1
        pmm.apply( working_pose )

# pack the Fc sugars and around them within 10 Angstroms
pack_rotamers_mover = make_pack_rotamers_mover( sf, working_pose, 
                                                apply_sf_sugar_constraints = False, 
                                                pack_branch_points = True, 
                                                residue_range = Fc_sugar_nums, 
                                                use_pack_radius = True, 
                                                pack_radius = PACK_RADIUS )
pack_rotamers_mover.apply( working_pose )
pmm.apply( working_pose )
print "After LCM and a 10 Ang" 
print "sphere pack/min:\t\t", sf( working_pose )
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
'''
