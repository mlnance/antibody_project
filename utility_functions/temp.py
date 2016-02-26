#!/usr/bin/python

from antibody_functions import *
from antibody_protocols import *
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file

native_pose_filename = "native_crystal_struct_3ay4_Fc_FcgRIII.pdb"
pose_filename = "crystal_struct_3ay4_Fc_FcgRIII_no_Fc_glycan.pdb"
#pose_filename = "just_ASN.pdb"

loops_filename = "3ay4_interface_loops.txt"

native_pose = load_pose( native_pose_filename )
pose = load_pose( pose_filename )

sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, pose )
print "Unmodified pose", sf( pose )

pmm.keep_history(True)

native_pose.pdb_info().name( "native" )
pmm.apply( native_pose )

pose.pdb_info().name( "orig" )
pmm.apply( pose )


#pose = make_loop_perturbations( loops_filename, sf, pose, trials = 100, model_loops = True, verbose = True, pmm = pmm )
#pose = make_rigid_body_moves( sf, pose, trials = 100, verbose = True, pmm = pmm )
#pose = make_base_pack_min_pose( sf, pose, outer_trials = 2, inner_trials = 2, verbose = True, pmm = pmm )


glyco_file = "/Users/Research/new_pyrosetta_git_repo/database/chemical/carbohydrates/common_glycans/N-glycan_core.iupac"
#glyco_file = "/Users/Research/new_pyrosetta_git_repo/database/chemical/carbohydrates/common_glycans/IgG1_Fc_fucosylated_N-glycan_core.iupac"
glycosylate_pose_by_file( pose, 69, "ND2", glyco_file )
glycosylate_pose_by_file( pose, 284, "ND2", glyco_file )
#glycosylate_pose_by_file( pose, 1, "ND2", glyco_file )

pose.pdb_info().name( "glycosylated" )
pmm.apply( pose )
print "After glycosylation", sf( pose )
print


# get res numbers of carbohydrates
carb_res_nums = []
for res in pose:
    if res.is_carbohydrate():
        carb_res_nums.append( res.seqpos() )


# get res numbers of carbohydrates in native
native_carb_res_nums = []
for res in native_pose:
    if res.is_carbohydrate():
        native_carb_res_nums.append( res.seqpos() )


## native to new dictionaries
# native to N-glycan_core.iupac on crystal
native_to_N_glycan_core = { 216: 603, 217: 604, 218: 605, 219: 606, 221: 607, 440: 608, 441: 609, 442: 610, 443: 611, 445: 612}


# reset the phi and psi values
dictionary = native_to_N_glycan_core

for key in dictionary:
    pose.set_phi( native_to_N_glycan_core[ key ], native_pose.phi( key ) )
    pose.set_psi( native_to_N_glycan_core[ key ], native_pose.psi( key ) )
pose.pdb_info().name( "reset_phi_psi" )
pmm.apply( pose )
print "After phi and psi reset", sf( pose )
print


# pack everything except branch points
# use a standard fa ScoreFunction
packer_mover = make_pack_rotamers_mover( sf, pose, apply_sf_sugar_constraints = False, pack_branch_points = False, verbose = True )
pose.pdb_info().name( "no_branch_pack" )
packer_mover.apply( pose )
pmm.apply( pose )
print "After no branch pack", sf( pose )
print


# minimize using a sf with only the sugar_bb term on
sugar_bb_sf = get_sugar_bb_only_sf()
mm = MoveMap()
for residue in pose:
    if residue.is_carbohydrate():
        mm.set_bb( residue.seqpos(), True )
min_mover = MinMover( mm, sugar_bb_sf, "dfpmin_strong_wolfe", 0.01, True )
min_mover.apply( pose )
pose.pdb_info().name( "first_min" )
pmm.apply( pose )
print "After min", sf( pose )
print


# do a regular pack and minimization round
pose = do_pack_min( sf, pose, apply_sf_sugar_constraints = False, pack_branch_points = True, verbose = True )
pose.pdb_info().name( "first_pack_min" )
pmm.apply( pose )
print "After first pack_min", sf( pose )
print

'''
# another pack and minimization
pose = do_pack_min( sf, pose, apply_sf_sugar_constraints = False, pack_branch_points = True, verbose = True )
pose.pdb_info().name( "second_pack_min" )
pmm.apply( pose )
print "After second pack_min", sf( pose )
print


# yet another pack and minimization
pose = do_pack_min( sf, pose, apply_sf_sugar_constraints = True, pack_branch_points = True, verbose = True )
pose.pdb_info().name( "third_pack_min" )
pmm.apply( pose )
print "After third pack_min", sf( pose )
print
'''
