#!/usr/bin/python

from antibody_functions import *
from antibody_protocols import *
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file

native_pose_filename = "pdb_copies_dont_touch/native_crystal_struct_3ay4_Fc_FcgRIII.pdb"
#pose_filename = "pdb_copies_dont_touch/crystal_struct_3ay4_Fc_FcgRIII_no_Fc_glycan.pdb"
pose_filename = "just_ASN.pdb"

loops_filename = "3ay4_interface_loops.txt"

native_pose = load_pose( native_pose_filename )
pose = load_pose( pose_filename )
pose_name = pose.pdb_info().name()
n_res_no_Fc_glycan = pose.n_residue()

sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, pose )
print "Unmodified pose", sf( pose )
print

pmm.keep_history(True)

#native_pose.pdb_info().name( "native" )
#pmm.apply( native_pose )

#pose.pdb_info().name( "orig" )
pmm.apply( pose )


#pose = make_loop_perturbations( loops_filename, sf, pose, trials = 100, model_loops = True, verbose = True, pmm = pmm )
#pose = make_rigid_body_moves( sf, pose, trials = 100, verbose = True, pmm = pmm )
#pose = make_base_pack_min_pose( sf, pose, outer_trials = 2, inner_trials = 2, verbose = True, pmm = pmm )


#glyco_file = "/Users/Research/pyrosetta_dr/database/chemical/carbohydrates/common_glycans/bisected_fucosylated_N-glycan_core.iupac"
#glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/N-glycan_core.iupac"
glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/bisected_N-glycan_core.iupac"
#glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/2_6-NSCT_CW_Lin.iupac"
#glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/G4_CW_Lin.iupac"

#glycosylate_pose_by_file( pose, 69, "ND2", glyco_file )
#glycosylate_pose_by_file( pose, 284, "ND2", glyco_file )
glycosylate_pose_by_file( pose, 1, "ND2", glyco_file )

#pose.pdb_info().name( "glycosylated" )
pose.pdb_info().name( pose_name )
pmm.apply( pose )
print "After glycosylation", sf( pose )
print


'''
# update the phi and psi values of the N-glycan core
# numbers retrived from native 3ay4 pose
core_Glc1_phi = -102.65971079696487
core_Glc1_psi = 178.68659795163614
core_Glc2_phi = -88.56117917014228
core_Glc2_psi = 92.61694177698288
core_Man3_phi = -95.3254524164334
core_Man3_psi = 90.98902629335966
core_Man4_phi = 68.73593249372874
core_Man4_psi = -118.92435416833992
core_Man1_phi = 65.35885873542475
core_Man1_psi = 174.8356913294995

# different dictionaries depending on the sugar being used
# { residue number : core_sugar }
A_N_glycan_core_phi = { 603 : core_Glc1_phi, 604 : core_Glc2_phi,
                      605 : core_Man3_phi, 606 : core_Man4_phi,
                      607 : core_Man1_phi }
A_N_glycan_core_psi = { 603 : core_Glc1_psi, 604 : core_Glc2_psi,
                      605 : core_Man3_psi, 606 : core_Man4_psi,
                      607 : core_Man1_psi }

B_N_glycan_core_phi = { 608 : core_Glc1_phi, 609 : core_Glc2_phi,
                      610 : core_Man3_phi, 611 : core_Man4_phi,
                      612 : core_Man1_phi }
B_N_glycan_core_psi = { 608 : core_Glc1_psi, 609 : core_Glc2_psi,
                      610 : core_Man3_psi, 611 : core_Man4_psi,
                      612 : core_Man1_psi }

n_res_Fc_glycan = pose.n_residue()
num_sugars_added = n_res_Fc_glycan - n_res_no_Fc_glycan

# fix phi and psi on one side
for res_num in range( n_res_no_Fc_glycan + 1, n_res_no_Fc_glycan + ( num_sugars_added / 2 ) + 1 ):
    pose.set_phi( res_num, A_N_glycan_core_phi[ res_num ] )
    pose.set_psi( res_num, A_N_glycan_core_psi[ res_num ] )
# fix phi and psi on the other side
for res_num in range( n_res_no_Fc_glycan + ( num_sugars_added / 2 ) + 1, n_res_Fc_glycan + 1 ):
    pose.set_phi( res_num, B_N_glycan_core_phi[ res_num ] )
    pose.set_psi( res_num, B_N_glycan_core_psi[ res_num ] )

pose.pdb_info().name( pose_name )
pmm.apply( pose )
print "After phi/psi reset", sf( pose )
print
'''


grm = GlycanRelaxMover()
grm.apply( pose )
pmm.apply( pose )
print "After GRM", sf( pose )
print


'''
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
native_to_N_glycan_core = { 216: 603, 217: 604, 218: 605, 219: 606, 221: 607, 440: 608, 441: 609, 442: 610, 443: 611, 445: 612 }


# reset the phi and psi values
dictionary = native_to_N_glycan_core

for key in dictionary:
    pose.set_phi( native_to_N_glycan_core[ key ], native_pose.phi( key ) )
    pose.set_psi( native_to_N_glycan_core[ key ], native_pose.psi( key ) )
#pose.pdb_info().name( "reset_phi_psi" )
pmm.apply( pose )
print "After phi and psi reset", sf( pose )
print
'''

'''
# pack everything except branch points
# use a standard fa ScoreFunction
packer_mover = make_pack_rotamers_mover( sf, pose, apply_sf_sugar_constraints = False, pack_branch_points = False, verbose = True )
packer_mover.apply( pose )
#pose.pdb_info().name( "no_branch_pack" )
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
pmm.apply( pose )
print "After min", sf( pose )
print


# do a regular pack and minimization round
pose = do_pack_min( sf, pose, apply_sf_sugar_constraints = False, pack_branch_points = True, verbose = True )
pmm.apply( pose )
print "After first pack_min", sf( pose )
print


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
