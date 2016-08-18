#!/usr/bin/python

from antibody_functions import *
from antibody_protocols import *
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file

native_pose_filename = "pdb_copies_dont_touch/native_crystal_struct_3ay4_Fc_FcgRIII.pdb"

loops_filename = "3ay4_interface_loops.txt"

native_pose = load_pose( native_pose_filename )
pose_name = native_pose.pdb_info().name()

sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, native_pose )
print "Native pose", sf( native_pose )
print

pmm.keep_history(True)

native_pose.pdb_info().name( "native" )
pmm.apply( native_pose )


#pose = make_loop_perturbations( loops_filename, sf, pose, trials = 100, model_loops = True, verbose = True, pmm = pmm )
#pose = make_rigid_body_moves( sf, pose, trials = 100, verbose = True, pmm = pmm )
#pose = make_base_pack_min_pose( sf, pose, outer_trials = 2, inner_trials = 2, verbose = True, pmm = pmm )


# pack and min
pose = do_pack_min( sf, pose, apply_sf_sugar_constraints = True, pack_branch_points = True, verbose = True )
pose.pdb_info().name( "pack_min" )
pmm.apply( pose )
print "After pack_min", sf( pose )
print
