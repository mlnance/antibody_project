#!/usr/bin/python

from antibody_functions import load_pose, initialize_rosetta, \
    decoy_to_native_res_map, decoy_Fc_glycan
from rosetta import PyMOL_Mover, Pose
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file
import sys
sys.path.append( "/Users/Research/pyrosetta_dir/project_utility_files" )
from pose_metrics_util import Fc_glycan_rmsd


initialize_rosetta()

pmm = PyMOL_Mover()
pmm.keep_history( True )

#working = load_pose( "just_native_sugar.pdb" )
native = load_pose( "/Users/Research/pyrosetta_dir/project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb" )
native.pdb_info().name( "native" )
#pmm.apply( native )

#working = load_pose( "just_native_ASN.pdb" )
working = load_pose( "/Users/Research/pyrosetta_dir/project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII_removed_Fc_sugar.pdb" )
working.pdb_info().name( "working" )

decoy = Pose()
decoy.assign( working )
decoy.pdb_info().name( "decoy" )
#pmm.apply( decoy )


glycosylate_these_ASN = [ 69, 284 ]
for ASN in glycosylate_these_ASN:
    glycosylate_pose_by_file( decoy, ASN, "ND2", "/Users/Research/antibody_project/send_to_louis/project_glyco_files/3ay4_Fc_Glycan.iupac" )
decoy.pdb_info().name( "decoy" )
#pmm.apply( decoy )


for decoy_res in glycosylate_these_ASN:
    native_res = decoy_to_native_res_map[ decoy_res ]
    print decoy_res, native_res
    decoy.set_chi( 1, decoy_res, native.chi( 1, native_res ) )
    decoy.set_chi( 2, decoy_res, native.chi( 2, native_res ) )
    decoy.set_phi( decoy_res, native.phi( native_res ) )
    decoy.set_psi( decoy_res, native.psi( native_res ) )
    decoy.set_omega( decoy_res, native.omega( native_res ) )

for decoy_res in decoy_Fc_glycan:
    native_res = decoy_to_native_res_map[ decoy_res ]
    print decoy_res, native_res
    decoy.set_phi( decoy_res, native.phi( native_res ) )
    decoy.set_psi( decoy_res, native.psi( native_res ) )
    decoy.set_omega( decoy_res, native.omega( native_res ) )
pmm.apply( decoy )


print Fc_glycan_rmsd( decoy, native, [ 'H', 'I', 'J', 'K' ], [ 'D', 'E', 'F', 'G' ], 1, "/Users/Research/pyrosetta_dir" )

