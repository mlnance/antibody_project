#!/usr/bin/python

import argparse

# parse and store input arguments
parser = argparse.ArgumentParser(description="Use PyRosetta to glycosylate a pose and find a low E structure")
parser.add_argument("native_pdb_file", type=str, help="the filename of the native PDB structure.")
parser.add_argument("working_pdb_file", type=str, help="the filename of the PDB structure to be glycosylated.")
parser.add_argument("glyco_file", type=str, help="/path/to/the .iupac glycan file to be used.")
input_args = parser.parse_args()


from antibody_functions import *
from antibody_protocols import *
from rosetta.core.pose.carbohydrates import glycosylate_pose_by_file

# load up the poses given
native_pose = load_pose( input_args.native_pdb_file )
native_pose.pdb_info().name( "native" )
working_pose = load_pose( input_args.working_pdb_file )
working_pose_name = "glycosylated_pose"
working_pose.pdb_info().name( working_pose_name )

# this is used later for resetting the core glycan
n_res_no_Fc_glycan = working_pose.n_residue()

# don't add sugar constraints?
sf = get_fa_scorefxn()
print "Unmodified pose", sf( working_pose )
print

# pymol stuff
pmm.keep_history(True)
pmm.apply( native_pose )
pmm.apply( working_pose )


## glyco files of interest
#glyco_file = "/Users/Research/pyrosetta_dr/database/chemical/carbohydrates/common_glycans/bisected_fucosylated_N-glycan_core.iupac"
#glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/N-glycan_core.iupac"
#glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/3ay4_Fc_Glycan.iupac"
#glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/2_6-NSCT_CW_Lin.iupac"
#glyco_file = "/Users/Research/pyrosetta_dir/database/chemical/carbohydrates/common_glycans/G4_CW_Lin.iupac"


# glycosylate the given working pose
# 69 and 284 are the two ASN297 residues from 3ay4
glycosylate_pose_by_file( working_pose, 69, "ND2", input_args.glyco_file )
glycosylate_pose_by_file( working_pose, 284, "ND2", input_args.glyco_file )

#working_pose.pdb_info().name( "glycosylated" )
working_pose.pdb_info().name( working_pose_name )
pmm.apply( working_pose )
print "After glycosylation", sf( working_pose )
print



# reset the chibose core (GlcNAc1 and 2 and Man 3, 4, and 5)
# don't do reset for G9 and bisecting 2,6-NSCT because wouldn't be as relevant
n_res_Fc_glycan = working_pose.n_residue()
num_sugars_added = n_res_Fc_glycan - n_res_no_Fc_glycan
size_of_one_glycan = num_sugars_added / 2
A_core_GlcNAc = n_res_no_Fc_glycan + 1
B_core_GlcNAc = n_res_no_Fc_glycan + size_of_one_glycan + 1

# numbers collected from native 3ay4 PDB
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

working_pose.pdb_info().name( working_pose_name )
pmm.apply( working_pose )
print "After core GlcNAc reset", sf( working_pose )
print



# use the Glycan Relax Mover to find a local sugar minima
grm = GlycanRelaxMover()
grm.apply( working_pose )
pmm.apply( working_pose )
print "After GRM", sf( working_pose )


'''
# pack around the Fc sugars
Fc_res_range = range( n_res_no_Fc_glycan + 1, n_res_Fc_glycan + 1 )
pack_rotamers_mover = make_pack_rotamers_mover( sf, working_pose, False, True, Fc_res_range, True, PACK_RADIUS, True )
pack_rotamers_mover.apply( working_pose )
pmm.apply( working_pose )
print "After packing of Fc sugar and within 10 A of it", sf( working_pose )
print
'''



# do a regular pack and minimization round
working_pose = do_pack_min( sf, working_pose, apply_sf_sugar_constraints = False, pack_branch_points = True, verbose = True )
pmm.apply( working_pose )
print "After first pack_min", sf( working_pose )
print



'''
# pack everything except branch points
# use a standard fa ScoreFunction
packer_mover = make_pack_rotamers_mover( sf, working_pose, apply_sf_sugar_constraints = False, pack_branch_points = False, verbose = True )
packer_mover.apply( working_pose )
#working_pose.pdb_info().name( "no_branch_pack" )
pmm.apply( working_pose )
print "After no branch pack", sf( working_pose )
print
'''



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
print "After min", sf( working_pose )
print
'''
