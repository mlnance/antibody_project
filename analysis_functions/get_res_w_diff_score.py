from antibody_functions import *

low_E_3ay4_dir = "/Users/Research/new_pyrosetta_git_repo/structures_from_jazz/test_dir/"
base_pdb = "/Users/Research/antibody_project/pdb_copies_dont_touch/lowest_E_single_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb"
base = load_pose( base_pdb )
sf = get_fa_scorefxn()

for mut_name in os.listdir( low_E_3ay4_dir ):
    mut = load_pose( low_E_3ay4_dir + mut_name )

    count, diff_residues = compare_pose_energy_per_residue( sf, base, mut, detailed_analysis = True )
    
    print '*' * 8, mut.pdb_info().name()
    for res in diff_residues:
        print res.name(), base.pdb_info().pose2pdb( res.seqpos() ), "native E:", base.energies().residue_total_energy( res.seqpos() ), "mut E:", mut.energies().residue_total_energy( res.seqpos() )
    print
    print


#base_energy = base.energies()
#mut_energy = mut.energies()

#for ii in range( 1, base.total_residue() + 1 ):
#    diff = mut_energy.residue_total_energy(ii) - base_energy.residue_total_energy(ii)
#    if not -1 < diff < 1:
#        print base_energy.residue_total_energy(ii), mut_energy.residue_total_energy(ii), ii
