#Relevant files
native_crystal_struct_3ay4_Fc_FcgRIII.pdb is the PDB file of 3ay4 with nothing done to it other than making the LINK and HETNAM records correct.<br /><br />
3ay4_interface_loops.txt defines the two loops of the Fc side at the interface of the Fc-FcgRIII.<br /><br />
antibody_functions.py holds all the relevant worker functions that are used in the protocols.<br /><br />
antibody_protocols.py holds the relevant protocols (combination of functions from antibody_functions) that are used in the files described below.<br /><br />

#make_low_E_struct_through_pack_min.py

This takes a pose and uses packing and minimization runs to find a low-energy structure.<br /><br />
This uses the make_base_pack_min_pose function in antibody_protocols, which is one initial pack and minimization, followed by 2 outer and 2 inner rounds of pack min.<br /><br />
I used this on native_crystal_struct_3ay4_Fc_FcgRIII.pdb and made 1000 structures on Louis.<br /><br />
I took the lowest energy structure from that run (which was native_crystal_struct_3ay4_Fc_FcgRIII_956.pdb) to make the mutations from the 2001 Shields paper.<br /><br />
I ran this protocol on each mutant structure to make 1000 low energy structures of each mutant.<br /><br />
I WILL TAKE the lowest energy structure from each mutant run and ANALYZE it to see how effective just packing and minimization can capture the structural and energetic differences seen in the mutants and their binding affinity.

#make_low_E_struct_through_interface_loop_perturb.py
This takes a pose and uses interface loop modeling to find a low-energy structure.<br /><br />
This uses the make_loop_perturbations function in antibody_protocols, which is 100 trials of loop modeling (not 100 per loop, 100 total), followed by a single pack and minimization round.<br /><br />
I WILL use this on the lowest energy structure of 1000 of native_crystal_struct_3ay4_Fc_FcgRIII that has gone through ONE round of pack and minimization.<br /><br />

#make_low_E_struct_through_sample_domain_motions.py
This takes a pose and uses rigid body motions to find a low-energy structure.<br /><br />
This uses the make_rigid_body_moves function from antibody_protocols, which is 100 trials of rigid body moves defined by the jumps in the FoldTree followed by a single pack and minimization round.<br /><br />
I WILL use this on the lowest energy structure of 1000 of native_crystal_struct_3ay4_Fc_FcgRIII that has gone through ONE round of pack and minimization.<br /><br />