#!/bin/sh

cd /home/mlnance
source /home/mlnance/SetPyRosettaEnvironment.sh

native_pdb_file="/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/fresh_start/get_base_pack_min_pose/std_sf/native_crystal_struct_3ay4_Fc_FcgRIII.pdb"
structure_dir="/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/fresh_start/get_base_pack_min_pose/std_sf"
utility_files_dir="/home/mlnance/project_utility_files"
base_nstruct=1000


/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/fresh_start/get_base_pack_min_pose/std_sf/./get_base_pack_min_pose.py $native_pdb_file $structure_dir $utility_files_dir $base_nstruct
