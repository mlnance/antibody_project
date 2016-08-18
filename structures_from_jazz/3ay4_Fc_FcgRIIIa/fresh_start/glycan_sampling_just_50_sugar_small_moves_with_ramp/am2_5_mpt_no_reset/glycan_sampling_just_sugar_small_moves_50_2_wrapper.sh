#!/bin/sh

cd /home/mlnance
source /home/mlnance/SetPyRosettaEnvironment.sh

native_pdb_file="/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/fresh_start/glycan_sampling_just_50_sugar_small_moves_with_ramp/am2_5_mpt_no_reset/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb"

glyco_file="/home/mlnance/project_glyco_files/3ay4_Fc_Glycan.iupac"
utility_files_dir="/home/mlnance/project_utility_files"
structure_dir="/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/fresh_start/glycan_sampling_just_50_sugar_small_moves_with_ramp/am2_5_mpt_no_reset"

nstruct=1000
num_sugar_small_moves_trials=50
num_sugar_small_moves_per_trial=5

angle_multiplier="--angle_multiplier=2"
ramp_sf="--ramp_sf"


/home/mlnance/project_created_structs/3ay4_Fc_FcgRIIIa/fresh_start/glycan_sampling_just_50_sugar_small_moves_with_ramp/am2_5_mpt_no_reset/./glycan_sampling_just_sugar_small_moves_on_native.py $native_pdb_file $glyco_file $utility_files_dir $structure_dir $nstruct $num_sugar_small_moves_trials $num_sugar_small_moves_per_trial $ramp_sf $angle_multiplier
