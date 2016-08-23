#!/usr/bin/python

from antibody_functions import initialize_rosetta, \
    PyMOL_Mover, load_pose, native_Fc_glycan_nums_except_core_GlcNAc
from antibody_functions import get_fa_scorefxn_with_given_weights, \
    get_ideal_SugarBB_phi_psi_info
from rosetta.core.scoring import score_type_from_name
import pandas as pd


initialize_rosetta()

#pmm = PyMOL_Mover()
#pmm.keep_history( True )

native = load_pose("project_structs/lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb")

sf = get_fa_scorefxn_with_given_weights( {"fa_intra_rep":0.44} )
sf( native )

# data collected from get_LCM_reset_avg.py
# glyc_num : average phi value from 1000 rounds of LCM reset on native
LCM_reset_phi_data_dict = {
    217: -75.57862704304912,
    218: -86.34757667121097,
    219: 71.65810626538641,
    220: -54.07113861114163,
    221: 64.49087672483115,
    222: -55.88631275286347,
    223: -71.68074664418894,
    441: -76.14390058548464,
    442: -86.07066156578118,
    443: 71.62491060257564,
    444: -56.21335685303532,
    445: 64.67063326321038,
    446: -54.8212208144709,
    447: -71.3670738332849 }

# glyc_num : average psi value from 1000 rounds of LCM reset on native
LCM_reset_psi_data_dict = {
    217: 119.06380191080791,
    218: 110.5856808525849,
    219: -120.39653090301904,
    220: -95.50588748028322,
    221: 15.904805387579117,
    222: -95.67118504137507,
    223: 132.4349256412706,
    441: 119.20836934592369,
    442: 110.98888813099123,
    443: -120.70893155147938,
    444: -96.01967601468648,
    445: 10.77809135695188,
    446: -95.66204572169633,
    447: 132.0304643911755 }

# glyc_num : average omega value from 1000 rounds of LCM reset on native
LCM_reset_omega_data_dict = {
    221: -18.691274994795254,
    445: -28.937828470613642 }


# skipping the beta phi minor as none of the native sugars are in that range anyway
#"native_within_ideal_phi_stdev_major", 
#"LCM_reset_within_phi_stdev_major", 
df_columns = [ "res_num", 
               "actual_sugar_bb", 
               "anomeric_position", 
               "actual_phi", 
               "ideal_phi_major", 
               "ideal_phi_stdev_major", 
               "delta_phi_major", 
               "LCM_phi_reset_avg_this_res", 
               "delta_LCM_phi_reset_from_native", 
               "delta_LCM_phi_reset_from_ideal_major", 
               "linkage_type", 
               "actual_psi", 
               "ideal_psi", 
               "ideal_psi_stdev", 
               "delta_psi", 
               "LCM_psi_reset_avg_this_res", 
               "delta_LCM_psi_reset_from_native", 
               "delta_LCM_psi_reset_from_ideal" ]
df = pd.DataFrame( columns = df_columns )

ii = 0
for glyc_num in native_Fc_glycan_nums_except_core_GlcNAc:
    # data returned as [ ideal_phi_major, ideal_phi_stdev_major, ideal_phi_minor, ideal_phi_stdev_minor, ideal_psi, ideal_psi_stdev, anomeric_position, linkage_type ]
    glyc_data = get_ideal_SugarBB_phi_psi_info( glyc_num, native )
    ideal_phi_major = glyc_data[ 0 ]
    ideal_phi_stdev_major = glyc_data[ 1 ]
    ideal_phi_minor = glyc_data[ 2 ]
    ideal_phi_stdev_minor = glyc_data[ 3 ]
    ideal_psi = glyc_data[ 4 ]
    ideal_psi_stdev = glyc_data[ 5 ]
    anomeric_position = glyc_data[ 6 ]
    linkage_type = glyc_data[ 7 ]


    # pull out LCM reset data
    LCM_reset_phi_avg_this_res = LCM_reset_phi_data_dict[ glyc_num ]
    LCM_reset_psi_avg_this_res = LCM_reset_psi_data_dict[ glyc_num ]

    # since there won't always be an ideal psi
    try:
        delta_psi = native.psi( glyc_num ) - ideal_psi
    except:
        delta_psi = "NA"
    try:
        delta_LCM_psi_reset_from_ideal = LCM_reset_psi_avg_this_res - ideal_psi
    except:
        delta_LCM_psi_reset_from_ideal = "NA"


    data_row = [ int( glyc_num ), # res num
                 round( native.energies().residue_total_energies( glyc_num ).get( score_type_from_name( "sugar_bb" ) ), 3 ), # actual sugar_bb
                 anomeric_position, 
                 native.phi( glyc_num ), # actual phi
                 ideal_phi_major, 
                 ideal_phi_stdev_major, 
                 native.phi( glyc_num ) - ideal_phi_major, # delta phi major
                 LCM_reset_phi_avg_this_res, 
                 native.phi( glyc_num ) - LCM_reset_phi_avg_this_res, # the diff b/w native and avg LCM phi reset
                 LCM_reset_phi_avg_this_res - ideal_phi_major, # the diff b/w avg LCM phi reset and ideal phi major
                 linkage_type, 
                 native.psi( glyc_num ), # actual psi
                 ideal_psi, 
                 ideal_psi_stdev, 
                 delta_psi, 
                 LCM_reset_psi_avg_this_res, 
                 native.psi( glyc_num ) - LCM_reset_psi_avg_this_res, # the diff b/w native and avg LCM psi reset
                 delta_LCM_psi_reset_from_ideal ] # the diff b/w avg LCM psi reset and ideal psi
                 
    df.loc[ ii ] = data_row
    ii += 1
