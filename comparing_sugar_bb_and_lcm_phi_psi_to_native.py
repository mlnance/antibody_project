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

# skipping the beta phi minor as none of the native sugars are in that range anyway
df_columns = [ "res_num", 
               "actual_sugar_bb", 
               "anomeric_position", 
               "actual_phi", 
               "ideal_phi_major", 
               "ideal_phi_stdev_major", 
               "delta_phi_major", 
               "linkage_type", 
               "actual_psi", 
               "ideal_psi", 
               "ideal_psi_stdev", 
               "delta_psi" ]
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


    # since there won't always be an ideal psi
    try:
        delta_psi = native.psi( glyc_num ) - ideal_psi
    except:
        delta_psi = "NA"


    data_row = [ int( glyc_num ), # res num
                 round( native.energies().residue_total_energies( glyc_num ).get( score_type_from_name( "sugar_bb" ) ), 3 ), # actual sugar_bb
                 anomeric_position, 
                 native.phi( glyc_num ), # actual phi
                 ideal_phi_major, 
                 ideal_phi_stdev_major, 
                 native.phi( glyc_num ) - ideal_phi_major, # delta phi major
                 linkage_type, 
                 native.psi( glyc_num ), # actual psi
                 ideal_psi, 
                 ideal_psi_stdev, 
                 delta_psi ]
                 
    df.loc[ ii ] = data_row
    ii += 1
