#!/usr/bin/python

from antibody_functions import initialize_rosetta, \
    load_pose, set_glycan_to_ideal_SugarBB_phi_psi, \
    native_Fc_glycan_nums_except_core_GlcNAc, \
    get_ideal_SugarBB_phi_psi_info, \
    set_3ay4_Fc_glycan_except_core_GlcNAc_to_ideal_LCM_phi_psi_omega, \
    get_3ay4_ideal_LCM_phi_psi_info
from rosetta import PyMOL_Mover
import pandas as pd
import matplotlib.pyplot as plt


initialize_rosetta()

pmm = PyMOL_Mover()

native = load_pose("project_structs/fa_intra_rep_lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb")
#native.pdb_info().name( "native" )
#pmm.apply( native )

native_df = pd.DataFrame()
native_res = []
native_phi = []
native_psi = []
for res in native_Fc_glycan_nums_except_core_GlcNAc:
    native_res.append( res )
    native_phi.append( native.phi( res ) )
    native_psi.append( native.psi( res ) )
native_df[ "res_num" ] = native_res
native_df[ "phi" ] = native_phi
native_df[ "psi" ] = native_psi



# data pulled from the RosettaCarbohydrates Tutorial/Demo 2 paper from Jason (the research paper, specifically)
#ideal_SugarBB_native = native.clone()
#ideal_SugarBB_native.assign( set_glycan_to_ideal_SugarBB_phi_psi( native_Fc_glycan_nums_except_core_GlcNAc, ideal_SugarBB_native ) )
#ideal_SugarBB_native.pdb_info().name( "ideal_SugarBB_native" )
#pmm.apply( ideal_SugarBB_native )
SugarBB_df = pd.DataFrame()
SugarBB_res = []
SugarBB_phi = []
SugarBB_psi = []
for res in native_Fc_glycan_nums_except_core_GlcNAc:
    # returns list( ideal_phi_major, ideal_phi_stdev_major, ideal_phi_minor, ideal_phi_stdev_minor, ideal_psi, ideal_psi_stdev, anomeric_position, linkage_type )
    phi_psi_data = get_ideal_SugarBB_phi_psi_info( res, native )
    phi = phi_psi_data[0]
    psi = phi_psi_data[4]
    SugarBB_res.append( res )
    SugarBB_phi.append( phi )
    SugarBB_psi.append( psi )
SugarBB_df[ "res_num" ] = SugarBB_res
SugarBB_df[ "phi" ] = SugarBB_phi
SugarBB_df[ "psi" ] = SugarBB_psi



# data pulled from database/chemical/carbohydrates/linkage_conformers/default.table
#ideal_LCM_native = native.clone()
#ideal_LCM_native.assign( set_3ay4_Fc_glycan_except_core_GlcNAc_to_ideal_LCM_phi_psi_omega( ideal_LCM_native ) )
#ideal_LCM_native.pdb_info().name( "ideal_LCM_native" )
#pmm.apply( ideal_LCM_native )
LCM_phi_dict, phi_stdev, LCM_psi_dict, psi_stdev, omega_data, omega_stdev = get_3ay4_ideal_LCM_phi_psi_info()    
LCM_df = pd.DataFrame()
LCM_res = []
LCM_phi = []
LCM_psi = []
keys = LCM_phi_dict.keys()
keys.sort()
for key in keys:
    LCM_res.append( key )
    LCM_phi.append( LCM_phi_dict[ key ] )
    LCM_psi.append( LCM_psi_dict[ key ] )
LCM_df[ "res_num" ] = LCM_res
LCM_df[ "phi" ] = LCM_phi
LCM_df[ "psi" ] = LCM_psi



fig = plt.figure(figsize=(12,12))
for ii in range( 1, len( keys ) + 1 ):
    ax = fig.add_subplot( 7, 2, ii )
    res_num = str( native.pdb_info().pose2pdb( native_df["res_num"][ ii - 1 ] ) ).strip()
    ax.set_title( "Residue %s" %res_num )
    #ax.set_title( "Residue %s" %native_df["res_num"][ ii - 1 ] )

    sc = plt.scatter( native_df["phi"][ ii - 1 ], native_df["psi"][ ii - 1 ], marker='o', c="red", s=50 )
    sc = plt.scatter( LCM_df["phi"][ ii - 1 ], LCM_df["psi"][ ii - 1 ], marker='o', c="green", s=50 )
    # some SugarBB residues (mainly the res that connects the Man branch at 6) don't have psi
    if SugarBB_df["psi"][ ii - 1 ] != "NA":
        sc = plt.scatter( SugarBB_df["phi"][ ii - 1 ], SugarBB_df["psi"][ ii - 1 ], marker='o', c="blue", s=50 )
    plt.xlim( [ -190, 190 ] )
    plt.xlabel( "phi" )
    plt.ylim( [ -190, 190 ] )
    plt.ylabel( "psi" )

plt.tight_layout()
plot_title = "Native vs SugarBB vs ideal LCM phi and psi"
plt.suptitle( plot_title, fontsize = 28 )
plt.subplots_adjust(top=0.93)
plt.savefig( plot_title, dpi=120, transparent=True )
plt.close()
