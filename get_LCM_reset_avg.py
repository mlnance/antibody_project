#!/usr/bin/python

from antibody_functions import initialize_rosetta, \
    PyMOL_Mover, load_pose, native_Fc_glycan_nums_except_core_GlcNAc, \
    get_fa_scorefxn_with_given_weights, calc_avg_of_list
from rosetta import MoveMap
from rosetta.protocols.carbohydrates import LinkageConformerMover
import pandas as pd


initialize_rosetta()

#pmm = PyMOL_Mover()
#pmm.keep_history( True )

native = load_pose("project_structs/fa_intra_rep_lowest_E_double_pack_and_min_only_native_crystal_struct_3ay4_Fc_FcgRIII.pdb")

sf = get_fa_scorefxn_with_given_weights( {"fa_intra_rep":0.44} )
sf( native )

phi_data_dict = {}
phi_avg_data_dict = {}
psi_data_dict = {}
psi_avg_data_dict = {}
omega_data_dict = {}
omega_avg_data_dict = {}
for ii in range( 1, 1000 + 1 ):
    testing_pose = native.clone()
    if ii % 10 == 0:
        print ii

    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        res_mm = MoveMap()
        res_mm.set_bb( res_num, True )
        res_mm.set_chi( res_num, True )

        lcm = LinkageConformerMover()
        lcm.set_movemap( res_mm )
        lcm.set_x_standard_deviations( 2 )

        lcm.apply( testing_pose )

        # data
        if res_num in phi_data_dict.keys():
            phi_data_dict[ res_num ].append( testing_pose.phi( res_num ) )
        else:
            phi_data_dict[ res_num ] = []
        if res_num in psi_data_dict.keys():
            psi_data_dict[ res_num ].append( testing_pose.psi( res_num ) )
        else:
            psi_data_dict[ res_num ] = []
        if res_num in omega_data_dict.keys():
            omega_data_dict[ res_num ].append( testing_pose.omega( res_num ) )
        else:
            omega_data_dict[ res_num ] = []

for res_num in phi_data_dict.keys():
    phi_avg_data_dict[ res_num ] = calc_avg_of_list( phi_data_dict[ res_num ] )
for res_num in psi_data_dict.keys():
    psi_avg_data_dict[ res_num ] = calc_avg_of_list( psi_data_dict[ res_num ] )
for res_num in omega_data_dict.keys():
    omega_avg_data_dict[ res_num ] = calc_avg_of_list( omega_data_dict[ res_num ] )
