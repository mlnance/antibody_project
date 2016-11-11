# instantiate the proper Protocol_X object
from native_3ay4_glycan_modeling_protocol import Model3ay4Glycan
if input_args.protocol_num == 0:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_0 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_0 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_5A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 1:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_1 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_1 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_5A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 2:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_2 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_2 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 3:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_3 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_3 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 100
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = True
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 4:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_4 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_4 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 100
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 5:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_5 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_5 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.5
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 6:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_6 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_6 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 7:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_7 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, False )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_7 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose
    GlycanModelProtocol.dump_reset_pose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 8:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_8 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_8 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 5
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )


elif input_args.protocol_num == 9:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_9 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_9 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 100
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 10:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_10 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_10 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 11:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_11 version
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_11 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose
    GlycanModelProtocol.dump_reset_pose = True

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 12:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_12 version
    #######################################################################
    #### CORE GlcNAc CAN MOVE IN THIS PROTOCOL BUT IT STARTS AT NATIVE ####
    #######################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_12 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 13:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_13 version
    ########################################################################################################
    #### CORE GlcNAc CAN MOVE IN THIS PROTOCOL BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    ########################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_13 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 14:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_14 version
    ###############################################
    #### CORE GlcNAc CAN MOVE IN THIS PROTOCOL ####
    ###############################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, False )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_14 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 15:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_15 version
    ####################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN ####
    ####################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_15 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 16:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_16 version
    ##################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN AND CORE GlcNAc CAN MOVE BUT STARTS AT NATIVE ####
    ##################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_16 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 17:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_17 version
    #############################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    #############################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_17 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 18:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_18 version
    #############################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN AND CORE GlcNAc CAN MOVE ####
    #############################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_18 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 19:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_19 version
    #############################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    #############################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_19 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.spin_using_ideal_omegas = True  # meaning use only 180, 60, -60. Don't sample within +/- something of these values
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 20:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_20 version
    #############################################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SPUN AROUND INSTEAD OF USING LCM DATA ####
    #############################################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_20 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = True  # at time of run, omega1 and omega2 can independently be 180, 60, or -60
    GlycanModelProtocol.spin_using_ideal_omegas = False  # meaning use 180, 60, -60 but sample within +/- 0-15 of these values
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = False
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 21:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_21 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_21 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 3,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = False
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 22:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_22 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_22 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 23:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_23 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44 } )

    # Protocol_23 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    # testing if the constraint is necessary
    GlycanModelProtocol.constraint_file = None
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 24:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_24 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    # fa_sol is 10% of normal value
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0, "fa_sol" : 0.09375 } )

    # Protocol_24 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 25:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_25 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_25 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 300
    GlycanModelProtocol.moves_per_trial = 3
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )

elif input_args.protocol_num == 26:
    # create the necessary minimization (and overall movement) MoveMap for Protocol_26 version
    ##########################################################################################
    #### THIS PROTOCOL INVOLVES PACKING AND CHI MIN BUT OMEGA 1 & 2 IS SET TO IgG FC DATA ####
    ##########################################################################################
    mm = MoveMap()
    for res_num in native_Fc_glycan_nums:
        mm.set_bb( res_num, True )
        mm.set_chi( res_num, True )
        if native_pose.residue( res_num ).is_branch_point():
            mm.set_branches( res_num, True )

    # create the desired scorefxn
    sf = get_fa_scorefxn_with_given_weights( { "fa_intra_rep" : 0.44, "atom_pair_constraint" : 1.0 } )

    # Protocol_26 creation and argument setting
    GlycanModelProtocol = Model3ay4Glycan( mm_in = mm, 
                                           sf_in = sf, 
                                           angle_max = 6.0 * 5,  # 6.0 comes from default angle_max from SmallMover and ShearMover
                                           dump_dir = input_args.structure_dir, 
                                           pmm = pmm )
    GlycanModelProtocol.trials = 200
    # what happens if we just make one move per trial?
    GlycanModelProtocol.moves_per_trial = 1
    GlycanModelProtocol.LCM_reset = True
    GlycanModelProtocol.use_population_ideal_LCM_reset = False
    GlycanModelProtocol.spin_carb_connected_to_prot = False
    GlycanModelProtocol.spin_using_ideal_omegas = False
    GlycanModelProtocol.set_native_omega = False
    GlycanModelProtocol.set_native_core = False
    GlycanModelProtocol.set_native_core_omegas_to_stats = True
    GlycanModelProtocol.ramp_sf = True
    GlycanModelProtocol.ramp_angle_max = True
    GlycanModelProtocol.angle_min = 6.0 * 3
    GlycanModelProtocol.fa_atr_ramp_factor = 2.0
    GlycanModelProtocol.fa_rep_ramp_factor = 0.01
    GlycanModelProtocol.minimize_each_round = True
    GlycanModelProtocol.pack_after_x_rounds = 3
    GlycanModelProtocol.make_small_moves = True
    GlycanModelProtocol.make_shear_moves = False
    GlycanModelProtocol.constraint_file = "project_constraint_files/native_3ay4_Gal_6A_1A_tol.cst"
    GlycanModelProtocol.verbose = input_args.verbose

    # write information to file (also prints to screen)
    GlycanModelProtocol.write_protocol_info_file( native_pose, input_args.protocol_num )
