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

# else I haven't made this protocol number yet
else:
    print "\nI haven't created the protocol number you gave me yet.\n"
    sys.exit()
