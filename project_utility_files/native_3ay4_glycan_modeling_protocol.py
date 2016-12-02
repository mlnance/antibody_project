#!/usr/bin/python
__author__="morganlnance"


'''
The intent for this file is to hold all of the single-call functions for each protocol that is desired to be tested in modeling the 3ay4 glycan. For this to work appropriately, a file or two should hold all of the necessary "worker" functions, and this file should combine those functions into step to create a protocol. Then another script will have all of the needed instantiation information (making dirs, checking files, making the JobDistributor, etc), and will call a single protocol function from here in the JobDistributor
'''


# global imports
import sys, os



class Model3ay4Glycan:
    #def __init__( self, mm_in, sf_in, argument_file = None, pmm = None ):
    #:param argument_file: str( /path/to/file that contains needed setter arguments for this protocol class
    def __init__( self, mm_in, sf_in, angle_max, dump_dir, pmm = None ):
        '''
        No argument_file passed means this will construct the base class protocol
        :param mm_in: MoveMap
        :param sf_in: ScoreFunction
        :param angle_max: int( or float( what is the largest range you would like the torsions to be able to sample within? This means the biggest, single most possible is +/- angle_max/2 ) ) Default should be 6 for now
        :param dump_dir: str( /path/to/structure dump directory )
        :param pmm: PyMOL_Mover
        '''
        # pandas isn't on Jazz, so use the csv module if you can't use pandas
        # this is for watching energy per trial number
        try:
            import pandas as pd
            self.df = pd.DataFrame()
            self.pandas = True
        except:
            self.df = None
            self.pandas = False
            pass

        # default arguments
        self.name = "Model3ay4Glycan"

        # required arguments
        self.mm = mm_in
        self.sf = sf_in
        self.angle_max = angle_max
        self.dump_dir = dump_dir
        self.watch_sf = sf_in.clone()  # used when looking for convergence and verbose print outs
        self.pmm = pmm
        if self.pmm is not None:
            pmm.keep_history( True )

        # default arguments
        #self.glyco_file = None  # not used for native protocol right now
        self.outer_trials = 2
        self.inner_trials = 50
        self.moves_per_trial = 1
        #self.light_reset = False  # not used at the moment
        self.random_reset = False
        self.LCM_reset = True
        #self.LCM_main_ideal_reset = False  # not a bad idea, just need to make sure I compiled the data correctly
        self.use_population_ideal_LCM_reset = False
        self.set_native_omega = False
        self.set_native_core = False
        self.set_native_core_omegas_to_stats = False
        self.spin_carb_connected_to_prot = False
        self.spin_using_ideal_omegas = True
        self.ramp_sf = True
        self.ramp_angle_max = False
        self.angle_min = 6.0
        self.fa_atr_ramp_factor = 2.0
        self.fa_rep_ramp_factor = 0.5
        self.minimize_each_round = True
        self.pack_after_x_rounds = 0
        self.make_small_moves = True
        self.make_shear_moves = False
        self.move_all_torsions = True
        self.moveable_residues = []  # gets filled when .apply() is called
        self.constraint_file = None
        self.kT = 0.8
        self.mc = None  # MonteCarlo object
        self.mc_acceptance = None
        self.native_pose = None
        self.reset_pose = None
        self.pmm_name = None
        self.decoy_name = None

        # dumping arguments
        self.dump_reset_pose = False
        self.zip_dump_poses = False

        # watch and listen arguments
        self.verbose = False
        self.watch_E_vs_trial = False
        self.lowest_score_seen = None
        self.make_movie = False
        self.movie_poses = []
        self.movie_poses_dir = None


        '''
        ##############
        #### TODO ####   at some point make this be able to take command line arguments. Or the modeling.py program takes them and sets them
        ##############
        # read the arguments from the argument_file
        if argument_file is not None:
            try:
                with open( argument_file, "rb" ) as fh:
                    arg_lines = fh.readlines()
            except:
                print "\nThere is something wrong with your argument_file. I tried reading it in the %s class and had trouble.\n" %self.name
                sys.exit()

            # set the arguments
            for arg_line in arg_lines:
                arg = arg_line.split( '=' )[0].strip()
                if arg == "native_pdb_file":
                    self.native_pdb_file = arg
                elif arg == "structure_dir":
                    self.dump_dir = arg
                elif arg == "make_small_moves":
                    if arg.lower() == "true":
                        self.make_small_moves = True
                    elif arg.lower() == "false":
                        self.make_small_moves = False
                    else:
                        print "\nI don't recognize your argument for make_small_moves in your argument_file. It should be True or False.\n"
                        sys.exit()
                elif arg == "make_shear_moves":
                    if arg.lower() == "true":
                        self.make_shear_moves = True
                    elif arg.lower() == "false":
                        self.make_shear_moves = False
                    else:
                        print "\nI don't recognize your argument for make_shear_moves in your argument_file. It should be True or False.\n"
                        sys.exit()
                elif arg == "trials":
                    try:
                        self.inner_trials = int( trials )
                    except ValueError:
                        print "\nDid you give me a number for your trials argument in your argument_file?\n"
                        sys.exit()
                elif arg == "moves_per_trial":
                    try:
                        self.moves_per_trial = int( moves_per_trial )
                    except ValueError:
                        print "\nDid you give me a number for your moves_per_trial argument in your argument_file?\n"
                        sys.exit()
                elif arg == "nstruct":
                    try:
                        self.nstruct = int( nstruct )
                    except ValueError:
                        print "\nDid you give me a number for your nstruct argument in your argument_file?\n"
                        sys.exit()
                elif arg == "LCM_reset":
                    if arg.lower() == "true":
                        self.LCM_reset = True
                    elif arg.lower() == "false":
                        self.LCM_reset = False
                    else:
                        print "\nI don't recognize your argument for LCM_reset in your argument_file. It should be True or False.\n"
                        sys.exit()
                elif arg == "use_population_ideal_LCM_reset":
                    if arg.lower() == "true":
                        self.use_population_ideal_LCM_reset = True
                    elif arg.lower() == "false":
                        self.use_population_ideal_LCM_reset = False
                    else:
                        print "\nI don't recognize your argument for use_population_ideal_LCM_reset in your argument_file. It should be True or False.\n"
                        sys.exit()
                elif arg == "set_native_omega":
                    if arg.lower() == "true":
                        self.set_native_omega = True
                    elif arg.lower() == "false":
                        self.set_native_omega = False
                    else:
                        print "\nI don't recognize your argument for set_native_omega in your argument_file. It should be True or False.\n"
                        sys.exit()
                elif arg == "ramp_sf":
                    if arg.lower() == "true":
                        self.ramp_sf = True
                    elif arg.lower() == "false":
                        self.ramp_sf = False
                    else:
                        print "\nI don't recognize your argument for ramp_sf in your argument_file. It should be True or False.\n"
                        sys.exit()
                elif arg == "constraint_file":
                    self.constraint_file = arg
                elif arg == "protocol_num":
                    try:
                        self.protocol_num = int( arg )
                    except ValueError:
                        print "\nI don't recognize your argument for protocol_num in your argument_file. Did you give me a number?\n"
                        sys.exit()
                elif arg == "verbose":
                    if arg.lower() == "true":
                        self.verbose = True
                    elif arg.lower() == "false":
                        self.verbose = False
                    else:
                        print "\nI don't recognize your argument for verbose in your argument_file. It should be True or False.\n"
                        sys.exit()
                else:
                    print "\nI don't recognize this argument ( %s ). Did you add this in? Did you forget to add it into the %s class?" %( arg, self.name )
        '''


        ###############################
        #### SETUP THE DIRECTORIES ####
        ###############################
        # add the trailing '/' if needed
        if not self.dump_dir.endswith( '/' ):
            self.dump_dir += '/'

        # base_structs and lowest_E_structs
        self.base_structs_dir = self.dump_dir + "base_structs/"
        self.lowest_E_structs_dir = self.dump_dir + "lowest_E_structs/"

        # make the needed dump directories
        if not os.path.isdir( self.base_structs_dir ):
            try:
                os.mkdir( self.base_structs_dir )
            except:
                pass
        if not os.path.isdir( self.lowest_E_structs_dir ):
            try:
                os.mkdir( self.lowest_E_structs_dir )
            except:
                pass

        # create a working directory to be used for metric calculations
        self.metrics_dump_dir = self.dump_dir + "metrics_dump_dir"
        try:
            os.mkdir( self.metrics_dump_dir )
        except:
            pass


    def write_protocol_info_file( self, pose, protocol_num ):
        info_file_details = []
        info_file_details.append( "Protocol Number Used:\t\t\t%s\n" %protocol_num )
        info_file_details.append( "Native PDB filename:\t\t\t%s\n" %pose.pdb_info().name().split( '/' )[-1] )
        #info_file_details.append( "Sugar filename:\t\t\t\t%s\n" %input_args.glyco_file.split( '/' )[-1] )
        info_file_details.append( "Number of Protocol outer trials:\t%s\n" %self.outer_trials )
        info_file_details.append( "Number of Protocol inner trials:\t%s\n" %self.inner_trials )
        info_file_details.append( "Number of SugarSmallMoves per trial:\t%s\n" %self.moves_per_trial )
        if self.ramp_angle_max:
            angle_max_txt = "%s (meaning a max of +/- %s to either side of current which ramps down to %s throughout protocol)" %( self.angle_max, self.angle_max/2, self.angle_min/2 )
        else:
            angle_max_txt = "%s (meaning a max of +/- %s to either side of current)" %( self.angle_max, self.angle_max/2 )
        info_file_details.append( "Angle max (arc range available):\t%s\n" %angle_max_txt )
        #info_file_details.append( "Light reset of Fc glycan?:\t\t%s\n" %self.light_reset )
        info_file_details.append( "Random reset of Fc glycan?:\t\t%s\n" %self.random_reset )
        info_file_details.append( "LCM reset of Fc glycan?:\t\t%s\n" %self.LCM_reset )
        info_file_details.append( "Use population ideals in LCM reset?:\t%s\n" %self.use_population_ideal_LCM_reset )
        #info_file_details.append( "Use main ideal in LCM reset?:\t\t%s\n" %self.LCM_main_ideal_reset )
        info_file_details.append( "Move all torsions on a residue?:\t%s\n" %self.move_all_torsions )
        info_file_details.append( "Reset omega torsion back to native?:\t%s\n" %self.set_native_omega )
        info_file_details.append( "Reset core GlcNAc torsions to native?:\t%s\n" %self.set_native_core )
        info_file_details.append( "Reset core GlcNAc omegas to Fc data?:\t%s\n" %self.set_native_core_omegas_to_stats )
        info_file_details.append( "Spin carb connected to the protein?:\t%s\n" %self.spin_carb_connected_to_prot )
        info_file_details.append( "Spin carb using ideal omegas?:\t\t%s\n" %self.spin_using_ideal_omegas )
        info_file_details.append( "Using score ramping?:\t\t\t%s\n" %self.ramp_sf )
        info_file_details.append( "Using angle_max ramping?:\t\t%s\n" %self.ramp_angle_max )
        info_file_details.append( "Minimize after each move?:\t\t%s\n" %self.minimize_each_round )
        if self.pack_after_x_rounds is not 0:
            info_file_details.append( "Local pack after X rounds:\t\tX = %x\n" %self.pack_after_x_rounds )
        info_file_details.append( "Native constraint file used?:\t\t%s\n" %self.constraint_file )
        info_file_details.append( "Main structure directory:\t\t%s\n" %self.dump_dir )
        info_file_details.append( "Base structure directory:\t\t%s\n" %self.base_structs_dir )
        info_file_details.append( "Lowest E structure directory:\t\t%s\n" %self.lowest_E_structs_dir )
        # ScoreFunction information
        score_types = self.sf.get_nonzero_weighted_scoretypes()
        info_file_details.append( "\nScore weights used in ScoreFunction:\n" )
        for score_type in score_types:
            if str( score_type ) == "fa_atr":
                info_file_details.append( "%s: %s * %s in ramp\n" %( str( score_type ), str( self.sf.get_weight( score_type ) ), self.fa_atr_ramp_factor ) )
            elif str( score_type ) == "fa_rep":
                info_file_details.append( "%s: %s * %s in ramp\n" %( str( score_type ), str( self.sf.get_weight( score_type ) ), self.fa_rep_ramp_factor ) )
            else:
                info_file_details.append( "%s: %s\n" %( str( score_type ), str( self.sf.get_weight( score_type ) ) ) )
        self.info_file = ''.join( info_file_details )
        print "\n", self.info_file, "\n"

        # write out the info file with the collected info from above
        info_filename = self.dump_dir + "protocol_%s_run.info" %protocol_num
        with open( info_filename, "wb" ) as fh:
            fh.write( "Info for this run of %s\n\n" %self.name )
            fh.write( self.info_file )


    def apply( self, pose ):
        '''
        This is an adjustable 3ay4-glycan modeling protocol (as of my PyRosetta-3's most recent abilities 9/26/16). Use all the adjustable arguments to figure out what the best protocol actually is
        :param pose: Pose
        '''
        #################
        #### IMPORTS ####
        #################
        from os import popen
        from random import choice
        from rosetta import MoveMap, MinMover, MonteCarlo, PyMOL_Mover
        from rosetta.core.scoring import fa_atr, fa_rep
        from native_3ay4_glycan_modeling_protocol_functions import native_3ay4_Fc_glycan_LCM_reset, \
            add_constraints_to_pose, get_ramp_score_weight, get_ramp_angle_max, \
            native_3ay4_Fc_glycan_random_reset, SugarSmallMover, \
            make_RotamerTrialsMover, get_res_nums_within_radius


        # for watching the energy per trial
        if self.watch_E_vs_trial:
            self.trial_nums = []
            self.energies = []
            self.lowest_seen_energies = []

        ##########################
        #### COPY INPUT POSES ####
        ##########################
        # get the working and native pose (for this particular script, they are the same thing)
        self.native_pose = pose.clone()
        working_pose = pose.clone()
        self.decoy_name = working_pose.pdb_info().name()


        #######################################
        #### PREPARE FOR MOVIE, IF DESIRED ####
        #######################################
        # create a movie_poses_dir, if desired
        if self.make_movie:
            # each apply should have an empty self.movie_poses list
            self.movie_poses = []
            # make the directory, if needed
            self.movie_poses_dir = self.dump_dir + "movie_poses_dir/"
            if not os.path.isdir( self.movie_poses_dir ):
                try:
                    os.mkdir( self.movie_poses_dir )
                except:
                    pass


        ###############################
        #### GET MOVEABLE RESIDUES ####
        ###############################
        # get the moveable residues from passed MoveMap
        # moveable means the BackBone in the MoveMap was set to True and it is a carbohydrate
        self.moveable_residues = [ res_num for res_num in range( 1, working_pose.n_residue() + 1 ) if self.mm.get_bb( res_num ) and working_pose.residue( res_num ).is_carbohydrate() ]


        ###################
        #### LCM RESET ####
        ###################
        # reset the glycan using the data from the default.table used for the LinkageConformerMover
        if self.LCM_reset:
            working_pose.assign( native_3ay4_Fc_glycan_LCM_reset( mm = self.mm, 
                                                                  input_pose = working_pose, 
                                                                  use_population_ideal_LCM_reset = self.use_population_ideal_LCM_reset ) )
            if self.verbose:
                print "score of LCM reset:", self.watch_sf( working_pose )


        ##########################
        #### OR, RANDOM RESET ####
        ##########################
        # reset the glycan to random values from -180 to 180
        if self.random_reset:
            working_pose.assign( native_3ay4_Fc_glycan_random_reset( mm = self.mm, 
                                                                     input_pose = working_pose ) )
            if self.verbose:
                print "score of random reset:", self.watch_sf( working_pose )


        ######################################################
        #### SPIN CARB OF PROTEIN-CARBOHYDRATE CONNECTION ####
        ######################################################
        # the intent of this spin is to set the carbohydrate involved in the protein-carbohydrate connection into
        # a reasonable starting position after the reset. This is because current LCM data found in the default.table
        # was collected for surface glycans. These data are not reflective of the glycans found in the Ig system
        # use the MoveMap to see which residues have a parent connection to a protein
        if self.spin_carb_connected_to_prot:
            from native_3ay4_glycan_modeling_protocol_functions import spin_carbs_connected_to_prot

            # the spin takes residue 216 and 440 of 3ay4 and assigns either 180, 60, or -60 to omega1 and omega2 individually
            # can sample within +/- 15 of those three values as well, if specified
            working_pose.assign( spin_carbs_connected_to_prot( self.mm, 
                                                               input_pose = working_pose, 
                                                               spin_using_ideal_omegas = self.spin_using_ideal_omegas) )
            if self.verbose:
                print "score of core GlcNAc spin:", self.watch_sf( working_pose )


        #########################################
        #### HARDCODED OMEGA RESET TO NATIVE ####
        #########################################
        # reset the native omega torsion, if desired
        # hardcoded for now as this is not something that would be done in a real protocol
        if self.set_native_omega:
            working_pose.set_omega( 221, self.native_pose.omega( 221 ) )
            working_pose.set_omega( 445, self.native_pose.omega( 445 ) )
            if self.verbose:
                print "score of omega branch reset:", self.watch_sf( working_pose )


        ########################################
        #### HARDCODED CORE RESET TO NATIVE ####
        ########################################
        # reset the native torsions for the core GlcNAcs, if desired
        # hardcoded for now as this is not something that would be done in a real protocol
        if self.set_native_core:
            from rosetta.core.id import omega2_dihedral
            from rosetta.core.pose.carbohydrates import get_glycosidic_torsion, set_glycosidic_torsion

            # reset the phi, psi, omega, and omega2 torsions of residues 216 and 440 (core GlcNAc to ASN-69)
            # 216
            working_pose.set_phi( 216, self.native_pose.phi( 216 ) )
            working_pose.set_psi( 216, self.native_pose.psi( 216 ) )
            working_pose.set_omega( 216, self.native_pose.omega( 216 ) )
            set_glycosidic_torsion( omega2_dihedral, working_pose, 216, 
                                    get_glycosidic_torsion( omega2_dihedral, self.native_pose, 216 ) )
            # 440
            working_pose.set_phi( 440, self.native_pose.phi( 440 ) )
            working_pose.set_psi( 440, self.native_pose.psi( 440 ) )
            working_pose.set_omega( 440, self.native_pose.omega( 440 ) )
            set_glycosidic_torsion( omega2_dihedral, working_pose, 440, 
                                    get_glycosidic_torsion( omega2_dihedral, self.native_pose, 440 ) )
            if self.verbose:
                print "score of core GlcNAc reset:", self.watch_sf( working_pose )


        # no longer needed because I changed the default.table to store these data
        # thus, I am letting the LCM choose from IgG1 Fc stats to reset GlcNAc core omega1 and omega2
        '''
        ####################################################
        #### HARDCODED CORE RESET TO IgG1 Fc STATISTICS ####
        ####################################################
        # reset the torsions for the core GlcNAcs to values seen in IgG Fcs, if desired
        # hardcoded for now as this is not something that would be done in a real protocol
        if self.set_native_core_omegas_to_stats:
            from rosetta.core.id import omega_dihedral, omega2_dihedral
            from rosetta.core.pose.carbohydrates import set_glycosidic_torsion
            from rosetta.numeric.random import gaussian
            from native_3ay4_glycan_modeling_protocol_functions import calc_mean_degrees, calc_stddev_degrees

            # compile the data I have found from native crystal structures
            # subtracting 360 from positive numbers because...I think that's how I should do this
            phi_data = []
            psi_data = []
            omega1_data = [ -154.566, -162.504,  # 3ay4 - IgG1 Fc G2 to FcgRIIIa
                             -163.850, 176.923,  # 3ave - IgG1 Fc no paper, but it comes out as G0F2
                             -164.387, -146.038, # 5d4q - IgG1 Fc G2F1 ( check )
                             -146.449, -146.340, # 5d6d - IgG1 Fc G2F2 to FcgRIIIa ( check )
                             -166.996, -171.113  # 1h3x - IgG1 Fc G0F2
                             ]
            omega1_data = [ deg - 360 if deg > 0 else deg for deg in omega1_data ]
            omega2_data = [ 59.198, 59.055,  # 3ay4
                            53.823, 47.082,  # 3ave
                            48.590, 63.976,  # 5d4q
                            64.005, 63.988,  # 5d6d
                            45.997, 59.196   # 1h3x
                            ]
            omega2_data = [ deg - 360 if deg > 0 else deg for deg in omega2_data ]

            # pick a number within a gaussian distribution (I think??) of the torsions
            # mean + ( stddev * gaussian() )
            # not doing it the SmallMover way because I would want to be able to sample out more than 1 standard deviation
            # SmallMover way would be periodic_range( mean_torsion - torsion_stddev * rg().uniform(), 360.0 )
            # mean
            base_omega1 = calc_mean_degrees( omega1_data )
            base_omega2 = calc_mean_degrees( omega2_data )
            # standard deviation
            omega1_stddev = calc_stddev_degrees( omega1_data )
            omega2_stddev = calc_stddev_degrees( omega2_data )

            # 216
            new_omega1_216 = base_omega1 + ( omega1_stddev * gaussian() )
            new_omega2_216 = base_omega2 + ( omega2_stddev * gaussian() )
            # 440
            new_omega1_440 = base_omega1 + ( omega1_stddev * gaussian() )
            new_omega2_440 = base_omega2 + ( omega2_stddev * gaussian() )

            # reset the omega1 and omega2 torsions of residues 216 and 440 (core GlcNAc to ASN-69)
            # 216
            set_glycosidic_torsion( omega_dihedral, working_pose, 216, new_omega1_216 )
            set_glycosidic_torsion( omega2_dihedral, working_pose, 216, new_omega2_216 )
            # 440
            set_glycosidic_torsion( omega_dihedral, working_pose, 440, new_omega1_440 )
            set_glycosidic_torsion( omega2_dihedral, working_pose, 440, new_omega2_440 )

            if self.verbose:
                print "score of core GlcNAc reset using IgG Fc statistics:", self.watch_sf( working_pose )
        '''


        ############################################
        #### VISUALIZE AND STORE THE RESET POSE ####
        ############################################
        # visualize the pose that has had all components of reset completed
        try:
            if self.pmm_name is not None:
                working_pose.pdb_info().name( self.pmm_name )
                self.pmm.apply( working_pose )
                working_pose.pdb_info().name( self.decoy_name )
            else:
                self.pmm.apply( working_pose )
        except:
            pass

        # store a copy of the reset pose object
        # storing it here so that all types of resets could have been done, including the reset back to natives
        self.reset_pose = working_pose.clone()

        # dump the reset pose, if desired
        if self.dump_reset_pose:
            reset_name = self.reset_pose.pdb_info().name().split( ".pdb" )[0] + "_reset.pdb"
            self.reset_pose.dump_file( reset_name )
            # zip the dump file, if desired
            if self.zip_dump_poses:
                os.popen( "gzip %s" %reset_name )

        # add the reset pose to the list of poses for the movie, if desired
        if self.make_movie:
            # change the name to just "protocol_X_decoy_Y_0.pdb" without path location
            # X is protocol number, Y is decoy number, 0 for the movie number
            # state 0 is the reset_pose for when making a movie
            orig_name = self.reset_pose.pdb_info().name()
            reset_name = self.reset_pose.pdb_info().name().split( '/' )[-1].split( ".pdb" )[0] + "_0.pdb"
            self.reset_pose.pdb_info().name( reset_name )
            self.movie_poses.append( self.reset_pose.clone() )
            self.reset_pose.pdb_info().name( orig_name )

        # append all the information to the lists
        if self.watch_E_vs_trial:
            # 0th trial is the reset pose
            self.trial_nums.append( 0 )
            self.energies.append( self.watch_sf( working_pose ) )
            self.lowest_seen_energies.append( self.watch_sf( working_pose ) )


        #########################
        #### ADD CONSTRAINTS ####
        #########################
        if self.constraint_file is not None:
            from rosetta.core.scoring import atom_pair_constraint
            try:
                # set the constraints from the constraint_file
                working_pose.assign( add_constraints_to_pose( self.constraint_file, working_pose ) )
                # add atom_pair_constraint to the ScoreFunction, if needed
                self.sf.set_weight_if_zero( atom_pair_constraint, 1.0 )
            except:
                print "\nThere was something wrong with your constraint file. Are you sure it exists? Are you sure you used the correct names for the constraints? Check on %s\n" %self.constraint_file
                sys.exit()


        ####################################
        #### MAKE THE MONTECARLO OBJECT ####
        ####################################
        # create the MonteCarlo object
        if self.mc is None:
            self.mc = MonteCarlo( working_pose, self.sf, self.kT )
        # reset everything just to be safe since MonteCarlo is a mess
        self.mc.reset_counters()
        self.mc.clear_poses()
        self.mc.reset( working_pose )
        self.mc.reset_scorefxn( working_pose, self.sf )


        #############################
        #### PREPARE FOR RAMPING ####
        #############################
        # store the original fa_atr and fa_rep
        # this will be accessed as the outer_trials pass and the sf gets re-ramped
        FA_ATR_ORIG = self.sf.get_weight( fa_atr )
        FA_REP_ORIG = self.sf.get_weight( fa_rep )
        # set the local angle_max to the passed angle_max value
        # this will be adjusted if ramp_angle_max is True
        angle_max = self.angle_max


        ########################################
        #### MAKE AND APPLY THE SUGAR MOVER ####
        ########################################
        # for as many trials as specified
        num_mc_accepts = 0
        num_mc_checks = 0
        mc_acceptance = -1  # -1 so that it will play nice in the metrics file if the number doesn't get updated somehow
        movie_num = 1
        # for each set of outer trials
        for outer_trial in range( 1, self.outer_trials + 1 ):
            # reset the MonteCarlo object with the working_pose before each chunk of inner_trials
            self.mc.reset_counters()
            self.mc.clear_poses()
            self.mc.reset( working_pose )
            # whether ramp_sf is True or False, this should reset the MC to the latest pose with the passed sf weights
            self.mc.reset_scorefxn( working_pose, self.sf )

            # inner_trial must be 1 to N because decoy 0 will be the reset_pose when creating a movie
            for inner_trial in range( 1, self.inner_trials + 1 ):
                # print current score
                if self.verbose:
                    print "\nstarting score:", self.watch_sf( working_pose )

            	####################################
            	#### RAMP WITH ANGLE MAX UPDATE ####
            	####################################
                if self.ramp_angle_max:
                    # the last ~10% of moves in the inner_trials will be the angle_min
                    # the first move should use the passed angle_max
                    if inner_trial == 1:
                        angle_max = self.angle_max
                    # otherwise, adjust the angle_max according to which inner_trial number we are on
                    else:
                        # by the end of the protocol our angle_max should be angle_min
                        # angle_min is 6.0 by default which is the angle_max for loop residues in the BackboneMover
                        angle_max = get_ramp_angle_max( current_angle_max = angle_max, 
                                                        target_angle_max = self.angle_min, 
                                                        current_step = inner_trial, 
                                                        total_steps = self.inner_trials )


            	#########################################
            	#### RAMP WITH SCORE FUNCTION UPDATE ####
            	#########################################
                if self.ramp_sf:
                    # if this is the first move, adjust the fa_atr and fa_rep terms by the corresponding factors
                    # after the first move, we will ramp the score up or down accordingly back to the original value
                    # the last ~10% of moves in the inner_trials will be the original values
                    if inner_trial == 1:
                        # adjust the fa_atr weight by the passed fa_atr_ramp_factor
                        FA_ATR_NEW = FA_ATR_ORIG * self.fa_atr_ramp_factor
                        self.sf.set_weight( fa_atr, FA_ATR_NEW )
                        # adjust the fa_rep weight by the passed fa_rep_ramp_factor
                        FA_REP_NEW = FA_REP_ORIG * self.fa_rep_ramp_factor
                        self.sf.set_weight( fa_rep, FA_REP_NEW )
                    # else, adjust the score weight based on the current inner_trial step
                    else:
                        # ramp up or down the appropriate scoring terms
                        self.sf.set_weight( fa_atr, get_ramp_score_weight( current_weight = self.sf.get_weight( fa_atr ), 
                                                                           target_weight = FA_ATR_ORIG, 
                                                                           current_step = inner_trial, 
                                                                           total_steps = self.inner_trials ) )
                        self.sf.set_weight( fa_rep, get_ramp_score_weight( current_weight = self.sf.get_weight( fa_rep ), 
                                                                           target_weight = FA_REP_ORIG, 
                                                                           current_step = inner_trial, 
                                                                           total_steps = self.inner_trials ) )
                    # DEBUG
                    #print "\n".join( [ "%s %s" %( str( score_type ), str( self.sf.get_weight( score_type ) ) ) for score_type in self.sf.get_nonzero_weighted_scoretypes() ] )
                    # give ramped sf back to MonteCarlo object
                    # this is because PyRosetta apparently doesn't do the whole pointer thing with the sf
                    self.mc.score_function( self.sf )


            	##########################
            	#### MAKE SUGAR MOVES ####
            	##########################
                # make as many moves per trial as desired
                # SugarSmallMover
                if self.make_small_moves:
                    working_pose.assign( SugarSmallMover( self.mm, self.moves_per_trial, angle_max, working_pose, 
                                                          move_all_torsions = self.move_all_torsions ) )
                # SugarShearMover
                elif self.make_shear_moves:
                    pass

                # relay score information
                if self.verbose:
                    print "score after sugar moves:", self.watch_sf( working_pose )


                ##################################
                #### PREPARE FOR PACK AND MIN ####
                ##################################
                # if you aren't going to pack this round, then the pack_and_min_residues will mirror the self.moveable_residues
                # if you are packing this round, then you need residues surrounding the self.moveable_residues as well
                pack_and_min_residues = [ res_num for res_num in self.moveable_residues ]
                # if you aren't going to pack this round, then the minimizer's MoveMap should only include the self.moveable_residues (bb and chi)
                # if you are packing this round, then the MoveMap has bb and chi for self.moveable_residues and chi for surrounding_residues
                # default MoveMap creation has everything set to False, so set the appropriate residues to True
                min_mm = MoveMap()
                for moveable_res in pack_and_min_residues:
                    min_mm.set_bb( moveable_res, True )
                    min_mm.set_chi( moveable_res, True )
                # check if we are packing, adjust the pack_and_min_residues and min_mm accordingly
                if self.pack_after_x_rounds > 0:
                    if inner_trial % self.pack_after_x_rounds == 0:
                        # this function I wrote uses the nbr_atom to calculate distances, which is about the C4 atom
                        # 10A should be enough to include Tyr296 if the fucose is in the pose
                        surrounding_residues = get_res_nums_within_radius( self.moveable_residues, working_pose, 
                                                                           radius = 10, 
                                                                           include_passed_res_nums = False )
                        # add the surrounding_residues to the list of residues that will be packed and minimized
                        pack_and_min_residues.extend( surrounding_residues )
                        # set chi to True in the min_mm
                        for surrounding_res in surrounding_residues:
                            min_mm.set_chi( surrounding_res, True )


            	##############
            	#### PACK ####
            	##############
                # pack the sugars and surrounding residues, if desired
                if self.pack_after_x_rounds > 0:
                    if inner_trial % self.pack_after_x_rounds == 0:
                        # have to make the packer task each time because the residues surrounding the sugars
                        # will likely change after each move, thus need to make a new RotamerTrialsMover
                        # pack_radius of None means that only the residues specified will be packed
                        rotamer_trials_mover = make_RotamerTrialsMover( pack_and_min_residues, self.sf, working_pose, 
                                                                        pack_radius = None )
                        rotamer_trials_mover.apply( working_pose )
                        if self.verbose:
                            print "score after local pack:", self.watch_sf( working_pose )


            	###################
                #### MINIMIZE #####
            	###################
                # make the MinMover from the pack_and_min_residues, if desired
                # if a pack was not just done, this will only include the self.moveable_residues (bb and chi)
                # if a pack was just done, this will include self.moveable_residues (bb and chi) and the surrounding_residues (chi)
                if self.minimize_each_round:
                    # make and apply the min_mover
                    min_mover = MinMover( movemap_in = min_mm, 
                                          scorefxn_in = self.sf,
                                          min_type_in = "dfpmin_strong_wolfe",
                                          #min_type_in = "lbfgs_armijo_nonmonotone",
                                          tolerance_in = 0.01,
                                          use_nb_list_in = True )
                    min_mover.apply( working_pose )
                    if self.verbose:
                        print "score after min:", self.watch_sf( working_pose )


            	###############################
            	#### ACCEPT OR REJECT MOVE ####
            	###############################
                # accept or reject the total move using the MonteCarlo object
                if self.mc.boltzmann( working_pose ):
                    # add the accepted-move pose to the list of poses for the movie, if desired
                    if self.make_movie:
                        # change the name to just "protocol_X_decoy_Y_Z.pdb" without path location
                        # X is protocol number, Y is decoy number, Z is movie number
                        # ex name) protocol_10_decoy_5_23.pdb
                        working_pose.pdb_info().name( self.decoy_name.split( '/' )[-1].split( ".pdb" )[0] + "_%s.pdb" %movie_num )
                        self.movie_poses.append( working_pose.clone() )
                        working_pose.pdb_info().name( self.decoy_name )
                        movie_num += 1
                    # up the counters and send to pymol
                    num_mc_accepts += 1
                    try:
                        if self.pmm_name is not None:
                            working_pose.pdb_info().name( self.pmm_name )
                            self.pmm.apply( working_pose )
                            working_pose.pdb_info().name( self.decoy_name )
                        else:
                            self.pmm.apply( working_pose )
                    except:
                        pass
                num_mc_checks += 1

                # for watching energy during a protocol
                if self.watch_E_vs_trial:
                    # collect the lowest-scoring decoy seen using the watch_sf
                    # can't use the mc.lowest_score() because that uses the ramped sf, so the score is variable
                    # if this is the first trial done, update the lowest_score_seen to the value
                    # of the decoy after making the first move
                    if self.lowest_score_seen is None:
                        self.lowest_score_seen = self.watch_sf( working_pose )
                    # otherwise, update the lowest_score_seen if the current working_pose has a lower score
                    elif self.watch_sf( working_pose ) < self.lowest_score_seen:
                        self.lowest_score_seen = self.watch_sf( working_pose )
                    # append all the information to the lists
                    self.trial_nums.append( inner_trial + self.inner_trials * ( outer_trial - 1 ) )
                    self.energies.append( self.watch_sf( working_pose ) )
                    self.lowest_seen_energies.append( self.lowest_score_seen )

                # print out the MC acceptance rate every 3 trials and on the last trial
                mc_acceptance = round( ( float( num_mc_accepts ) / float( num_mc_checks ) * 100 ), 2 )
                if self.verbose:
                    if inner_trial % 3 == 0 or inner_trial == self.inner_trials:
                        print "Moves made so far:", num_mc_checks,
                        print "  Moves accepted:", num_mc_accepts,
                        print "  Acceptance rate:", mc_acceptance

        # finished outer_trials and inner_trials
        # add any relevant data to the class object
        self.mc_acceptance = mc_acceptance

        # for watching energy during a protocol
        if self.watch_E_vs_trial:
            E_vs_trial_filename = self.dump_dir + self.decoy_name.split( '/' )[-1].split( ".pdb" )[0] + "_E_vs_trial.csv"
            # use a pandas DataFrame, if possible
            if self.pandas:
                self.df[ "trial_nums" ] = self.trial_nums
                self.df[ "total_score" ] = self.energies
                self.df[ "lowest_score" ] = self.lowest_seen_energies
                self.df.to_csv( E_vs_trial_filename )
            # otherwise I'm on jazz and can't use a DataFrame, so use a csv file instead
            else:
                import csv
                data_rows = zip( self.trial_nums, self.energies, self.lowest_seen_energies )
                with open( E_vs_trial_filename, "wb" ) as fh:
                    csvwriter = csv.writer( fh )
                    csvwriter.writerow( [ "trial_num", "score", "lowest_score" ] )
                    [ csvwriter.writerow( row ) for row in data_rows ]

        return working_pose
