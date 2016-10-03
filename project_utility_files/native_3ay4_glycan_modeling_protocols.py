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
        # default arguments
        self.name = "Model3ay4Glycan"

        # required arguments
        self.mm = mm_in
        self.sf = sf_in
        self.angle_max = angle_max
        self.dump_dir = dump_dir
        self.watch_sf = sf_in.clone()  # used when looking for convergence and verbose print outs
        self.pmm = pmm

        # default arguments
        #self.glyco_file = None  # not used for native protocol right now
        self.trials = 50
        self.moves_per_trial = 1
        #self.light_reset = False  # not used at the moment
        self.random_reset = False
        self.LCM_reset = True
        #self.LCM_main_ideal_reset = False  # not a bad idea, just need to make sure I compiled the data correctly
        self.use_population_ideal_LCM_reset = False
        self.set_native_omega = True
        self.ramp_sf = True
        self.fa_atr_ramp_factor = 2.0
        self.fa_rep_ramp_factor = 0.5
        self.minimize_each_round = True
        self.make_small_moves = True
        self.make_shear_moves = False
        self.constraint_file = None
        self.kT = 0.8
        self.mc = None  # MonteCarlo object
        self.mc_acceptance = None
        self.min_mover = None
        self.reset_pose_obj = None

        # "optional" arguments
        self.verbose = False


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
                        self.trials = int( trials )
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
        info_file_details.append( "Native PDB filename:\t\t\t%s\n" %pose.pdb_info().name().split( '/' )[-1] )
        #info_file_details.append( "Sugar filename:\t\t\t\t%s\n" %input_args.glyco_file.split( '/' )[-1] )
        info_file_details.append( "Number of SugarSmallMove trials:\t%s\n" %self.trials )
        info_file_details.append( "Number of SugarSmallMoves per trial:\t%s\n" %self.moves_per_trial )
        angle_max_txt = "%s (meaning a max of +/- %s to either side of current)" %( self.angle_max, self.angle_max/2 )
        info_file_details.append( "Angle max (arc range available):\t%s\n" %angle_max_txt )
        #info_file_details.append( "Light reset of Fc glycan?:\t\t%s\n" %self.light_reset )
        info_file_details.append( "Random reset of Fc glycan?:\t\t%s\n" %self.random_reset )
        info_file_details.append( "LCM reset of Fc glycan?:\t\t%s\n" %self.LCM_reset )
        info_file_details.append( "Use population ideals in LCM reset?:\t%s\n" %self.use_population_ideal_LCM_reset )
        #info_file_details.append( "Use main ideal in LCM reset?:\t\t%s\n" %self.LCM_main_ideal_reset )
        info_file_details.append( "Reset omega torsion back to native?:\t%s\n" %self.set_native_omega )
        info_file_details.append( "Using score ramping?:\t\t\t%s\n" %self.ramp_sf )
        info_file_details.append( "Minimize after each move?:\t\t%s\n" %self.minimize_each_round )
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
        from random import choice
        from rosetta import MinMover, MonteCarlo
        from rosetta.core.scoring import fa_atr, fa_rep
        #from antibody_functions import native_Fc_glycan_nums_except_core_GlcNAc  # shouldn't need this as a MoveMap is passed to create this class
        from native_3ay4_glycan_modeling_protocol_functions import native_3ay4_Fc_glycan_LCM_reset, \
            add_constraints_to_pose, get_ramp_score_weight, native_3ay4_Fc_glycan_random_reset


        ########################################
        #### COPY POSES AND MAKE VISUALIZER ####
        ########################################
        # get the working and native pose (for this particular script, they are the same thing)
        native_pose = pose.clone()
        working_pose = pose.clone()


        ###################
        #### LCM RESET ####
        ###################
        # reset the glycan using the data from the default.table used for the LinkageConformerMover
        if self.LCM_reset:
            working_pose.assign( native_3ay4_Fc_glycan_LCM_reset( mm = self.mm, 
                                                                  input_pose = working_pose, 
                                                                  use_population_ideal_LCM_reset = self.use_population_ideal_LCM_reset ) )
            # store a copy of the reset pose object
            self.reset_pose_obj = working_pose.clone()

            # visualize and relay score information
            try:
                self.pmm.apply( working_pose )
            except:
                pass
            if self.verbose:
                print "score of LCM reset:", self.watch_sf( working_pose )


        ##########################
        #### OR, RANDOM RESET ####
        ##########################
        # reset the glycan to random values from -180 to 180
        if self.random_reset:
            working_pose.assign( native_3ay4_Fc_glycan_random_reset( mm = self.mm, 
                                                                     input_pose = working_pose ) )
            # store a copy of the reset pose object
            self.reset_pose_obj = working_pose.clone()

            # visualize and relay score information
            try:
                self.pmm.apply( working_pose )
            except:
                pass
            if self.verbose:
                print "score of random reset:", self.watch_sf( working_pose )


        #########################################
        #### HARDCODED OMEGA RESET TO NATIVE ####
        #########################################
        # reset the native omega torsion, if desired
        # hardcoded for now as this is not something that would be done in a real protocol
        if self.set_native_omega:
            working_pose.set_omega( 221, native_pose.omega( 221 ) )
            working_pose.set_omega( 445, native_pose.omega( 445 ) )


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
                print "\nThere was something wrong with your constraint file. Are you sure it exists? Are you sure you used the correct names for the constraints?\n"
                sys.exit()


        ############################
        #### MIN MOVER CREATION ####
        ############################
        # make the MinMover from the passed MoveMap, if needed
        if self.min_mover is None:
            self.min_mover = MinMover( movemap_in = self.mm, 
                                       scorefxn_in = self.sf,
                                       min_type_in = "dfpmin_strong_wolfe",
                                       tolerance_in = 0.01,
                                       use_nb_list_in = True )


        ####################################
        #### MAKE THE MONTECARLO OBJECT ####
        ####################################
        # create the MonteCarlo object, if needed
        if self.mc is None:
            self.mc = MonteCarlo( working_pose, self.sf, self.kT )


        ########################################
        #### MAKE AND APPLY THE SUGAR MOVER ####
        ########################################
        # imports 
        from native_3ay4_glycan_modeling_protocol_functions import SugarSmallMover

        # for as many trials as specified
        num_mc_accepts = 0
        num_mc_checks = 0
        for trial_num in range( 1, self.trials + 1 ):
            # print current score
            if self.verbose:
                print "\nstarting score:", self.watch_sf( working_pose )

            #########################################
            #### RAMP WITH SCORE FUNCTION UPDATE ####
            #########################################
            if self.ramp_sf:
                # if this is the first move, adjust the fa_atr and fa_rep terms by the corresponding factors
                if trial_num == 1:
                    # store the original fa_atr and fa_rep
                    FA_ATR_ORIG = self.sf.get_weight( fa_atr )
                    FA_REP_ORIG = self.sf.get_weight( fa_rep )
                    # adjust the fa_atr weight
                    FA_ATR_NEW = FA_ATR_ORIG * self.fa_atr_ramp_factor
                    self.sf.set_weight( fa_atr, FA_ATR_NEW )
                    # adjust the fa_rep weight
                    FA_REP_NEW = FA_REP_ORIG * self.fa_rep_ramp_factor
                    self.sf.set_weight( fa_rep, FA_REP_NEW )

                # else, adjust the score weight
                else:
                    # ramp up or down the appropriate scoring terms
                    self.sf.set_weight( fa_atr, get_ramp_score_weight( self.sf.get_weight( fa_atr ), 
                                                                       FA_ATR_ORIG, 
                                                                       trial_num, 
                                                                       self.trials ) )
                    self.sf.set_weight( fa_rep, get_ramp_score_weight( self.sf.get_weight( fa_rep ), 
                                                                       FA_REP_ORIG, 
                                                                       trial_num, 
                                                                       self.trials ) )
                # give ramped sf back to MC and MinMover
                # this is because PyRosetta apparently doesn't do the whole pointer thing with the sf
                self.mc.score_function( self.sf )
                self.min_mover.score_function( self.sf )


            ##########################
            #### MAKE SUGAR MOVES ####
            ##########################
            # make as many moves per trial as desired
            # SugarSmallMover
            if self.make_small_moves:
                working_pose.assign( SugarSmallMover( self.mm, self.moves_per_trial, self.angle_max, working_pose, 
                                                      set_phi = True, 
                                                      set_psi = True, 
                                                      set_omega = True ) )
            # SugarShearMover
            elif self.make_shear_moves:
                pass

            # relay score information
            if self.verbose:
                print "score after sugar moves:", self.watch_sf( working_pose )


            ###################
            #### MINIMIZE #####
            ###################
            # minimize the backbone of the sugars, if desired
            if self.minimize_each_round:
                self.min_mover.apply( working_pose )
                if self.verbose:
                    print "score after min:", self.watch_sf( working_pose )


            ###############################
            #### ACCEPT OR REJECT MOVE ####
            ###############################
            # accept or reject the total move using the MonteCarlo object
            if self.mc.boltzmann( working_pose ):
                # up the counters and send to pymol
                num_mc_accepts += 1
                try:
                    self.pmm.apply( working_pose )
                except:
                    pass
            num_mc_checks += 1

            # print out the MC acceptance rate every 3 trials and on the last trial
            mc_acceptance = round( ( float( num_mc_accepts ) / float( num_mc_checks ) * 100 ), 2 )
            if self.verbose:
                if trial_num % 3 == 0 or trial_num == self.trials:
                    print "Moves made so far:", num_mc_checks,
                    print "  Moves accepted:", num_mc_accepts,
                    print "  Acceptance rate:", mc_acceptance

        # add any relevant data to the class object
        self.mc_acceptance = mc_acceptance

        return working_pose
