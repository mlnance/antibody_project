#!/usr/bin/python
__author__="morganlnance"


'''
The intent for this file is to hold all of the single-call functions for each protocol that is desired to be tested in modeling the 3ay4 glycan. For this to work appropriately, a file or two should hold all of the necessary "worker" functions, and this file should combine those functions into step to create a protocol. Then another script will have all of the needed instantiation information (making dirs, checking files, making the JobDistributor, etc), and will call a single protocol function from here in the JobDistributor

protocol_0: Base protocol. LCM reset starting from ideal for each linkage and sampling within +/- 1 stdev from there. SSM-200 am3 3mpt. am3 applies to each linkage and torsion equally. Core GlcNAc does not move. After LCM reset and before any movements, native omega at branch residues is restored. native_3ay4_Gal_5A_1A_tol.cst is used. Standard fa_intra_rep (0.44) sf. No packing at all, but minimization of just Fc glycan bb (no chi, no branch) before each MonteCarlo call

protocol_1: Comparison protocol. LCM reset starting from ideal for each linkage (no stdev sampling). SSM-200 am3 3mpt. am3 applies to each linkage and torsion equally. Core GlcNAc does not move. After LCM reset and before any movements, native omega at branch residues is restored. native_3ay4_Gal_5A_1A_tol.cst is used. Standard fa_intra_rep (0.44) sf. No packing at all, but minimization of just Fc glycan bb (no chi, no branch) before each MonteCarlo call
'''


# global imports
import sys, os



class ProtocolBase:
    def __init__( self, mm_in, sf_in, dump_dir, util_dir ):
        '''
        :param mm_in: MoveMap
        :param sf_in: ScoreFunction
        :param dump_dir: str( /path/to/dump location )
        :param util_dir: str( /path/to/utility files )
        '''
        # default arguments
        self.name = "Protocol"

        # required arguments
        self.mm = mm_in
        self.sf = sf_in
        self.watch_sf = sf_in.clone()  # used when looking for convergence and verbose print outs
        self.dump_dir = dump_dir
        self.util_dir = util_dir

        # protocol-specific information
        # setup up to emulate default Protocol_0
        self.glyco_file = None  # not used for native protocol right now
        self.trials = 10
        #self.trials = 200
        self.moves_per_trial = 3
        self.angle_multiplier = 3
        self.light_reset = False
        self.LCM_reset = True
        self.LCM_main_ideal_reset = False
        self.LCM_population_ideal_reset = False
        self.set_native_omega = True
        self.ramp_sf = True
        self.constraint_file = None
        #self.constraint_file = "project_constraint_files/native_3ay4_Gal_5A_1A_tol.cst"
        self.scorefxn_file = None
        self.kT = 0.8
        self.mc = None  # MonteCarlo object
        self.make_small_moves = True
        self.make_shear_moves = False

        # "optional" arguments
        self.verbose = False


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


    def write_protocol_info_file( self, pose ):
        info_file_details = []
        info_file_details.append( "Native PDB filename:\t\t\t%s\n" %pose.pdb_info().name().split( '/' )[-1] )
        #info_file_details.append( "Sugar filename:\t\t\t\t%s\n" %input_args.glyco_file.split( '/' )[-1] )
        info_file_details.append( "Number of SugarSmallMove trials:\t%s\n" %self.trials )
        info_file_details.append( "Number of SugarSmallMoves per trial:\t%s\n" %self.moves_per_trial )
        info_file_details.append( "Angle multiplier used:\t\t\t%s\n" %self.angle_multiplier )
        info_file_details.append( "Light reset of Fc glycan?:\t\t%s\n" %self.light_reset )
        info_file_details.append( "LCM reset of Fc glycan?:\t\t%s\n" %self.LCM_reset )
        info_file_details.append( "Use population ideals in LCM reset?:\t%s\n" %self.LCM_population_ideal_reset )
        info_file_details.append( "Use main ideal in LCM reset?:\t\t%s\n" %self.LCM_main_ideal_reset )
        info_file_details.append( "Reset omega torsion back to native?:\t%s\n" %self.set_native_omega )
        info_file_details.append( "Using score ramping?:\t\t\t%s\n" %self.ramp_sf )
        info_file_details.append( "Native constraint file used?:\t\t%s\n" %self.constraint_file )
        info_file_details.append( "ScoreFunction file used?:\t\t%s\n" %self.scorefxn_file )
        info_file_details.append( "Main structure directory:\t\t%s\n" %self.dump_dir )
        info_file_details.append( "Base structure directory:\t\t%s\n" %self.base_structs_dir )
        info_file_details.append( "Lowest E structure directory:\t\t%s\n" %self.lowest_E_structs_dir )
        # ScoreFunction information
        score_types = self.sf.get_nonzero_weighted_scoretypes()
        info_file_details.append( "\nScore weights used in main_sf:\n" )
        for score_type in score_types:
            if str( score_type ) == "fa_atr":
                info_file_details.append( "%s: %s * 2 in ramp\n" %( str( score_type ), str( self.sf.get_weight( score_type ) ) ) )
            elif str( score_type ) == "fa_rep":
                info_file_details.append( "%s: %s * 0.5 in ramp\n" %( str( score_type ), str( self.sf.get_weight( score_type ) ) ) )
            else:
                info_file_details.append( "%s: %s\n" %( str( score_type ), str( self.sf.get_weight( score_type ) ) ) )
        self.info_file = ''.join( info_file_details )
        print "\n", self.info_file, "\n"

        # write out the info file with the collected info from above
        info_filename = self.dump_dir + "protocol_0_run.info"
        with open( info_filename, "wb" ) as fh:
            fh.write( "Info for this run of %s\n\n" %self.name )
            fh.write( self.info_file )




####################
#### PROTOCOL_0 ####
####################

class Protocol_0( ProtocolBase ):
    '''
    Protocol_0: am3, 3mpt, native_3ay4_Gal_5A_1A_tol, SSM-200, LCM reset with +/- stdev, core GlcNAc unmoved, fa_intra_rep sf used (given from outside call), no packing, minimization of Fc glycan before each MC call
    '''
    def __init__( self, mm_in, sf_in, dump_dir, util_dir ):
        '''
        :param mm_in: MoveMap
        :param sf_in: ScoreFunction
        :param dump_dir: str( /path/to/dump location )
        :param util_dir: str( /path/to/utility files )
        '''
        # inherited arguments
        ProtocolBase.__init__( self, mm_in, sf_in, dump_dir, util_dir )
        # keep default arguments of base class Protocol for this Protocol_0, except add the constraint file (default is None)
        self.name = "Protocol_0"
        self.constraint_file = "project_constraint_files/native_3ay4_Gal_5A_1A_tol.cst"


    def apply( self, pose ):
        '''
        Base protocol. This is the most generic protocol (as of my PyRosetta-3's most recent abilities 9/26/16). Use to compare how other more "unique" protocols do
        :param pose: Pose
        '''
        #################
        #### IMPORTS ####
        #################
        from random import choice
        from rosetta import PyMOL_Mover, MinMover, MonteCarlo
        #from antibody_functions import native_Fc_glycan_nums_except_core_GlcNAc
        from native_3ay4_glycan_modeling_protocol_functions import native_3ay4_Fc_glycan_LCM_reset, \
            add_constraints_to_pose, ramp_score_weight


        ########################################
        #### COPY POSES AND MAKE VISUALIZER ####
        ########################################
        # get the working and native pose (for this particular script, they are the same thing)
        native_pose = pose.clone()
        working_pose = pose.clone()

        # make the PyMOL_Mover
        pmm = PyMOL_Mover()
        pmm.keep_history( True )


        ################################################
        #### GET RESIDUES TO WORK WITH FROM MOVEMAP ####
        ################################################
        # residues that are allowed to move are based on residues with BB set to True in the MoveMap and the residue is a carbohydrate
        # I'm copying this idea from GlycanRelaxMover, but overall it makes sense anyway. Only things with BB freedom should be sampled
        moveable_res_nums = [ res_num for res_num in range( 1, working_pose.n_residue() + 1 ) if self.mm.get_bb( res_num ) and working_pose.residue( res_num ).is_carbohydrate() ]


        ###################
        #### LCM RESET ####
        ###################
        working_pose.assign( native_3ay4_Fc_glycan_LCM_reset( input_pose = working_pose, 
                                                              residue_numbers = moveable_res_nums, 
                                                              use_population_ideal_LCM_reset = False ) )
        pmm.apply( working_pose )
        if self.verbose:
            print "score of LCM reset:", self.watch_sf( working_pose )

        # reset the native omega torsion, if desired
        # hardcoded for now as this is not something that would be done in a real protocol
        if self.set_native_omega:
            working_pose.set_omega( 221, native_pose.omega( 221 ) )
            working_pose.set_omega( 445, native_pose.omega( 445 ) )


        #########################
        #### ADD CONSTRAINTS ####
        #########################
        if self.constraint_file is not None:
            try:
                # set the constraints
                working_pose.assign( add_constraints_to_pose( self.constraint_file, working_pose ) )
                # add atom_pair_constraint to the ScoreFunction
                from rosetta import score_type_from_name
                self.sf.set_weight_if_zero( score_type_from_name( "atom_pair_constraint" ), 1.0 )
            except:
                print "\nThere was something wrong with your constraint file. Are you sure it exists? Are you sure you used the correct names for the constraints?\n"
                sys.exit()


        ############################
        #### MIN MOVER CREATION ####
        ############################
        # make the MinMover from the passed MoveMap
        Fc_glycan_min_mover = MinMover( movemap_in = self.mm, 
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


        #################################
        #### RAMP THE SCORE FUNCTION ####
        #################################
        # raise the fa_atr term and lower the fa_rep term in the ScoreFunction for ramping, if desired
        if self.ramp_sf:
            # store the original fa_atr and fa_rep
            FA_ATR_ORIG = self.sf.get_weight( score_type_from_name( "fa_atr" ) )
            FA_REP_ORIG = self.sf.get_weight( score_type_from_name( "fa_rep" ) )
            # double the fa_atr weight
            FA_ATR_NEW = FA_ATR_ORIG * 2
            self.sf.set_weight( score_type_from_name( "fa_atr" ), FA_ATR_NEW )
            # half the fa_rep weight
            FA_REP_NEW = FA_REP_ORIG * 0.5
            self.sf.set_weight( score_type_from_name( "fa_rep" ), FA_REP_NEW )


        ########################################
        #### MAKE AND APPLY THE SUGAR MOVER ####
        ########################################
        if self.make_small_moves:
            # imports 
            from rosetta import SmallMover
            from native_3ay4_glycan_modeling_protocol_functions import SugarSmallMover

            # make a SmallMover as to get the angle_max value
            sh = SmallMover()
            angle_max = sh.get_angle_max( 'L' ) * self.angle_multiplier

            # for as many trials as specified
            num_mc_accepts = 0
            num_mc_checks = 0
            for trial_num in range( 1, self.trials + 1 ):
                #################
                #### RAMPING ####
                #################
                if self.ramp_sf:
                    # ramp up or down the appropriate scoring terms and get it back to the MonteCarlo object
                    self.sf = ramp_score_weight( self.sf,
                                                 "fa_atr",
                                                 FA_ATR_ORIG,
                                                 trial_num - 1,
                                                 self.trials )
                    self.sf = ramp_score_weight( self.sf,
                                                 "fa_rep",
                                                 FA_REP_ORIG,
                                                 trial_num - 1,
                                                 self.trials )
                    self.mc.score_function( self.sf )
                # print current score
                if self.verbose:
                    print "\nstarting score:", self.watch_sf( working_pose )

                ################################
                #### MAKE SUGAR SMALL MOVES ####
                ################################
                # for as many moves per trial as desired
                for ii in range( self.moves_per_trial ):
                    # pick a random moveable residue (Fc glycan residue except the core GlcNAc)
                    res_num = choice( moveable_res_nums )
                    # apply the SugarSmallMover
                    working_pose.assign( SugarSmallMover( res_num, working_pose, angle_max ) )
                    pmm.apply( working_pose )
                if self.verbose:
                    print "score after SugarSmallMover:", self.watch_sf( working_pose )

                ###################
                #### MINIMIZE #####
                ###################
                # minimize the backbone of the Fc sugars
                #Fc_glycan_min_mover.apply( working_pose )
                #if self.verbose:
                #    print "score after min:", self.watch_sf( working_pose )

                ###############################
                #### ACCEPT OR REJECT MOVE ####
                ###############################
                # accept or reject the total move using the MonteCarlo object
                #if self.mc.boltzmann( working_pose ):
                #    # up the counters and send to pymol
                #    num_mc_accepts += 1
                #    pmm.apply( working_pose )
                #num_mc_checks += 1

                # print out the MC acceptance rate every 3 trials and on the last trial
                #mc_acceptance = round( ( float( num_mc_accepts ) / float( num_mc_checks ) * 100 ), 2 )
                #if trial_num % 3 == 0 or trial_num == self.trials:
                #if self.verbose:
                    #print "Moves made so far:", num_mc_checks,
                    #print "  Moves accepted:", num_mc_accepts,
                    #print "  Acceptance rate:", mc_acceptance

        elif self.make_shear_moves:
            pass
        else:
            pass


        return working_pose
