#!/usr/bin/python
__author__="morganlnance"


'''
The intent for this file is to hold all of the single-call functions for each protocol that is desired to be tested in modeling the 3ay4 glycan. For this to work appropriately, a file or two should hold all of the necessary "worker" functions, and this file should combine those functions into step to create a protocol. Then another script will have all of the needed instantiation information (making dirs, checking files, making the JobDistributor, etc), and will call a single protocol function from here in the JobDistributor

protocol_0: Base protocol. LCM reset starting from ideal for each linkage and sampling within +/- 1 stdev from there. SSM-200 am3 3mpt. am3 applies to each linkage and torsion equally. Core GlcNAc does not move. After LCM reset and before any movements, native omega at branch residues is restored. native_3ay4_Gal_5A_1A_tol.cst is used. Standard fa_intra_rep (0.44) sf. No packing at all, but minimization of just Fc glycan bb (no chi, no branch) before each MonteCarlo call

protocol_1: Comparison protocol. LCM reset starting from ideal for each linkage (no stdev sampling). SSM-200 am3 3mpt. am3 applies to each linkage and torsion equally. Core GlcNAc does not move. After LCM reset and before any movements, native omega at branch residues is restored. native_3ay4_Gal_5A_1A_tol.cst is used. Standard fa_intra_rep (0.44) sf. No packing at all, but minimization of just Fc glycan bb (no chi, no branch) before each MonteCarlo call
'''


# global variables and imports
import sys, os
kT = 0.8



class ProtocolBase:
    def __init__( self, sf, nstruct, dump_dir, util_dir ):
        '''
        :param sf: ScoreFunction for protocol
        :param nstruct: int( number of decoys to be made )
        :param dump_dir: str( /path/to/dump location )
        :param util_dir: str( /path/to/utility files )
        '''
        # default arguments
        self.name = "Protocol"

        # required arguments
        self.sf = sf
        self.nstruct = nstruct
        self.dump_dir = dump_dir
        self.util_dir = util_dir

        # "optional" arguments
        #self.glyco_file = None
        self.verbose = False
        self.watch_for_convergence = False


        '''
        Setup the dump directory
        '''
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



class Protocol_0( ProtocolBase ):
    '''
    Protocol_0: am3, 3mpt, native_3ay4_Gal_5A_1A_tol, SSM-200, LCM reset with +/- stdev, core GlcNAc unmoved, fa_intra_rep sf used (given from outside call), no packing, minimization of Fc glycan before each MC call
    '''
    def __init__( self, sf, nstruct, dump_dir, util_dir ):
        '''
        :param decoy_num: int( which decoy number is this? )
        :param sf: ScoreFunction for protocol
        :param nstruct: int( number of decoys to be made )
        :param dump_dir: str( /path/to/dump location )
        :param util_dir: str( /path/to/utility files )
        '''
        # inherited arguments
        ProtocolBase.__init__( self, sf, nstruct, dump_dir, util_dir )

        # default arguments
        self.name = "Protocol_0"
        self.nmoves = 200
        self.angle_multiplier = 3
        self.moves_per_trial = 3
        self.LCM_reset = True
        self.use_population_ideal_LCM_reset = True
        self.set_native_omega = True
        self.ramp_sf = True
        self.native_constraint_file = "~/pyrosetta_dir/project_constraint_files/native_3ay4_Gal_5A_1A_tol.cst"


    def apply( self, pose ):
        '''
        Base protocol. This is the most generic protocol (as of my PyRosetta-3's most recent abilities 9/26/16). Use to compare how other more "unique" protocols do
        :param pose: Pose
        '''
        #########################
        #### protocol_0 info ####
        #########################
        info_file_details = []
        info_file_details.append( "Native PDB filename:\t\t\t%s\n" %pose.pdb_info().name().split( '/' )[-1] )
        #info_file_details.append( "Sugar filename:\t\t\t\t%s\n" %input_args.glyco_file.split( '/' )[-1] )
        info_file_details.append( "Creating this many decoys:\t\t%s\n" %self.nstruct )
        info_file_details.append( "Number of SugarSmallMove trials:\t%s\n" %self.nmoves )
        info_file_details.append( "Number of SugarSmallMoves per trial:\t%s\n" %self.moves_per_trial )
        info_file_details.append( "LCM reset of Fc glycan?:\t\t%s\n" %self.LCM_reset )
        info_file_details.append( "Use population ideals in LCM reset?:\t%s\n" %self.use_population_ideal_LCM_reset )
        info_file_details.append( "Reset omega torsion back to native?:\t%s\n" %self.set_native_omega )
        info_file_details.append( "Using score ramping?:\t\t\t%s\n" %self.ramp_sf )
        info_file_details.append( "Native constraint file used?:\t\t%s\n" %self.native_constraint_file )
        info_file_details.append( "Angle multiplier used:\t\t\t%s\n" %self.angle_multiplier )
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
        info_file = ''.join( info_file_details )
        print "\n", info_file, "\n"

        # write out the info file with the collected info from above
        info_filename = self.dump_dir + "protocol_0_run.info"
        with open( info_filename, "wb" ) as fh:
            fh.write( "Info for this run of %s\n\n" %__file__ )
            fh.write( info_file )


        #################
        #### IMPORTS ####
        #################
        from rosetta import PyMOL_Mover, MoveMap, MinMover
        from antibody_functions import native_Fc_glycan_nums_except_core_GlcNAc
        from native_3ay4_glycan_modeling_protocol_functions import native_3ay4_Fc_glycan_LCM_reset


        ############################
        #### START THE PROTOCOL ####
        ############################
        # make the PyMOL_Mover
        pmm = PyMOL_Mover()
        pmm.keep_history( True )

        # get the working and native pose (for this script, they are the same thing)
        # native
        native_pose = pose.clone()
        native_pose.pdb_info().name( "protocol_0_native" )
        pmm.apply( native_pose )
        # working
        working_pose = pose.clone()
        working_pose.pdb_info().name( "protocol_0_decoy" )
        pmm.apply( working_pose )


        ###################
        #### LCM RESET ####
        ###################
        working_pose.assign( native_3ay4_Fc_glycan_LCM_reset( input_pose = working_pose, 
                                                              residue_numbers = native_Fc_glycan_nums_except_core_GlcNAc, 
                                                              bb_freedom = True, 
                                                              chi_freedom = False, 
                                                              use_population_ideal_LCM_reset = False ) )
        pmm.apply( working_pose )
        if self.verbose:
            print "score of LCM reset:", self.sf( working_pose )

        # reset the native omega torsion, if desired
        if self.set_native_omega:
            working_pose.set_omega( 221, native_pose.omega( 221 ) )
            working_pose.set_omega( 445, native_pose.omega( 445 ) )


        ############################
        #### MIN MOVER CREATION ####
        ############################
        # make a MoveMap
        min_mm = MoveMap()
        for res_num in native_Fc_glycan_nums_except_core_GlcNAc:
            min_mm.set_bb( res_num, True )
            min_mm.set_chi( res_num, False )
        for branch_point in native_Fc_glycan_nums_except_core_GlcNAc:
            min_mm.set_branches( branch_point, False )

        # make the MinMover
        Fc_glycan_min_mover = MinMover( movemap_in = min_mm,
                                        scorefxn_in = self.sf,
                                        min_type_in = "dfpmin_strong_wolfe",
                                        tolerance_in = 0.01,
                                        use_nb_list_in = True )

        return working_pose
