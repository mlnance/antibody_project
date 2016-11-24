#!/usr/bin/python
__author__ = 'morganlnance'



#####################
#### ALL IMPORTS ####
#####################

# main PyMOL_Mover watcher
from rosetta import PyMOL_Mover

# import extras
import os, sys

# for sugar work
#from rosetta.protocols.carbohydrates import GlycanRelaxMover
#from rosetta.protocols.carbohydrates import LinkageConformerMover
#from rosetta.core.chemical.rings import RingConformer

# for extra scoring functionality
#from rosetta.core.scoring import CA_rmsd



# create global pymol object
# TODO-add visualization options to each relevant function using this global PyMOL_Mover
pmm = PyMOL_Mover()


# global variables
AA_name1_list = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
AA_name3_list = [ "ALA", "CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR" ]
AA_name1_to_name3 = { 'A':"ALA", 'C':"CYS", 'D':"ASP", 'E':"GLU", 'F':"PHE", 'G':"GLY", 'H':"HIS", 'I':"ILE", 'K':"LYS", 'L':"LEU", 'M':"MET", 'N':"ASN", 'P':"PRO", 'Q':"GLN", 'R':"ARG", 'S':"SER", 'T':"THR", 'V':"VAL", 'W':"TRP", 'Y':"TYR" }
AA_name3_to_name1 = { "ALA":'A', "CYS":'C', "ASP":'D', "GLU":'E', "PHE":'F', "GLY":'G', "HIS":'H', "ILE":'I', "LYS":'K', "LEU":'L', "MET":'M', "ASN":'N', "PRO":'P', "GLN":'Q', "ARG":'R', "SER":'S', "THR":'T', "VAL":'V', "TRP":'W', "TYR":'Y' }
all_letters_list = [ 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z' ]
CUTOFF_DISTANCE = 5.0  # used when calculating the number of residue contacts at the interface
PACK_RADIUS = 20.0  # used when making mutations to structures, repacks in this area
PROBE_RADIUS = 1.4  # used for calculating total SASA
kT = 0.8  # used in MonteCarlo and small and shear movers



## Pose numbering information ONLY relevant to PDB 3ay4
# native
native_Fc_chain_A_nums = range( 1, 215 + 1 )
native_Fc_glycan_A_nums = range( 216, 223 + 1 )
native_Fc_glycan_A_nums_except_core_GlcNAc = range( 217, 223 + 1 )
native_Fc_chain_B_nums = range( 224, 439 + 1 )
native_Fc_glycan_B_nums = range( 440, 447 + 1 )
native_Fc_glycan_B_nums_except_core_GlcNAc = range( 441, 447 + 1 )
native_FcR_protein_nums = range( 448, 607 + 1 )
native_FcR_main_glycan_nums = range( 608, 615 + 1 )
native_FcR_three_mer_nums = range( 616, 618 + 1 )
native_FcR_glycan_nums = range( 608, 618 + 1 )
native_branch_points = [ 69, 218, 292, 442, 478, 595, 608, 610 ]
native_Fc_glycan_branch_point_nums = [ 218, 442 ]
native_Fc_glycan_branch_point_nums_with_ASN = [ 69, 218, 292, 442 ]
native_Fc_protein_nums = []
native_Fc_protein_nums.extend( native_Fc_chain_A_nums )
native_Fc_protein_nums.extend( native_Fc_chain_B_nums )
native_Fc_glycan_nums = []
native_Fc_glycan_nums.extend( native_Fc_glycan_A_nums )
native_Fc_glycan_nums.extend( native_Fc_glycan_B_nums )
native_Fc_glycan_nums_except_core_GlcNAc = []
native_Fc_glycan_nums_except_core_GlcNAc.extend( native_Fc_glycan_A_nums_except_core_GlcNAc )
native_Fc_glycan_nums_except_core_GlcNAc.extend( native_Fc_glycan_B_nums_except_core_GlcNAc )
native_order_nums = range( 1, 618 + 1 )
native_Fc_protein_chains = [ 'A', 'B' ]
native_FcR_protein_chains = [ 'C' ]
native_Fc_glycan_chains = [ 'D', 'E', 'F', 'G' ]
native_Fc_glycan_A_chains = [ 'D', 'E' ]
native_Fc_glycan_B_chains = [ 'F', 'G' ]
native_FcR_glycan_chains = [ 'H', 'I', 'J', 'K' ]

# glycosylated decoy
decoy_Fc_chain_A_nums = range( 1, 215 + 1 )
decoy_Fc_glycan_A_nums = range( 603, 610 + 1)
decoy_Fc_chain_B_nums = range( 216, 431 + 1 )
decoy_Fc_glycan_B_nums = range( 611, 618 + 1)
decoy_FcR_protein_nums = range( 432, 591 + 1)
decoy_FcR_main_glycan_nums = range( 592, 598 + 1 )
decoy_FcR_three_mer_nums = range( 599, 602 + 1 )
decoy_FcR_glycan_nums = range( 592, 602 + 1 )
decoy_Fc_protein_nums = []
decoy_Fc_protein_nums.extend( decoy_Fc_chain_A_nums )
decoy_Fc_protein_nums.extend( decoy_Fc_chain_B_nums )
decoy_Fc_glycan_nums = []
decoy_Fc_glycan_nums.extend( decoy_Fc_glycan_A_nums )
decoy_Fc_glycan_nums.extend( decoy_Fc_glycan_B_nums )
decoy_order_nums = []
decoy_order_nums.extend( decoy_Fc_chain_A_nums )
decoy_order_nums.extend( decoy_Fc_glycan_A_nums )
decoy_order_nums.extend( decoy_Fc_chain_B_nums )
decoy_order_nums.extend( decoy_Fc_glycan_B_nums )
decoy_order_nums.extend( decoy_FcR_protein_nums )
decoy_order_nums.extend( decoy_FcR_main_glycan_nums )
decoy_order_nums.extend( decoy_FcR_three_mer_nums )
decoy_Fc_protein_chains = [ 'A', 'B' ]
decoy_FcR_protein_chains = [ 'C' ]
decoy_Fc_glycan_chains = [ 'H', 'I', 'J', 'K' ]
decoy_Fc_glycan_A_chains = [ 'H', 'I' ]
decoy_Fc_glycan_B_chains = [ 'J', 'K' ]
decoy_FcR_glycan_chains = [ 'D', 'E', 'F', 'G' ]



# make an appropriate dictionary map
native_to_decoy_res_map = {}
for ii in range( len( native_order_nums ) ):
    native_to_decoy_res_map[ native_order_nums[ii] ] = decoy_order_nums[ii]
decoy_to_native_res_map = {}
for ii in range( len( decoy_order_nums ) ):
    decoy_to_native_res_map[ decoy_order_nums[ii] ] = native_order_nums[ii]




# adds these two functions to FoldTree so that now FoldTree has the functions and can use elsewhere
def _new_loop( loop, loop_buffer = 1):
    begin = loop.start() - loop_buffer
    end = loop.stop() + loop_buffer
    cut = loop.cut()
    new_jump( begin, end, cut )


def _new_loops( loops ):
    for ii in range( 1, loops.num_loop() + 1 ):
        loop = loops[ ii ]
        new_loop( loop )

from rosetta import FoldTree
FoldTree.new_loop = _new_loop
FoldTree.new_loops = _new_loops




##########################
#### WORKER FUNCTIONS ####
##########################

def initialize_rosetta( constant_seed = False, debug = False ):
    """
    Initialize Rosetta and mute basic, core, and protocols.
    If constant_seed == True, use default constant seed 1111111
    If debug == True, use default constant seed and do not mute Rosetta
    """
    from rosetta import init


    print "Initializing Rosetta with sugar flags"

    # makes Rosetta quiet and sugar I/O ready
    #init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records" )
    if constant_seed:
        init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -write_pdb_link_records -constant_seed" )
    elif debug:
        init( extra_options="-include_sugars -override_rsd_type_limit -write_pdb_link_records -constant_seed -out:level 400" )
    else:
        init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -write_pdb_link_records" )



def load_pose( pose_filename ):
    """
    Load pose from a filename
    :param pose_filename: str( /path/to/pose/filename )
    :return: a Rosetta Pose
    """
    # imports
    from rosetta import Pose, pose_from_file, FoldTree
    
    # create Pose object from filename
    pose = Pose()
    pose_from_file( pose, pose_filename )
    
    # clean up the name of the pose
    pose_name = pose.pdb_info().name()
    pose_name = pose_name.split( '/' )[-1]
    pose.pdb_info().name( pose_name )

    # store the original FoldTree and add empty loops for use later
    pose.orig_fold_tree = FoldTree( pose.fold_tree() )
    pose.loops = None
    pose.loops_file = None

    return pose



class DataHolder:
    pass



class hold_chain_and_res_designations_3ay4:
    # store residue and chain designation information for both native and decoy for 3ay4
    def __init__( self ):
        self.name = "3ay4"

    def native( self ):
        # native information
        self.native_Fc_chain_A_nums = native_Fc_chain_A_nums
        self.native_Fc_glycan_A_nums = native_Fc_glycan_A_nums
        self.native_Fc_chain_B_nums = native_Fc_chain_B_nums
        self.native_Fc_glycan_B_nums = native_Fc_glycan_B_nums
        self.native_FcR_protein_nums = native_FcR_protein_nums
        self.native_FcR_main_glycan_nums = native_FcR_main_glycan_nums
        self.native_FcR_three_mer_nums = native_FcR_three_mer_nums
        self.native_FcR_glycan_nums = native_FcR_glycan_nums
        self.native_Fc_protein_nums = native_Fc_protein_nums
        self.native_Fc_glycan_nums = native_Fc_glycan_nums
        self.native_Fc_glycan_nums_except_core_GlcNAc = native_Fc_glycan_nums_except_core_GlcNAc
        self.native_Fc_glycan_branch_point_nums = native_Fc_glycan_branch_point_nums
        self.native_Fc_glycan_branch_point_nums_with_ASN = native_Fc_glycan_branch_point_nums_with_ASN
        self.native_order_nums = native_order_nums
        self.native_Fc_protein_chains = native_Fc_protein_chains
        self.native_FcR_protein_chains = native_FcR_protein_chains
        self.native_Fc_glycan_chains = native_Fc_glycan_chains
        self.native_Fc_glycan_A_chains = native_Fc_glycan_A_chains
        self.native_Fc_glycan_B_chains = native_Fc_glycan_B_chains
        self.native_FcR_glycan_chains = native_FcR_glycan_chains

    # decoy designations are commented out for now as my protocol does not use the glycosylate_pose function
    '''
    def decoy( self ):
        # decoy information
        self.Fc_protein_nums = decoy_Fc_protein_nums
        self.FcR_main_glycan_nums = decoy_FcR_main_glycan_nums
        self.Fc_glycan_A_nums = decoy_Fc_glycan_A_nums
        self.Fc_glycan_A_chains = decoy_Fc_glycan_A_chains
        self.Fc_glycan_B_chains = decoy_Fc_glycan_B_chains
        """
        self.decoy_Fc_chain_A_nums = decoy_Fc_chain_A_nums
        self.decoy_Fc_glycan_A_nums = decoy_Fc_glycan_A_nums
        self.decoy_Fc_chain_B_nums = decoy_Fc_chain_B_nums
        self.decoy_Fc_glycan_B_nums = decoy_Fc_glycan_B_nums
        self.decoy_FcR_protein_nums = decoy_FcR_protein_nums
        self.decoy_FcR_main_glycan_nums = decoy_FcR_main_glycan_nums
        self.decoy_FcR_three_mer_nums = decoy_FcR_three_mer_nums
        self.decoy_FcR_glycan_nums = decoy_FcR_glycan_nums
        self.decoy_Fc_protein_nums = decoy_Fc_protein_nums
        self.decoy_Fc_glycan_nums = decoy_Fc_glycan_nums
        self.decoy_order_nums = decoy_order_nums
        self.decoy_Fc_protein_chains = decoy_Fc_protein_chains
        self.decoy_FcR_protein_chains = decoy_FcR_protein_chains
        self.decoy_Fc_glycan_chains = decoy_Fc_glycan_chains
        self.decoy_FcR_glycan_chains = decoy_FcR_glycan_chains
        """
    '''



def align_sugar_virtual_atoms( input_pose ):
    """
    Set the xyz coordinates of each sugar virtual atom to the corresponding atom it is shadowing
    :param input_pose: Pose
    :return: Pose
    """
    from rosetta import AtomID


    # copy the input pose
    pose = input_pose.clone()

    # for each residue in the pose
    for res in pose:
        # if it is a sugar
        if res.is_carbohydrate():
            # if this sugar resiude has shadow atoms ( they all should )
            if res.type().has_shadow_atoms():
                # for each atom in the sugar
                for atom_num in range( 1, res.natoms()+1 ):
                    # get the corresponding real atom it is following
                    atom_being_shadowed = res.type().atom_being_shadowed( atom_num )
                    # if this atom is being shadowed
                    # atom_being_shadowed returns 0 if the atom_num of interest isn't shadowing anyone
                    if atom_being_shadowed != 0:
                        # get the AtomID object for this virtual atom
                        atom_id = AtomID( atom_num, res.seqpos() )
                        # set the xyz of the virtual atom to that of the real atom being shadowed
                        pose.set_xyz( atom_id, res.atom( atom_being_shadowed ).xyz() )

                        # update the residue conformation in a hacky way
                        # calling it this way just makes Pose conformation update itself
                        pose.conformation().residue( res.seqpos() )
    '''
    # check that it worked
    for res in pose:
        if res.is_carbohydrate():
            if res.type().has_shadow_atoms():
                for atom_num in range(1, res.natoms()+1 ):
                    atom_being_shadowed = res.type().atom_being_shadowed( atom_num )
                    if atom_being_shadowed != 0:
                        print list(res.atom(atom_num).xyz()) == list(res.atom(atom_being_shadowed).xyz()), atom_num, atom_being_shadowed, res.seqpos()
    '''

    return pose



'''
def stepwise_3ay4_torsion_sampler_using_LCM_ideals( sf, residues, input_pose, num_cycles = 3, num_stdev = 1 ):
    """
    TODO: make this randomly choose which torsion to sample based on how many it has ( this will determine the order. Phi then psi? Or vice versa?
    TODO: add a pack after each torsion change? Of just that single residue?
    TODO: make this check the current torsion and find the nearest ideal LCM population cluster
    TODO: make this work by just checking what kind of linkage pair is being sampled
    Iterate through <residues> one-by-one and adjust the phi, psi, and (if any) omega within +/- the standard devation of the the main LCM ideal population cluster
    Uses get_3ay4_ideal_LCM_phi_psi_omega_info to get the needed phi, psi, omega data
    :param sf: ScoreFunction
    :param residues: list( residue Pose positions )
    :param input_pose: Pose
    :param num_cycles: int( how many cycles of the torsion sampler do you want to have done? ) Default = 3
    :param num_stdev: int( how many standard deviations do you want to sample around the torsion value? ) Default = 1
    :return: Pose
    """
    # imports
    from rosetta.core.pose.carbohydrates import get_reference_atoms_for_1st_omega
    from random import choice


    # copy the <input_pose>
    pose = input_pose.clone()

    # score the Pose to gain access to its .energies() object
    sf( pose )

    # get the needed ideal LCM data
    phi_ideal_dict, phi_stdev_dict, psi_ideal_dict, psi_stdev_dict, omega_ideal_dict, omega_stdev_dict = get_3ay4_ideal_LCM_phi_psi_omega_info()

    # num cycles
    for ii in range( 1 ):
        for res_num in residues:
            # for each residue in the <residues> list given
            for jj in range( 1, num_cycles + 1 ):
                # depending on which iteration we're on, get the starting torsion value up to a certain number of sigfigs
                # first iteration == finding an integer
                if jj == 1:
                    # get integer version of the current torsions for the first iteration
                    # since we're doing a range above and below the start, the start will already be included in the range
                    start_phi = int( pose.phi( res_num ) )

                    # get the integer version of the ideal LCM standard deviations associated with this residue
                    phi_stdev = int( phi_stdev_dict[ res_num ] )
                    phi_range = range( start_phi - phi_stdev, start_phi + phi_stdev + 1 )

                # second iteration == finding tenths place
                elif jj == 2:
                    # an integer-valued torsion has been found ( or we're back at the starting torsion if that was the best )
                    # reaffirm we have an integer torsion, or, if we went back to the original, get the integer value of that torsion
                    start_phi = int( pose.phi( res_num ) )
                    phi_range = [ start_phi ]

                    # go over and below the starting phi in steps of 0.1 up until 0.9 above and below
                    phi_range.extend( [ round( start_phi - ( x / 10.0 ), 1 ) for x in range( 1, 10 ) ] )
                    phi_range.extend( [ round( start_phi + ( x / 10.0 ), 1 ) for x in range( 1, 10 ) ] )

                # third iteration == finding hundreths place
                else:
                    # a torsion up to the tenths place has been found ( or we're back at the starting torsion if that was the best )
                    start_phi = round( pose.phi( res_num ), 1 )
                    phi_range = [ start_phi ]

                    # go over and below the starting phi in steps of 0.01 from 0.01 to 0.99 above and below
                    phi_range.extend( [ round( start_phi - ( x / 100.0 ), 2 ) for x in range( 1, 100 ) ] )
                    phi_range.extend( [ round( start_phi + ( x / 100.0 ), 2 ) for x in range( 1, 100 ) ] )

                # get the starting residue total energy
                start_res_tot_E = pose.energies().residue_total_energy( res_num )

                # sample each integer torsion value within +/- the ideal standard devation while keeping residue_total_energy data
                phi_to_res_tot_E = {}

                # store a bool to use to see if this residue has an omega torsion
                # 0 = False = no omega, 4 = True = has omega
                res_has_omega = bool( len( get_reference_atoms_for_1st_omega( pose, res_num ) ) )

                # set the new torsion, score the Pose to gain access to its updated .energies() object, and store the new residue total energy
                # find the best torsion in the given range
                for torsion in phi_range:
                    pose.set_phi( res_num, torsion )
                    #sf( pose )
                    # phi_to_res_tot_E[ torsion ] = pose.energies().residue_total_energy( res_num )
                    phi_to_res_tot_E[ torsion ] = sf( pose )
                
                # determine which torsions resulted in the lowest energy for each residue
                # if more than one torsion (somehow) resulted in the same lowest energy, randomly choose one of the torsions to continue with
                # first, get the minimum(s) value(s) from each torsion-to-energy dict
                phi_min = min( phi_to_res_tot_E.values() )

                # get the torsion(s) associated with each minimum value(s)
                # instantiate the lists
                best_phi_torsion_list = []

                # use the minimum torsions to find the corresponding torsion key
                # phi
                for torsion in phi_to_res_tot_E.keys():
                    if phi_to_res_tot_E[ torsion ] == phi_min:
                        best_phi_torsion_list.append( torsion )

                # choose the best torsion from the lists. Again, doing this in case two or more torsions resulted in the same low E
                best_phi_torsion = choice( best_phi_torsion_list )

                # or keep the current torsion if that was a better total energy
                if phi_to_res_tot_E[ best_phi_torsion ] > start_res_tot_E:
                    best_phi_torsion = start_phi

                # set the Pose to the best torsions
                pose.set_phi( res_num, best_phi_torsion )
                print best_phi_torsion, sf( pose ), jj

    return pose



def get_3ay4_stepwise_SugarSmallMover_LCM_reset_phi_psi_omega_info( Bsubtractor = 0 ):
    """
    Stepwise means focusing on just one residue at a time where all upper residues are deleted ( this was done manually )
    LCM reset means resetting the target residue using the LCM reset data where all populations were sampled ( using population weights ) within 1 stdev
    Protocol used SugarSmallMover with am2 and 2mpt
    Each lower residue was set to the "ideal" value found by the protocol
    The data is collected from the .fasc file and printed out from the plotter script found in metrics/stepwise_by_removal
    :param Bsubtractor: this is terribly hacky, but please give me the number to subtract from the standard Pose number for the side B glycan
    :return: dict( phi_ideal ), dict( phi_stdev ), dict( psi_ideal ), dict( psi_stdev ), dict( omega_ideal ), dict( omega_stdev )
    """
    phi_data = { 217: -81.7946306335, 
                 218: -98.0859420888, 
                 219: 72.4906160885, 
                 220: -69.2819693849, 
                 221: 70.2338462561, 
                 222: -58.6313776519, 
                 441 - Bsubtractor: -79.3869827081, 
                 442 - Bsubtractor: -122.979289015, 
                 443 - Bsubtractor: 78.8230022219, 
                 444 - Bsubtractor: -97.7081094911, 
                 445 - Bsubtractor: 75.9149018719, 
                 446 - Bsubtractor: -66.571881752 }

    phi_stdev = { 217: 8.02969709033, 
                  218: 6.60818707069, 
                  219: 9.10612718064, 
                  220: 46.5191149619, 
                  221: 12.8199926237, 
                  222: 55.451733708, 
                  441 - Bsubtractor: 9.48701851985, 
                  442 - Bsubtractor: 11.2109158286, 
                  443 - Bsubtractor: 7.96196655884, 
                  444 - Bsubtractor: 51.203299921, 
                  445 - Bsubtractor: 13.5208014489, 
                  446 - Bsubtractor: 54.8486064208 }


    psi_data = { 217: 87.113986145, 
                 218: 89.2904142081, 
                 219: -134.691402013, 
                 220: -90.5627178154, 
                 221: 168.404809533, 
                 222: -138.336042013, 
                 441 - Bsubtractor: 102.86554493, 
                 442 - Bsubtractor: 87.8167745987, 
                 443 - Bsubtractor: -116.599767556, 
                 444 - Bsubtractor: -91.2041858112, 
                 445 - Bsubtractor: 157.966096193, 
                 446 - Bsubtractor: -101.66542411 }

    psi_stdev = { 217: 6.21243150763, 
                  218: 7.97750106414, 
                  219: 18.3655776045, 
                  220: 16.8077078643, 
                  221: 150.990946995, 
                  222: 37.0552148656, 
                  441 - Bsubtractor: 6.85640117404, 
                  442 - Bsubtractor: 16.6282496517, 
                  443 - Bsubtractor: 17.6863856501, 
                  444 - Bsubtractor: 16.7236056227, 
                  445 - Bsubtractor: 155.39189028, 
                  446 - Bsubtractor: 19.9466136736 }

    omega_data = { 217: 0.0, 
                   218: 0.0, 
                   219: 0.0, 
                   220: 0.0, 
                   221: 66.0738407506, 
                   222: 0.0, 
                   441 - Bsubtractor: 0.0, 
                   442 - Bsubtractor: 0.0, 
                   443 - Bsubtractor: 0.0,
                   444 - Bsubtractor: 0.0, 
                   445 - Bsubtractor: 169.819545766, 
                   446 - Bsubtractor: 0.0 }

    omega_stdev = { 217: 0.0, 
                    218: 0.0, 
                    219: 0.0, 
                    220: 0.0, 
                    221: 111.972066698, 
                    222: 0.0, 
                    441 - Bsubtractor: 0.0, 
                    442 - Bsubtractor: 0.0, 
                    443 - Bsubtractor: 0.0, 
                    444 - Bsubtractor: 0.0, 
                    445 - Bsubtractor: 121.68743512, 
                    446 - Bsubtractor: 0.0 }

    return phi_data, phi_stdev, psi_data, psi_stdev, omega_data, omega_stdev



def set_3ay4_Fc_glycan_except_core_GlcNAc_to_stepwise_SugarSmallMover_LCM_reset_phi_psi_omega( input_pose, use_stdev = False, Bsubtractor = 0 ):
    """
    Set the Fc glycan (intended for 3ay4) to the values collected from the stepwise SugarSmallMover with LCM reset protocol
    If you want to set phi/psi/omega to within +/- the standard deviation, set <use_stdev> to True
    :param input_pose: Pose
    :param use_stdev: bool( do you want to sample within +/- standard deviation of the phi/psi/omega values from the stepwise LCM reset  data? ) Default = False
    :param Bsubtractor: this is terribly hacky, but please give me the number to subtract from the standard Pose number for the side B glycan
    :return: Pose
    """
    # imports
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg


    # make a copy of the <input_pose>
    pose = input_pose.clone()

    # get the phi, psi, and omega data from the getter function
    phi_data, phi_stdev, psi_data, psi_stdev, omega_data, omega_stdev = get_3ay4_stepwise_SugarSmallMover_LCM_reset_phi_psi_omega_info( Bsubtractor = Bsubtractor )

    # for each residue in phi_data, psi_data, and omega_data, set the appropriate torsions
    # sample within the standard deviation too, if desired
    for res in phi_data.keys():
        if use_stdev:
            # this method is pulled from core.pose.carbohydrates.util::set_dihedrals_from_linkage_conformer_data
            # periodic_range( mean - sd + rg().uniform() * sd * 2, 360 ) not sure what it means though
            new_phi = periodic_range( phi_data[ res ] - phi_stdev[ res ] + rg().uniform() * phi_stdev[ res ] * 2, 360.0 )
            pose.set_phi( res, new_phi )
        else:
            pose.set_phi( res, phi_data[ res ] )
    for res in psi_data.keys():
        if use_stdev:
            # this method is pulled from core.pose.carbohydrates.util::set_dihedrals_from_linkage_conformer_data
            # periodic_range( mean - sd + rg().uniform() * sd * 2, 360 ) not sure what it means though
            new_psi = periodic_range( psi_data[ res ] - psi_stdev[ res ] + rg().uniform() * psi_stdev[ res ] * 2, 360.0 )
            pose.set_psi( res, new_psi )
        else:
            pose.set_psi( res, psi_data[ res ] )
    for res in omega_data.keys():
        if use_stdev:
            # this method is pulled from core.pose.carbohydrates.util::set_dihedrals_from_linkage_conformer_data
            # periodic_range( mean - sd + rg().uniform() * sd * 2, 360 ) not sure what it means though
            new_omega = periodic_range( omega_data[ res ] - omega_stdev[ res ] + rg().uniform() * omega_stdev[ res ] * 2, 360.0 )
            pose.set_omega( res, new_omega )
        else:
            pose.set_omega( res, omega_data[ res ] )

    return pose



def get_3ay4_ideal_LCM_phi_psi_omega_info():
    """
    TODO: make this work by just checking what kind of linkage pair is being checked
    Data is pulled from the highest population cluster from default.table in database/chemical/carbohydrates/linkage_conformers
    Returns data for native_Fc_glycan_nums_except_core_GlcNAc
    :return: dict( residue: [ phi, psi, omega ] ), dict( residue: [ phi stdev, psi stdev, omega stdev ] )
    """
    # data pulled from LCM table default.table
    # key: residue num, value: [ phi, psi, omega ]
    ideal_LCM_phi_psi_omega = { 217: [ -75.9, 119.0, 0.0 ], 
                                218: [ -86.5, 110.7, 0.0 ], 
                                219: [ 71.5, -120.6, 0.0 ], 
                                220: [ -80.1, -97.6, 0.0 ], 
                                221: [ 67.0, 178.5, 186 ], # actually the second-most ideal, but the omega is better
                                222: [ -80.1, -97.6, 0.0 ], 
                                223: [ -71.4, 132.2, 0.0 ], 
                                441: [ -75.9, 119.0, 0.0 ], 
                                442: [ -86.5, 110.7, 0.0 ], 
                                443: [ 71.5, -120.6, 0.0 ], 
                                444: [ -80.1, -97.6, 0.0 ], 
                                445: [ 67.0, 178.5, 186 ], # actually the second-most ideal, but the omega is better
                                446: [ -80.1, -97.6, 0.0 ], 
                                447: [ -71.4, 132.2, 0.0 ], 
                                }

    ideal_LCM_phi_psi_omega_stdev = { 217: [ 11.6, 15.4, 0.0 ], 
                                      218: [ 11.6, 19.4, 0.0 ], 
                                      219: [ 8.8, 16.8, 0.0 ], 
                                      220: [ 12.6, 22.3, 0.0 ], 
                                      221: [ 10.5, 13.7, 12.8 ], 
                                      222: [ 12.6, 22.3, 0.0 ], 
                                      223: [ 10.9, 7.4, 0.0 ], 
                                      441: [ 11.6, 15.4, 0.0 ], 
                                      442: [ 11.6, 19.4, 0.0 ], 
                                      443: [ 8.8, 16.8, 0.0 ], 
                                      444: [ 12.6, 22.3, 0.0 ], 
                                      445: [ 10.5, 13.7, 12.8 ], 
                                      446: [ 12.6, 22.3, 0.0 ], 
                                      447: [ 10.9, 7.4, 0.0 ], 
                                      }

    return ideal_LCM_phi_psi_omega, ideal_LCM_phi_psi_omega_stdev
'''


def set_3ay4_Fc_glycan_except_core_GlcNAc_to_ideal_LCM_phi_psi_omega( input_pose, use_ideal_stdev = False, set_3_D_and_F_phi_to_native = False ):
    """
    Set the Fc glycan (intended for 3ay4) to the ideal values from the LCM data found in default.table in database/chemical/carbohydrates/linkage_conformers
    If you want to use any native values (such as the Man branch), you must give a native pose as well
    If you want to set phi/psi/omega to within +/- the ideal standard deviation using LCM data, set <use_ideal_stdev> to True
    :param input_pose: Pose
    :param use_ideal_stdev: bool( do you want to sample within +/- <ideal_stdev> of the phi/psi/omega values from the LCM data? ) Default = False
    :param set_3_D_and_F_phi_to_native: bool( do you want to set residue 3 on chains D and F to the native phi value? ) Default = False
    :return: Pose
    """
    # imports
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg
    from rosetta.core.id import phi_dihedral, psi_dihedral, omega_dihedral
    from rosetta.core.pose.carbohydrates import set_glycosidic_torsion


    # data pulled from LCM table
    # returned as residue: phi, psi, omega and residue: phi, psi, omega stdevs
    ideal_LCM_phi_psi_omega, ideal_LCM_phi_psi_omega_stdev = get_3ay4_ideal_LCM_phi_psi_omega_info()

    # get a copy of the input_pose
    pose = input_pose.clone()

    # set the phi, psi, and omega to ideal data
    for glyc_num in native_Fc_glycan_nums_except_core_GlcNAc:
        # set residue 3 on chain D and F to native phi value, if desired
        if set_3_D_and_F_phi_to_native is True:
            # chain D Man branch
            if glyc_num == 218:
                # pulled from fa_intra_rep low E native
                set_glycosidic_torsion( phi_dihedral, pose, glyc_num, -109.917 )
            # chain F Man branch
            elif glyc_num == 442:
                # pulled from fa_intra_rep low E native
                set_glycosidic_torsion( phi_dihedral, pose, glyc_num, -122.893 )
            # the other Fc glycan residues
            else:
                # else pull from LCM phi data dict for this particular residue
                # determine the new phi value using standard deviation from ideal or not
                if use_ideal_stdev:
                    # this method is pulled from core.pose.carbohydrates.util::set_dihedrals_from_linkage_conformer_data
                    # periodic_range( mean - sd + rg().uniform() * sd * 2, 360 ) not sure what it means though
                    new_phi = periodic_range( ideal_LCM_phi_psi_omega[ glyc_num ][0] - ideal_LCM_phi_psi_omega_stdev[ glyc_num ][0] + rg().uniform() * ideal_LCM_phi_psi_omega_stdev[ glyc_num ][0] * 2, 360.0 )
                # otherwise, just use the ideal phi
                else:
                    new_phi = ideal_LCM_phi_psi_omega[ glyc_num ][0]
                # set the new phi
                set_glycosidic_torsion( phi_dihedral, pose, glyc_num, new_phi )
        # otherwise, use LCM phi data dict values for all residues
        else:
            # determine the new phi value using standard deviation from ideal or not
            if use_ideal_stdev:
                # this method is pulled from core.pose.carbohydrates.util::set_dihedrals_from_linkage_conformer_data
                # periodic_range( mean - sd + rg().uniform() * sd * 2, 360 ) not sure what it means though
                new_phi = periodic_range( ideal_LCM_phi_psi_omega[ glyc_num ][0] - ideal_LCM_phi_psi_omega[ glyc_num ][0] + rg().uniform() * ideal_LCM_phi_psi_omega_stdev[ glyc_num ][0] * 2, 360.0 )
            # otherwise, just use the ideal phi
            else:
                new_phi = ideal_LCM_phi_psi_omega[ glyc_num ][0]
            # set the new phi
            set_glycosidic_torsion( phi_dihedral, pose, glyc_num, new_phi )

        # use LCM data dict values for all residues for psi and omega
        # determine the new psi and omega value using standard deviation from ideal or not
        if use_ideal_stdev:
            # this method is pulled from core.pose.carbohydrates.util::set_dihedrals_from_linkage_conformer_data
            # periodic_range( mean - sd + rg().uniform() * sd * 2, 360 ) not sure what it means though
            new_psi = periodic_range( ideal_LCM_phi_psi_omega[ glyc_num ][1] - ideal_LCM_phi_psi_omega_stdev[ glyc_num ][1] + rg().uniform() * ideal_LCM_phi_psi_omega_stdev[ glyc_num ][1] * 2, 360.0 )
            new_omega = periodic_range( ideal_LCM_phi_psi_omega[ glyc_num ][2] - ideal_LCM_phi_psi_omega_stdev[ glyc_num ][2] + rg().uniform() * ideal_LCM_phi_psi_omega_stdev[ glyc_num ][2] * 2, 360.0 )
        # otherwise, just use the ideal psi and omega
        else:
            new_psi = ideal_LCM_phi_psi_omega[ glyc_num ][1]
            new_omega = ideal_LCM_phi_psi_omega[ glyc_num ][2]
        # set the new phi and omega
        set_glycosidic_torsion( psi_dihedral, pose, glyc_num, new_psi )
        # setting the omega to 0 doesn't do anything to glycans if they don't have an omega
        set_glycosidic_torsion( omega_dihedral, pose, glyc_num, new_omega )

    return pose



# TODO: finish this function! There is a function that must do this buried in Rosetta as CarbohydrateInfoManager can get this data from default.table
def get_ideal_LCM_phi_psi_omega_info( linkage_conformer_filename, verbose = False ):
    """
    Pulls out ideal phi/psi/omega according to linkage type from a <linkage_conformer_filename> data file
    This file should be column/tab delimited
    :param linkage_conformer_filename: str( /path/to/linkage conformer data. See default.table in database/chemical/carbohydrates/linkage_conformers
    :param verbose: bool( do you also want to print out the data? ) Default = False
    :return:
    """
    # imports
    import csv, os
    from rosetta.core.chemical.carbohydrates import CarbohydrateInfo


    # check that linkage_conformer_filename exists as a file
    if not os.path.isfile( linkage_conformer_filename ):
        print "\nArgument error. The path to <linkage_conformer_filename> ( %s ) does not exist. Returning None." %linkage_conformer_filename
        return None

    # pull out the linkage conformer data lines
    # the file should be delimited by tabs (\t)
    linkage_conformer_lines = []
    with open( linkage_conformer_filename, "rb" ) as fh:
        # read the file separated by columns
        reader = csv.reader( fh, delimiter="\t" )
        for line in reader:
            # skip empty lines and lines that aren't separated by columns (likely commented lines)
            # csv reader returns each line as a list
            if line != [] and len( line ) != 1:
                linkage_conformer_lines.append( line )

    # the first line in the list should be the headers, the rest should be the actual data
    header_info = linkage_conformer_lines[ 0 ]
    linkage_conformer_data = linkage_conformer_lines[ 1: ]

    # create a dictionary for the header info to the index number of the column
    index_to_header_dict = {}
    for ii in range( len( header_info ) ):
        # this is a commented line, so replace the '#' and any spaces with ''
        column_name = header_info[ ii ].replace( '#', '' ).replace( ' ', '' )

        # add the column name and its index number to the dictionary
        index_to_header_dict[ ii ] = column_name

    # create a dictionary for each line of data, adding to each key a list of the data
    # dict: key = non-reducing_reducing : value = [ [ population, phi_mean, phi_stdev, psi_mean, psi_stdev, omega_mean, omega_stdev, omega2_mean, omega2_stdev
    nonred_to_red_data_dict = {}
 
    return header_info, linkage_conformer_data



def get_ideal_SugarBB_phi_psi_info( sugar_num, input_pose, verbose = False ):
    """
    CHI Phi: determines if residue is axial or equatorial at its anomeric position (ie. if alpha or beta sugar) and returns appropriate statistic
    CHI Psi: determines if residue is axial or equatorial at the attachment position to its parent residue and returns appropriate statistic
    See core::pose::carbohydrates::util::get_linkage_type_for_residue_for_CHI and core::chemical::carbohydrates::LinkageType for more information
    :param sugar_num: int( Pose number for sugar residue of interest )
    :param input_pose: Pose
    :param verbose: bool( do you also want to print out the data? ) Default = False
    :return: list( ideal_phi_major, ideal_phi_stdev_major, ideal_phi_minor, ideal_phi_stdev_minor, ideal_psi, ideal_psi_stdev, anomeric_position, linkage_type )
    :return: None if error occured, NA if no information
    """
    # imports
    from rosetta.core.id import phi_dihedral, psi_dihedral
    from rosetta.core.pose.carbohydrates import get_linkage_type_for_residue_for_CHI
    # for C++ to Python value converter
    from rosetta.core.chemical.carbohydrates import LinkageType


    # instantiate the data holders
    # data taken from Jason's RosettaCarbohydrates Tutorial/Demo 2 (Scoring & Packing)
    # this describes the sugar_bb (SugarBB) data taken from the CHI energy function (Grant & Woods, Curr. Opin. Struct. Biol, 2014, 28C, 47-55
    # specific values taken from http://www.ncbi.nlm.nih.gov/pubmed/24375430 (Nivedha & Woods, J Comput Chem, 2014, 35(7), 526-39
    # -360 for values greater than 180 for the sake of Rosetta
    # phi information
    ideal_alpha_phi = 72.5
    ideal_alpha_phi_stdev = 17.5  # ( 55-90)
    ideal_beta_phi_major = 290 - 360
    ideal_beta_phi_stdev_major = 15 # (275-305) - 360
    ideal_beta_phi_minor = 55
    ideal_beta_phi_stdev_minor = 15 # (40-70)

    # psi information
    ideal_2ax_3eq_4ax_psi = 242.5 - 360
    ideal_2ax_3eq_4ax_psi_stdev = 42.5 # (200-285) - 360
    ideal_2eq_3ax_4eq_psi = 112.5
    ideal_2eq_3ax_4eq_psi_stdev = 37.5 # (75-150)

    # check that the residue number passed actually exists
    try:
        sugar_res = input_pose.residue( sugar_num )
    except RuntimeError:
        print "\nIt appears that residue number", sugar_num, "does not actually exist in the Pose. Check your input. Returning None."
        return None

    # check that the residue passed is actually a sugar
    if not sugar_res.is_carbohydrate():
        print "\nIt appears that residue number", sugar_num, "is not actually a carbohydrate residue. Check your input. Returning None."
        return None

    # return the different combinations of ideal phi/psi statistics based on the anomeric position (phi) and the linkage position (psi)
    anomeric_position = str( get_linkage_type_for_residue_for_CHI( phi_dihedral, sugar_res, input_pose ) )
    linkage_type = str( get_linkage_type_for_residue_for_CHI( psi_dihedral, sugar_res, input_pose ) )

    # Python isn't smart enough to store two values created in the source code, so replace linkage_type if necessary
    # N_LINK_TYPES = 4, which is also what _2EQ_3AX_4EQ_LINKS equals. So it gets overwritten apparently from C++ --> Python
    if linkage_type == "N_LINK_TYPES":
        linkage_type = "_2EQ_3AX_4EQ_LINKS"

    # alpha sugar with _2AX_3EQ_4AX_LINKS LinkageType
    if anomeric_position == "ALPHA_LINKS" and linkage_type == "_2AX_3EQ_4AX_LINKS":
        if verbose:
            print "Residue number:", sugar_num
            print "  ", anomeric_position, linkage_type
            print "   ideal phi major:", ideal_alpha_phi
            print "\tideal phi stdev major:", ideal_alpha_phi_stdev
            print "   ideal phi minor:", "NA", 
            print "\tideal phi stdev minor:", "NA"
            print "      actual phi:", input_pose.phi( sugar_num )
            print "   ideal psi:", ideal_2ax_3eq_4ax_psi, 
            print "\tideal psi stdev:", ideal_2ax_3eq_4ax_psi_stdev
            print "      actual psi:", input_pose.psi( sugar_num )
        return [ ideal_alpha_phi, ideal_alpha_phi_stdev, "NA", "NA", ideal_2ax_3eq_4ax_psi, ideal_2ax_3eq_4ax_psi_stdev, anomeric_position, linkage_type ]

    # alpha sugar with _2EQ_3AX_4EQ_LINKS LinkageType
    elif anomeric_position == "ALPHA_LINKS" and linkage_type == "_2EQ_3AX_4EQ_LINKS":
        if verbose:
            print "Residue number:", sugar_num
            print "  ", anomeric_position, linkage_type
            print "   ideal phi major:", ideal_alpha_phi
            print "\tideal phi stdev major:", ideal_alpha_phi_stdev
            print "   ideal phi minor:", "NA", 
            print "\tideal phi stdev minor:", "NA"
            print "      actual phi:", input_pose.phi( sugar_num )
            print "   ideal psi:", ideal_2eq_3ax_4eq_psi, 
            print "\tideal psi stdev:", ideal_2eq_3ax_4eq_psi_stdev
            print "      actual psi:", input_pose.psi( sugar_num )
        return [ ideal_alpha_phi, ideal_alpha_phi_stdev, "NA", "NA", ideal_2eq_3ax_4eq_psi, ideal_2eq_3ax_4eq_psi_stdev, anomeric_position, linkage_type ]

    # alpha sugar with LINKAGE_NA LinkageType
    elif anomeric_position == "ALPHA_LINKS" and linkage_type == "LINKAGE_NA":
        if verbose:
            print "Residue number:", sugar_num
            print "  ", anomeric_position, linkage_type
            print "   ideal phi major:", ideal_alpha_phi
            print "\tideal phi stdev major:", ideal_alpha_phi_stdev
            print "   ideal phi minor:", "NA", 
            print "\tideal phi stdev minor:", "NA"
            print "      actual phi:", input_pose.phi( sugar_num )
            print "   ideal psi:", "NA", 
            print "\tideal psi stdev:", "NA"
            print "      actual psi:", input_pose.psi( sugar_num )
        return [ ideal_alpha_phi, ideal_alpha_phi_stdev, "NA", "NA", "NA", "NA", anomeric_position, linkage_type ]

    # beta sugar with _2AX_3EQ_4AX_LINKS LinkageType
    elif anomeric_position == "BETA_LINKS" and linkage_type == "_2AX_3EQ_4AX_LINKS":
        if verbose:
            print "Residue number:", sugar_num
            print "  ", anomeric_position, linkage_type
            print "   ideal phi major:", ideal_beta_phi_major, 
            print "\tideal phi stdev major:", ideal_beta_phi_stdev_major
            print "   ideal phi minor:", ideal_beta_phi_minor, 
            print "\t\tideal phi stdev minor:", ideal_beta_phi_stdev_minor
            print "      actual phi:", input_pose.phi( sugar_num )
            print "   ideal psi:", ideal_2ax_3eq_4ax_psi, 
            print "\t\tideal psi stdev:", ideal_2ax_3eq_4ax_psi_stdev
            print "      actual psi:", input_pose.psi( sugar_num )
        return [ ideal_beta_phi_major, ideal_beta_phi_stdev_major, ideal_beta_phi_minor, ideal_beta_phi_stdev_minor, ideal_2ax_3eq_4ax_psi, ideal_2ax_3eq_4ax_psi_stdev, anomeric_position, linkage_type ]

    # beta sugar with _2EQ_3AX_4EQ_LINKS LinkageType
    elif anomeric_position == "BETA_LINKS" and linkage_type == "_2EQ_3AX_4EQ_LINKS":
        if verbose:
            print "Residue number:", sugar_num
            print "  ", anomeric_position, linkage_type
            print "   ideal phi major:", ideal_beta_phi_major, 
            print "\tideal phi stdev major:", ideal_beta_phi_stdev_major
            print "   ideal phi minor:", ideal_beta_phi_minor, 
            print "\t\tideal phi stdev minor:", ideal_beta_phi_stdev_minor
            print "      actual phi:", input_pose.phi( sugar_num )
            print "   ideal psi:", ideal_2eq_3ax_4eq_psi, 
            print "\t\tideal psi stdev:", ideal_2eq_3ax_4eq_psi_stdev
            print "      actual psi:", input_pose.psi( sugar_num )
        return [ ideal_beta_phi_major, ideal_beta_phi_stdev_major, ideal_beta_phi_minor, ideal_beta_phi_stdev_minor, ideal_2eq_3ax_4eq_psi, ideal_2eq_3ax_4eq_psi_stdev, anomeric_position, linkage_type ]

    # beta sugar with LINKAGE_NA LinkageType
    elif anomeric_position == "BETA_LINKS" and linkage_type == "LINKAGE_NA":
        if verbose:
            print "Residue number:", sugar_num
            print "  ", anomeric_position, linkage_type
            print "   ideal phi major:", ideal_beta_phi_major, 
            print "\tideal phi stdev major:", ideal_beta_phi_stdev_major
            print "   ideal phi minor:", ideal_beta_phi_minor, 
            print "\t\tideal phi stdev minor:", ideal_beta_phi_stdev_minor
            print "      actual phi:", input_pose.phi( sugar_num )
            print "   ideal psi:", "NA",
            print "\tideal psi stdev:", "NA"
            print "      actual psi:", input_pose.psi( sugar_num )
        return [ ideal_beta_phi_major, ideal_beta_phi_stdev_major, ideal_beta_phi_minor, ideal_beta_phi_stdev_minor, "NA", "NA", anomeric_position, linkage_type ]

    # otherwise, there is no information for this carbohydrate anomeric/LinkageType information
    else:
        return [ "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA" ]



def set_glycan_to_ideal_SugarBB_phi_psi( sugar_nums, input_pose, verbose = False ):
    """
    Given a list of <sugar_nums>, reset the phi and psi of each residue to the ideal phi and psi determined by the SugarBB statistical data (the sugar_bb scoring term implements the CHI energy function using this data)
    Only changes phi and psi if statistical data is available for the given anomeric position and/or linkage type
    Can give this function a single residue number as well, instead of a list of residue numbers
    :param sugar_nums: list( carbohydrate residue Pose numbers ) or int( single carbohydrate residue Pose number )
    :param input_pose: Pose
    :param verbose: bool( do you want to print out the phi and psi reset information for each residue? ) Default = False
    :return: Pose
    """
    # imports
    from rosetta.core.id import phi_dihedral, psi_dihedral, omega_dihedral
    from rosetta.core.pose.carbohydrates import set_glycosidic_torsion


    # ensure that either a list of ints was passed, or a single int was passed
    # if the input argument is a list, ensure it is a list of valid ints
    if type( sugar_nums ) == list:
        for sugar_num in sugar_nums:
            if not type( sugar_num ) == int:
                print "\nArgument error. You gave me a list for <sugar_nums>, but it does not only contain integers. Returning input Pose."
                return input_pose
            # ensure the residue number is found within the <input_pose> as well
            if not sugar_num in range( 1, input_pose.n_residue() + 1 ):
                print "\nArgument error. You gave me a list of residue numbers where one or more numbers is not actually within the Pose range. Returning input Pose."
                return input_pose
            # also ensure it is a carbohydrate residue
            if not input_pose.residue( sugar_num ).is_carbohydrate():
                print "\nArgument error. You gave me a list of residue numbers where one or more numbers is not actually a carbohydrate. Returning input Pose."
                return input_pose                
    else:
        # if it's not a list, ensure it is a valid int
        if not type( sugar_nums ) == int:
            print "\nArgument error. You gave didn't give me a list or an int for <sugar_nums>. Returning input Pose."
            return input_pose
        # if an int was passed, ensure it is found within <input_pose>
        else:
            if not sugar_nums in range( 1, input_pose.n_residue() + 1 ):
                print "\nArgument error. You gave me a single residue number, but it is not actually within the Pose range. Returning input Pose."
                return input_pose
            # also ensure it is a carbohydrate residue
            if not input_pose.residue( sugar_nums ).is_carbohydrate():
                print "\nArgument error. You gave me a single residue number, but it is not actually a carbohydrate. Returning input Pose."
                return input_pose                

    # copy over the input_pose
    pose = input_pose.clone()

    # if a single residue was given, put it into a list
    if type( sugar_nums ) == int:
        sugar_nums = [ sugar_nums ]

    # use get_ideal_SugarBB_phi_psi_info to set sugar residue(s) to statistical ideal phi/psi
    if type( sugar_nums ) == list:
        for sugar_num in sugar_nums:
            # get ideal data
            # order of data in list: [ ideal_phi_major, ideal_phi_stdev_major, ideal_phi_minor, ideal_phi_stdev_minor, ideal_psi, ideal_psi_stdev, anomeric_position, linkage_type ]
            # returns "NA" for value if no statistical data is available
            ideal_data_list = get_ideal_SugarBB_phi_psi_info( sugar_num, pose )

            # reset phi and psi, if statistical data is available
            # ideal phi major (never NA, a sugar is always alpha or beta)
            ideal_phi_major = ideal_data_list[ 0 ]
            set_glycosidic_torsion( phi_dihedral, pose, sugar_num, ideal_phi_major )

            # ideal psi (could be NA for sugars whose connection is not on linkage number 2, 3, or 4
            ideal_psi = ideal_data_list[ 4 ]
            if ideal_psi != "NA":
                set_glycosidic_torsion( psi_dihedral, pose, sugar_num, ideal_psi )

            # different verbose possibilities based on statistical data available
            if verbose:
                if ideal_psi != "NA":
                    print "Residue", sugar_num
                    print "   old phi:", input_pose.phi( sugar_num ), "new phi:", ideal_phi_major
                    print "   old psi:", input_pose.psi( sugar_num ), "new psi:", ideal_psi
                else:
                    print "Residue", sugar_num
                    print "   old phi:", input_pose.phi( sugar_num ), "new phi:", ideal_phi_major
                    print "   old psi:", input_pose.psi( sugar_num ), "no new psi; no statistical data available for residue type"        

    return pose



def get_fa_scorefxn_with_given_weights( weights_dict, verbose = False ):
    """
    Return an sf from get_fa_scoretype but with adjusted weights <scoretypes> with given <weights>
    If <input_scoretype> is not already part of the <sf>, this function will add it to <sf> with a weight of <weight>, and then get the score
    Will exit if the string( <input_scoretype> ) is not a valid ScoreType
    :param weights_dict: dict( ScoreType or str of ScoreType name : int( or float( weight ) ) )
    :param verbose: bool( print the final weights of the returned ScoreFunction? ) Default = False
    "return: ScoreFunction( fa_scorefxn with adjusted weights of given scoretypes )
    """
    from rosetta import get_fa_scorefxn
    from rosetta.core.scoring import score_type_from_name


    # argument check - check the passed argument is a dict
    if not isinstance( weights_dict, dict ):
        print "You didn't give me a dictionary for your input. I need a dict of ScoreType (or name) : weight. Exiting."
        sys.exit()
        
    # get a standard fa_scorefxn to start with
    sf = get_fa_scorefxn()
    
    # make a dummy ScoreType to use for isinstance() checking
    fa_dun = score_type_from_name( "fa_dun" )

    # for each entry of the dictionary, change the weight
    for scoretype_name in weights_dict:
        # if the key is a string
        if isinstance( scoretype_name, str ):
            try:
                scoretype = score_type_from_name( scoretype_name )
            except:
                print "\nThe string name: '%s' does not appear to be a valid ScoreType. Exiting" %scoretype_name
                sys.exit()
            
            # get the corresponding weight
            weight = weights_dict[ scoretype_name ]

            # set the weight
            sf.set_weight( scoretype, weight )

        # if the argument is a ScoreType object
        elif isinstance( scoretype_name, type( fa_dun ) ):
            # adjust the weight in the scorefxn using the corresponding weight given
            sf.set_weight( scoretype_name, weights_dict[ scoretype_name ] )

        # else, I don't know what they gave me as a scoretype
        else: 
            print "I'm not sure what '%s' is from your ScoreType key in your <weights_dict> argument. Exiting" %scoretype_name
            sys.exit()

    # if verbose, print out the weights of the new ScoreFunction
    if verbose:
        print "\nNew score weights sf:\n%s\n" %( "\n".join( [ "%s: %s" %( str( name ), sf.get_weight( name ) ) for name in sf.get_nonzero_weighted_scoretypes() ] ) )

    # return the newly weighted fa_scorefxn
    return sf



def make_fa_scorefxn_from_file( scorefxn_file, verbose = False ):
    """
    Return an sf from get_fa_scoretype but with adjusted weights read from the passed <scorefxn_file>
    If a ScoreType is not already part of the sf, this function will add it to sf with the specified weight
    :param scorefxn_file: str( path/to/scorefxn.sf file with ScoreTypes and weights )
    :param verbose: bool( print the final weights of the returned ScoreFunction? ) Default = False
    "return: ScoreFunction( fa_scorefxn with adjusted weights of given scoretypes )
    """
    from rosetta import get_fa_scorefxn
    from rosetta.core.scoring import score_type_from_name


    # make a standard fa_scorefxn
    sf = get_fa_scorefxn()

    # edit the ScoreFunction
    try:
        with open( scorefxn_file, "rb" ) as fh:
            score_lines = fh.readlines()
    except:
        print "Something seems to be wrong with your scorefxn_file ( %s ). Please check your input" %scorefxn_file
        sys.exit()
    for score_line in score_lines:
        # ignore new line characters
        score_line = score_line.rstrip()

        # get the score type and weight from each line, which should be space delimited
        if score_line != '':
            # ignore commented lines
            if not score_line.startswith( '#' ):
                # check if the user followed the right format
                try:
                    score_line_split = score_line.split( ' ' )
                    score_type = str( score_line_split[0] )

                    # if user wants a multiplier
                    if score_line_split[1] == '*':
                        score_weight = sf.get_weight( score_type_from_name( score_type ) ) * float( score_line_split[2] )

                    # else user wants a specific value
                    else:
                        score_weight = float( score_line_split[1] )
                except:
                    print "\nIt seems that your scorefxn_file did not follow the proper format. Please follow 'ScoreType Weight' or 'ScoreType * Multiplier'"
                    sys.exit()
        try:
            sf.set_weight( score_type_from_name( score_type ), score_weight )
        except:
            print "It seems like you did not pass a valid ScoreType (or some other issue). Check your scorefxn_file (%s)" %input_args.scorefxn_file
            sys.exit()

    # if verbose, print the final weights
    if verbose:
        print "\nScore weights used in this sf:\n%s\n" %( "\n".join( [ "%s: %s" %( str( name ), sf.get_weight( name ) ) for name in sf.get_nonzero_weighted_scoretypes() ] ) )

    return sf




def show_score_breakdown( sf, pose ):
    """
    Shows the breakdown of the <pose>'s total score by printing the score of each nonzero weighted ScoreType in <sf>
    :param sf: ScoreFunction
    :param pose: Pose
    """
    # print out each score
    print "\n".join( [ "%s: %s" %( score_type, round( sf.score_by_scoretype( pose, score_type ), 3 ) ) for score_type in sf.get_nonzero_weighted_scoretypes() ] )
    print



def show_score_breakdown_by_res( sf, res_num, pose ):
    """
    Shows the breakdown of the <pose>'s total score by printing the score for the <res_num> of each nonzero weighted ScoreType in <sf>
    :param sf: ScoreFunction
    :param res_num: int( Pose residue number )
    :param pose: Pose
    """
    # get the string version of the residue's energies from Pose
    sf( pose )
    res_energies_obj = pose.energies().residue_total_energies( res_num )
    res_energies = res_energies_obj.show_nonzero().strip().split( ' ' )

    # print out each score
    for ii in range( 1, len( res_energies ) + 1, 2 ):
        print res_energies[ ii - 1 ], res_energies[ ii ]



def show_sf_weights_breakdown( scorefxn ):
    """
    Shows the breakdown of nonzeroweighted ScoreTypes in <scorefxn>
    :param sf: ScoreFunction
    :param pose: Pose
    """
    # print out each score
    print
    print "\n".join( [ "%s: %s" %( score_type, scorefxn.get_weight( score_type ) ) for score_type in scorefxn.get_nonzero_weighted_scoretypes() ] )
    print



def get_score_by_scoretype( sf, input_scoretype, pose, weight = 1.0, verbose = False ):
    """
    Return the specified <input_scoretype> value using <sf> on the <pose>
    If <input_scoretype> is not already part of the <sf>, this function will add it to <sf> with a weight of <weight>, and then get the score
    Will exit if the string( <input_scoretype> ) is not a valid ScoreType
    :param sf: ScoreFunction
    :param input_scoretype: str( name of a valid ScoreType ) or ScoreType object
    :param pose: Pose
    :param weight: int( or float( weight of the desired <input_scoretype> if it's not already a part of <sf> ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( <input_scoretype> score of <pose> )
    """
    from rosetta.core.scoring import score_type_from_name


    if verbose:
        print "Getting the total score of a specific ScoreType for the pose"

    # make a dummy ScoreType to use for isinstance() checking
    fa_dun = score_type_from_name( "fa_dun" )

    # check to make sure the string <input_scoretype> is a valid ScoreType
    if isinstance( input_scoretype, str ):
        try:
            scoretype = score_type_from_name( input_scoretype )
        except:
            print
            print
            print "The string name:", input_scoretype, "does not appear to be a valid ScoreType. Exiting"
            sys.exit()

    # if the input was a ScoreType, it must be valid. Turn its name into a string version for use
    elif isinstance( input_scoretype, type( fa_dun ) ):
        scoretype = input_scoretype
        input_scoretype = str( input_scoretype )

    # else, I don't know what they gave me
    else:
        print "I'm not sure what '%s'is. Exiting" %scoretype
        sys.exit()

    # tell user it was successful in turning str( <input_scoretype> ) to a ScoreType object
    if verbose:
        print "Successfully turned", input_scoretype, "into a valid ScoreType object"

    # check to see if <sf> already has the given <input_scoretype>
    has_scoretype = False
    current_scoretypes = sf.get_nonzero_weighted_scoretypes()
    for score_type in current_scoretypes:
        if str( score_type ) == input_scoretype:
            has_scoretype = True
            break

    # if has_scoretype is False, add the <input_scoretype> to the <sf> with weight <weight>
    if not has_scoretype:
        print "It appears the ScoreType is not already included in your ScoreFunction, adding", input_scoretype, "with a weight of", weight
        sf.set_weight( scoretype, weight )

    # score the <pose> with the specified <input_scoretype> and return the value
    score = sf.score_by_scoretype( pose, scoretype )

    return score



def get_residue_score_by_scoretype( sf, input_scoretype, seq_pos, pose, weight = 1.0, verbose = False ):
    """
    Given a <pose>, a specific <seq_pos> for a residue, a <input_scoretype> of type string or ScoreType, return the score of that residue using the ScoreFunction <sf>
    If the <sf> does not already have <input_scoretype> as a non-zero weighted ScoreType, then this function will add it with a weight of <weight> and then score the residue
    Will exit and print an error message if the string <input_scoretype> is invalid
    :param sf: ScoreType
    :param input_scoretype: str( a valid ScoreType name ) or a ScoreType object
    :param seq_pos: int( the residue of interest )
    :param pose: Pose
    :param weight: int( or float( the weight to be used if <input_scoretype> is not already in <sf> ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( the <input_scoretype> score of residue <seq_pos> )
    """
    from rosetta.core.scoring import score_type_from_name


    if verbose:
        print "Getting and comparing each residue between the two poses based on ScoreType"

    # make a dummy ScoreType to use for isinstance() checking
    fa_dun = score_type_from_name( "fa_dun" )

    # check to make sure the string <input_scoretype> is a valid ScoreType
    if isinstance( input_scoretype, str ):
        try:
            scoretype = score_type_from_name( input_scoretype )
        except:
            print
            print
            print "The string name:", input_scoretype, "does not appear to be a valid ScoreType. Exiting"
            sys.exit()

    # if the input was a ScoreType, it must be valid. Turn its name into a string version for use
    elif isinstance( input_scoretype, type( fa_dun ) ):
        scoretype = input_scoretype
        input_scoretype = str( input_scoretype )

    # else, I don't know what they gave me
    else:
        print "I'm not sure what", input_scoretype, "is. Exiting"
        sys.exit()

    # tell user it was successful going from str( <input_scoretype> ) to ScoreType
    if verbose:
        print "Successfully turned", input_scoretype, "into a valid ScoreType object"

    # check to see if <sf> already has the given <input_scoretype>
    has_scoretype = False
    current_scoretypes = sf.get_nonzero_weighted_scoretypes()
    for score_type in current_scoretypes:
        if str( score_type ) == input_scoretype:
            has_scoretype = True
            break

    # if has_scoretype is False, add the <input_scoretype> to the <sf> with weight <weight>
    if not has_scoretype:
        if verbose:
            print "It appears the ScoreType is not already included in your ScoreFunction, adding", input_scoretype, "it with a weight of", weight
        sf.set_weight( scoretype, weight )

    # get access to <pose>'s energies() object
    sf( pose )

    # get the <input_scoretype> score of the residue at <seq_pos>
    score_of_res_of_scoretype = pose.energies().residue_total_energies( seq_pos ).get( scoretype )
    
    
    # return the score of the residue given the ScoreType given
    return score_of_res_of_scoretype



def get_scoretype_with_biggest_score_diff( in_decoy, in_native, in_sf ):
    """
    Determine which ScoreType in <in_sf> corresponds to the largest total score difference between <in_decoy> and <in_native>. Returns an object that holds information for the biggest positive difference and the biggest negative difference.
    :param in_decoy: Pose 1
    :param in_native: Pose 2
    :param in_sf: ScoreFunction
    :return: obj( holds ScoreType name and the delta total score of that ScoreType )
    """
    # make clones of decoy, native, and sf
    decoy = in_decoy.clone()
    native = in_native.clone()
    sf = in_sf.clone()

    # score the poses to gain access to the energies object
    sf( decoy )
    sf( native )

    # holders for the biggest delta score and ScoreType
    biggest_pos_delta_score_by_scoretype = None
    biggest_pos_delta_score_type = None
    biggest_neg_delta_score_by_scoretype = None
    biggest_neg_delta_score_type = None

    # use the energies object to determine which ScoreType gives the biggest difference
    for score_type in sf.get_nonzero_weighted_scoretypes():
        # get the total scores by ScoreType
        decoy_score_by_scoretype = decoy.energies().total_energies().get( score_type )
        native_score_by_scoretype = native.energies().total_energies().get( score_type )
        delta_score_by_scoretype = decoy_score_by_scoretype - native_score_by_scoretype

        # update the prelimary data holders
        if biggest_pos_delta_score_by_scoretype is None:
            # if this holder is None, then the other ones are as well, so update them all with the first residue
            biggest_pos_delta_score_by_scoretype = delta_score_by_scoretype
            biggest_pos_delta_score_scoretype = score_type
            biggest_neg_delta_score_by_scoretype = delta_score_by_scoretype
            biggest_neg_delta_score_scoretype = score_type
        # otherwise, if data has been found, check to see if the new score is higher and update the holders
        else:
            if delta_score_by_scoretype > biggest_pos_delta_score_by_scoretype:
                biggest_pos_delta_score_by_scoretype = delta_score_by_scoretype
                biggest_pos_delta_score_scoretype = score_type
            if delta_score_by_scoretype < biggest_neg_delta_score_by_scoretype:
                biggest_neg_delta_score_by_scoretype = delta_score_by_scoretype
                biggest_neg_delta_score_scoretype = score_type

    # make an object to hold relevant data
    data = DataHolder()
    data.biggest_pos_delta_score_by_scoretype = biggest_pos_delta_score_by_scoretype
    data.score_type_pos = biggest_pos_delta_score_scoretype
    data.str_score_type_pos = str( biggest_pos_delta_score_scoretype )
    data.biggest_neg_delta_score_by_scoretype = biggest_neg_delta_score_by_scoretype
    data.score_type_neg = biggest_neg_delta_score_scoretype
    data.str_score_type_neg = str( biggest_neg_delta_score_scoretype )

    return data



def get_res_with_biggest_score_diff( in_decoy, in_native, in_sf, decoy_to_native_res_map = None ):
    """
    Iterate through <in_decoy> residues and compare each total energy to the corresponding residue in <in_native>. If the residues don't match up exactly in order from 1 - n, give me a <decoy_to_native_res_map> dictionary.
    If multiple residues have the same biggest delta score, then only the last residue seen with this value will be stored
    :param in_decoy: Pose 1
    :param in_native: Pose 2
    :param in_sf: ScoreFunction
    :param decoy_to_native_res_map: dict( if you're comparing a decoy to a native that has a different numbering scheme, pass in a dictionary that maps each resnum from the decoy to its corresponding resnum in the native )
    :return: obj( holds residue <in_decoy> number and score, <in_native> number and score, delta total score, corresponding ScoreType, delta ScoreType score )
    """
    # make clones of decoy, native, and sf
    decoy = in_decoy.clone()
    native = in_native.clone()
    sf = in_sf.clone()

    # score the poses to gain access to the energies object
    sf( decoy )
    sf( native )

    # holders for the biggest delta score and corresponding residue numbers
    biggest_delta_score = None
    biggest_delta_score_decoy_num = None
    biggest_delta_score_native_num = None

    # iterate through each residue in decoy and compare its total score to the corresponding residue found in native
    for decoy_num in range( 1, decoy.n_residue() + 1 ):
        # get the corresponding residue number in native, if needed
        if decoy_to_native_res_map is not None:
            native_num = decoy_to_native_res_map[ decoy_num ]
        else:
            native_num = decoy_num

        # get the total score of each residue and the delta
        decoy_score = decoy.energies().residue_total_energy( decoy_num )
        native_score = native.energies().residue_total_energy( native_num )
        delta_score = decoy_score - native_score

        # update the prelimary data holders
        if biggest_delta_score is None:
            # if this holder is None, then the other two are as well, so update them all with the first residue
            biggest_delta_score = delta_score
            biggest_delta_score_decoy_num = decoy_num
            biggest_delta_score_native_num = native_num
        # otherwise, if data has been found, check to see if the new score is higher and update the holders
        else:
            if delta_score > biggest_delta_score:
                biggest_delta_score = delta_score
                biggest_delta_score_decoy_num = decoy_num
                biggest_delta_score_native_num = native_num

    # make holders for the info related to the specific ScoreType
    biggest_delta_score_by_scoretype = None
    biggest_delta_score_scoretype = None

    # using the decoy and native residue numbers with the biggest delta total score, find the corresponding biggest ScoreType
    for score_type in sf.get_nonzero_weighted_scoretypes():
        # get the residue scores by ScoreType
        decoy_score_by_scoretype = decoy.energies().residue_total_energies( biggest_delta_score_decoy_num ).get( score_type )
        native_score_by_scoretype = native.energies().residue_total_energies( biggest_delta_score_native_num ).get( score_type )
        delta_score_by_scoretype = decoy_score_by_scoretype - native_score_by_scoretype

        # update the prelimary data holders
        if biggest_delta_score_by_scoretype is None:
            # if this holder is None, then the other one is well, so update them all with the first residue
            biggest_delta_score_by_scoretype = delta_score_by_scoretype
            biggest_delta_score_scoretype = score_type
        # otherwise, if data has been found, check to see if the new score is higher and update the holders
        else:
            if delta_score_by_scoretype > biggest_delta_score_by_scoretype:
                biggest_delta_score_by_scoretype = delta_score_by_scoretype
                biggest_delta_score_scoretype = score_type
                
    # make an object to hold relevant data
    data = DataHolder()
    data.biggest_delta_score = biggest_delta_score
    data.decoy_num = biggest_delta_score_decoy_num
    data.native_num = biggest_delta_score_native_num
    data.native_res_info = native.pdb_info().pose2pdb( biggest_delta_score_native_num ).strip()
    data.biggest_delta_score_by_scoretype = biggest_delta_score_by_scoretype
    data.score_type = biggest_delta_score_scoretype
    data.str_score_type = str( biggest_delta_score_scoretype )

    return data



def get_sugar_bb_only_sf( weight = 1 ):
    """
    Creates and returns a ScoreFunction with only the sugar_bb term as a non-zero. Can specify weight with <weight> param
    :param weight: int( or float( weight of the sugar_bb term in the sf ) ). Default = 1
    :return: ScoreFunction
    """
    from rosetta import get_fa_scorefxn
    from rosetta.core.scoring import score_type_from_name


    # instantiate a fa ScoreFunction
    sugar_sf = get_fa_scorefxn()
    
    # get a list of all the ScoreTypes in the fa_scorefxn
    scoretypes = sugar_sf.get_nonzero_weighted_scoretypes()

    # set all ScoreTypes to 0
    for scoretype in scoretypes:
        sugar_sf.set_weight( scoretype, 0 )
            
    # set sugar_bb ScoreType to given <weight>
    # doing this outside the loop in case "sugar_bb" stops becoming a part of fa_scorefxn
    sugar_bb = score_type_from_name( "sugar_bb" )
    sugar_sf.set_weight( sugar_bb, weight )
    
    
    # return the sugar_sf that now only has the sugar_bb ScoreType as non-zero
    return sugar_sf



def apply_sugar_constraints_to_sf( sf, pose, weight = 1.0, verbose = False ):
    """
    Applies bond distance and bond angle constraints to sugar branch points within <pose> with a weight of <weight>
    Uses a HarmonicFun for the bond distance, and a CircularHarmonicFun for the bond angle
    :param sf: ScoreFunction
    :param pose: Pose
    :param weight: int( or float( what do you want the weight to be for the bond distance and bond angle constraints? ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: ScoreFunction( sf including sugar constraints )
    """
    from rosetta import AtomID
    from rosetta.core.scoring.constraints import AtomPairConstraint, AngleConstraint
    from rosetta.core.scoring import score_type_from_name
    from rosetta.core.scoring.func import HarmonicFunc, CircularHarmonicFunc


    # get list of chemical edges ( marked as -2 in FoldTree f )
    edges = pose.fold_tree().get_chemical_edges()

    if verbose:
        print "Constraining the bond distance and bond angles for", len( edges ), "sugar residues"

    # need atom IDs from residue numbers
    for ii in range( 1, len( edges ) + 1 ):
        # extract info from edge
        start_res_num = edges[ii].start()
        start_atom_name = edges[ii].start_atom()
        stop_res_num = edges[ii].stop()
        stop_atom_name = edges[ii].stop_atom()

        # get atom id's of upper and lower branch point atoms
        upper_branch_atom_index = pose.residue( start_res_num ).atom_index( start_atom_name )
        lower_branch_atom_index = pose.residue( stop_res_num ).atom_index( stop_atom_name )

        # get atom index numbers of the first adjacent heavy atoms to each branch point atom
        upper_branch_angle_atom_index = pose.residue( start_res_num ).first_adjacent_heavy_atom( upper_branch_atom_index )
        lower_branch_angle_atom_index = pose.residue( stop_res_num ).first_adjacent_heavy_atom( lower_branch_atom_index )

        # make AtomID objects
        upper_branch_atom = AtomID( upper_branch_atom_index, start_res_num )
        lower_branch_atom = AtomID( lower_branch_atom_index, stop_res_num )
        upper_branch_angle_atom = AtomID( upper_branch_angle_atom_index, start_res_num )
        lower_branch_angle_atom = AtomID( lower_branch_angle_atom_index, stop_res_num )

        # get bond length of the atoms to constrain
        dist_upper_lower = pose.conformation().bond_length( upper_branch_atom, lower_branch_atom )

        # make and apply harmonic function
        harmonic_func = HarmonicFunc( dist_upper_lower, 0.5 )
        atm_pair_constraint = AtomPairConstraint( upper_branch_atom, lower_branch_atom, harmonic_func )
        pose.add_constraint( atm_pair_constraint )

        # get bond angle of the atoms to constrain protein C, N - C sugar
        ang_upper_angle_upper_lower = pose.conformation().bond_angle( upper_branch_angle_atom, upper_branch_atom, lower_branch_atom )

        # make and apply circular harmonic function
        circ_harmonic_func = CircularHarmonicFunc( ang_upper_angle_upper_lower, 5 )
        ang_constraint = AngleConstraint( upper_branch_angle_atom, upper_branch_atom, lower_branch_atom, circ_harmonic_func )
        pose.add_constraint( ang_constraint )

        # get bond angle of the atoms to constrain protein N - C, O sugar
        ang_upper_lower_lower_angle = pose.conformation().bond_angle( upper_branch_atom, lower_branch_atom, lower_branch_angle_atom )

        # make and apply circular harmonic function
        circ_harmonic_func = CircularHarmonicFunc( ang_upper_lower_lower_angle, 5 )
        ang_constraint = AngleConstraint( upper_branch_atom, lower_branch_atom, lower_branch_angle_atom, circ_harmonic_func )
        pose.add_constraint( ang_constraint )

    # add constraints to score function weights
    sf.set_weight( score_type_from_name( "atom_pair_constraint" ), weight )
    sf.set_weight( score_type_from_name( "angle_constraint" ), weight )

    return sf



def create_random_AtomPair_cst_file( residues_to_be_constrained, residues_to_be_constrained_against, pose, num_of_atoms_per_res, cst_filename, dist_max = 10.0, stdev = 1.0, tolerance = 0.5 ):
    """
    Constrain <num_of_atoms_per_res> heavy atoms from each residue from <residues_to_be_constrained> to a heavy atom from <residues_to_be_constrained_against> in <pose>
    THIS COULD ENTER AN INFINITE LOOP. TODO: NEED REVISIT THE LOGIC OF THIS CODE WHEN I MOVE ON FROM 3AY4
    Since this will be random, the only thing ensured is that the atoms are not H or virtual, and that the same atom-to-atom cst is not repeated. Otherwise, a single atom could have multiple constraints.
    Standard function call writes out a FLAT_HARMOINIC function with a <stdev> of 1.0 and a <tolerance> of 0.5
    :param residues_to_be_constrained: list( residues whose atom(s) should be constrained )
    :param residues_to_be_constrained_against: list( residues to whom the atom(s) should be constrained against )
    :param pose: Pose
    :param num_of_atoms_per_res: int( how many atoms per residue should be constrained )
    :param cst_filename: str( /path/to/cst filename you desire ) Default behavior is to add .cst if not already there
    :param dist_max: float( the max distance between the two atoms in the constraint ). Default = 10.0
    :param stdev: float( standard deviation for the FLAT_HARMONIC function )
    :param tolerance: float( tolerance for the FLAT_HARMONIC function )
    :return: bool( True if successful, False if not )
    """
    # imports
    from random import choice


    # change the input var names to something shorter
    cst_residues = residues_to_be_constrained
    cst_residues.sort()
    cst_options = residues_to_be_constrained_against
    cst_options.sort()
    natoms = num_of_atoms_per_res

    # constrain as many atoms as specified by user to atoms found in residues_to_be_constrained_against
    # if no atom-to-atom distances satisfy <dist_max>, then that residue just isn't being constrained
    atom_pair_cst_lines = []
    for cst_res_num in cst_residues:
        cst_res = pose.residue( cst_res_num ) 

        # create empty lists to keep info for each atom that passes <dist_max> filter
        # dictionary makes more sense, but the logic is difficult. List positions stay the same anyway
        cst_atom_num_list = []
        to_cst_atom_num_list = []
        to_cst_res_num_list = []
        atom_pair_dist_list = []
        for cst_atom_num in range( 1, cst_res.natoms() + 1 ):
            cst_atom = cst_res.atom( cst_atom_num )

            # if this is a non-hydrogen, non-virtual atom
            if (not cst_res.atom_is_hydrogen( cst_atom_num )) and (not cst_res.is_virtual( cst_atom_num )):
                # collect atom indices of atoms within <dist_max>
                for to_cst_res_num in cst_options:
                    to_cst_res = pose.residue( to_cst_res_num )

                    # check every heavy atom within the other residues given
                    for to_cst_atom_num in range( 1, to_cst_res.natoms() + 1 ):
                        # get the to_cst_atom info
                        to_cst_atom = to_cst_res.atom( to_cst_atom_num )

                        # if this is a non-hydrogen, non-virtual atom
                        if (not to_cst_res.atom_is_hydrogen( to_cst_atom_num )) and (not to_cst_res.is_virtual( to_cst_atom_num )):
                            # if this atom is within <dist_max> of the cst_atom
                            atom_pair_dist = cst_atom.xyz().distance( to_cst_atom.xyz() )
                            if atom_pair_dist <= dist_max:
                                # update the lists with this atom-pair cst info
                                cst_atom_num_list.append( cst_atom_num )
                                to_cst_atom_num_list.append( to_cst_atom_num )
                                to_cst_res_num_list.append( to_cst_res_num )
                                atom_pair_dist_list.append( str( round( atom_pair_dist, 2 ) ) )

        # grab the index numbers of the atom-pairs to be constrained
        atom_pair_cst_index_numbers = []

        # if there are fewer atoms than atom-pair constraints desired, but not zero, then add them all to the AtomPair cst file
        if ( len( cst_atom_num_list ) != 0 ) and (len( cst_atom_num_list ) < natoms):
            atom_pair_cst_index_numbers = range( len( cst_atom_num_list ) )
        # else, pick as many random index numbers as desired
        elif ( len( cst_atom_num_list ) != 0 ) and (len( cst_atom_num_list ) >= natoms):
            index_options = range( len( cst_atom_num_list ) )
            for ii in range( natoms ):
                index_choice = choice( index_options )
                atom_pair_cst_index_numbers.append( index_choice )
                index_options.remove( index_choice )
        # otherwise, this residue has no atoms passing the <dist_max> filter, so skip it
        else:
            pass

        # create the AtomPair constraints using the index numbers
        for index_num in atom_pair_cst_index_numbers:
            # pull out the needed information by index number
            cst_atom_num = cst_atom_num_list[ index_num ]
            to_cst_atom_num = to_cst_atom_num_list[ index_num ]
            to_cst_res_num = to_cst_res_num_list[ index_num ]

            # get the corresponding atom names
            cst_atom_name = cst_res.atom_name( cst_atom_num ).strip()
            to_cst_atom_name = pose.residue( to_cst_res_num ).atom_name( to_cst_atom_num ).strip()

            # get the PDB numbering info the tag in the cst file
            # this splits on the space between the PDB num and the chain
            pdb_name = pose.pdb_info().pose2pdb( cst_res_num ).strip().split( ' ' )
            cst_res_pdb_name = pdb_name[1] + pdb_name[0]
            pdb_name = pose.pdb_info().pose2pdb( to_cst_res_num ).strip().split( ' ' )
            to_cst_res_pdb_name = pdb_name[1] + pdb_name[0]
            tag = "_to_".join( [ cst_res_pdb_name, to_cst_res_pdb_name ] )

            # create the AtomPair cst line for the cst file
            atom_pair_cst_line = "AtomPair %s %s %s %s FLAT_HARMONIC %s %s %s %s\n" %( cst_atom_name, 
                                                                                       str( cst_res_num ), 
                                                                                       to_cst_atom_name, 
                                                                                       str( to_cst_res_num ), 
                                                                                       str( atom_pair_dist_list[ index_num ] ), 
                                                                                       str( stdev ), 
                                                                                       str( tolerance ), 
                                                                                       tag )
            atom_pair_cst_lines.append( atom_pair_cst_line )

    # add .cst to the end of the cst_filename, if not already there
    if not cst_filename.endswith( ".cst" ):
        cst_filename += ".cst"

    # write out the new cst file
    try:
        with open( cst_filename, "wb" ) as fh:
            fh.writelines( atom_pair_cst_lines )
        return True
    except:
        return False
        


def SugarSmallMover( seqpos, input_pose, angle_max, set_phi = True, set_psi = True, set_omega = True ):
    """
    Randomly resets the phi, psi, and omega values of the sugar residue <seqpos> in <pose> to old_value +/- angle_max/2
    Emulates the SmallMover but with the additional omega mover
    :param seqpos: int( the pose number for the residue )
    :param input_pose: Pose
    :param angle_max: int( or float( the max angle around which the phi/psi/omega could move ) )
    :param set_phi: bool( do you want to change the phi angle? ) Default = True
    :param set_psi: bool( do you want to change the psi angle? ) Default = True
    :param set_omega: bool( do you want to change the omega angle? ) Default = True
    :return: Pose
    """
    # imports
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg
    from rosetta.core.id import phi_dihedral, psi_dihedral, omega_dihedral
    from rosetta.core.pose.carbohydrates import set_glycosidic_torsion


    # copy the input pose
    pose = input_pose.clone()

    # from rosetta.protocols.simple_moves:BackboneMover.cc file for SmallMover
    big_angle = angle_max
    small_angle = big_angle / 2.0

    # get current phi, psi, and omega
    old_phi = pose.phi( seqpos )
    old_psi = pose.psi( seqpos )
    old_omega = pose.omega( seqpos )

    # get random values for phi, psi, and omega
    # this specific format is pulled from rosetta.protocols.simple_moves:ShearMover::make_move
    new_phi = periodic_range( old_phi - small_angle + rg().uniform() * big_angle, 360.0 )
    new_psi = periodic_range( old_psi - small_angle + rg().uniform() * big_angle, 360.0 )
    new_omega = periodic_range( old_omega - small_angle + rg().uniform() * big_angle, 360.0 )

    # set the new values
    if set_phi:
        set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
    if set_psi:
        set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
    if set_omega:
        set_glycosidic_torsion( omega_dihedral, pose, seqpos, new_omega )

    return pose



def SugarShearMover_3ay4( seqpos, input_pose, angle_max ):
    """
    HARD-CODED TO WORK WITH 3AY4
    Emulates the ShearMover but with an additional omega consideration
    :param seqpos: int( the pose number for the residue )
    :param input_pose: Pose
    :param angle_max: int( or float( the max angle around which the phi/psi/omega could move ) )
    :return: Pose
    """
    # imports
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg, random_range
    from rosetta.core.id import phi_dihedral, psi_dihedral, omega_dihedral, \
        omega2_dihedral, MainchainTorsionType
    from rosetta.core.pose.carbohydrates import find_seqpos_of_saccharides_parent_residue, \
        get_reference_atoms, get_glycosidic_torsion, set_glycosidic_torsion


    # copy the input pose and get the residue object of interest
    pose = input_pose.clone()
    residue = pose.residue( seqpos )

    # from rosetta.protocols.simple_moves:BackboneMover.cc file for SmallMover
    big_angle = angle_max
    small_angle = big_angle / 2.0

    # get the residue properties
    residue_properties = pose.residue_type( seqpos ).properties()

    # get the string versions of the VariantType(s) of this residue
    # this is a Vector1 list, so index starts at 1!
    variant_types = residue_properties.get_list_of_variants()

    # get current phi, psi, omega, and omega2 of the current residue
    old_phi = get_glycosidic_torsion( phi_dihedral, pose, seqpos )
    old_psi = get_glycosidic_torsion( psi_dihedral, pose, seqpos )
    old_omega = get_glycosidic_torsion( omega_dihedral, pose, seqpos )
    old_omega2 = get_glycosidic_torsion( omega2_dihedral, pose, seqpos )

    # get current phi, psi, omega, and omega2 of the following residue, if it exists
    # ie. check if this is not a tail residue by checking if it is not an UPPER_TERMINUS_VARIANT
    if not "UPPER_TERMINUS_VARIANT" in variant_types:
        following_seqpos = seqpos + 1
        old_phi_following = get_glycosidic_torsion( phi_dihedral, pose, following_seqpos )
        old_psi_following = get_glycosidic_torsion( psi_dihedral, pose, following_seqpos )
        old_omega_following = get_glycosidic_torsion( omega_dihedral, pose, following_seqpos )
        old_omega2_following = get_glycosidic_torsion( omega2_dihedral, pose, following_seqpos )

    # get the parent of this residue
    parent_seqpos = find_seqpos_of_saccharides_parent_residue( residue )
    parent = pose.residue( parent_seqpos )

    # get current phi, psi, omega, and omega2 of the parent residue
    old_phi_parent = get_glycosidic_torsion( phi_dihedral, pose, parent_seqpos )
    old_psi_parent = get_glycosidic_torsion( psi_dihedral, pose, parent_seqpos )
    old_omega_parent = get_glycosidic_torsion( omega_dihedral, pose, parent_seqpos )
    old_omega2_parent = get_glycosidic_torsion( omega2_dihedral, pose, parent_seqpos )

    # create the shear_delta to be used
    # this specific format is pulled from rosetta.protocols.simple_moves:ShearMover::make_move
    shear_delta = small_angle - rg().uniform() * big_angle

    # hard coding based on sequence position for 3ay4
    # D and F 2 and 4, E and G 2
    if seqpos == 217 or seqpos == 219 or seqpos == 222 or seqpos == 441 or seqpos == 443 or seqpos == 446:
        if random_range( 1, 2 ) == 1:
            # make and set the new phi value of the residue
            new_phi = periodic_range( old_phi - shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
            # make and set the negating new psi value of the following residue
            new_psi_following = periodic_range( old_psi_following + shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, following_seqpos, new_psi_following )
        else:
            # make and set the new psi value of the residue
            new_psi = periodic_range( old_psi - shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
            # make and set the negating new phi value of the following residue
            new_phi_following = periodic_range( old_phi_following + shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, following_seqpos, new_phi_following )
    # D and F 3, the branching Man
    elif seqpos == 218 or seqpos == 442: 
        if random_range( 1, 2 ) == 1:
            # make and set the new phi value of the residue
            new_phi = periodic_range( old_phi - shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
            # make and set the negating new psi value of the PARENT residue
            new_psi_parent = periodic_range( old_psi_parent + shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, parent_seqpos, new_psi_parent )
        else: 
            # make and set the new psi value of the residue
            new_psi = periodic_range( old_psi - shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
            # make and set the negating new phi value of the PARENT residue
            new_phi_parent = periodic_range( old_phi_parent + shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, parent_seqpos, new_phi_parent )
    # E and G 1, the residue with the omega connected to the branching Man
    elif seqpos == 221 or seqpos == 445:
        if random_range( 1, 2 ) == 1:
            # make and set the new phi value of the residue
            new_phi = periodic_range( old_phi - shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
            # make and set the negating new omega value of the same residue
            new_omega = periodic_range( old_omega + shear_delta, 360.0 )
            set_glycosidic_torsion( omega_dihedral, pose, seqpos, new_omega )
        else:
            # make and set the new psi value of the residue
            new_psi = periodic_range( old_psi - shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
            # make and set the negating new psi value of the following residue
            new_psi_following = periodic_range( old_psi_following + shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, following_seqpos, new_psi_following )
    # tail residues on chain D/E and F/G
    elif seqpos == 220 or seqpos == 223 or seqpos == 444 or seqpos == 447:
        if random_range( 1, 2 ) == 1:
            # make and set the new phi value of the residue
            new_phi = periodic_range( old_phi - shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
        else:
            # make and set the new psi value of the residue
            new_psi = periodic_range( old_psi - shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
    else:
        print "\nWhat residue did you give me?"
        sys.exit()

    return pose



def SugarShearMover( seqpos, input_pose, angle_max ):
    """    
    Emulates the ShearMover but with an additional omega consideration
    :param seqpos: int( the pose number for the residue )
    :param input_pose: Pose
    :param angle_max: int( or float( the max angle around which the phi/psi/omega could move ) )
    :return: Pose
    """
    # imports
    from rosetta.basic import periodic_range
    from rosetta.numeric.random import rg, random_range
    from rosetta.core.id import phi_dihedral, psi_dihedral, omega_dihedral, \
        omega2_dihedral, MainchainTorsionType
    from rosetta.core.pose.carbohydrates import find_seqpos_of_saccharides_parent_residue, \
        get_reference_atoms, get_glycosidic_torsion, set_glycosidic_torsion


    # copy the input pose and get the residue object of interest
    pose = input_pose.clone()
    residue = pose.residue( seqpos )

    # from rosetta.protocols.simple_moves:BackboneMover.cc file for SmallMover
    big_angle = angle_max
    small_angle = big_angle / 2.0

    # get the residue properties
    residue_properties = pose.residue_type( seqpos ).properties()

    # get the string versions of the VariantType(s) of this residue ( this is a Vector1 list, so index starts at 1 )
    variant_types = residue_properties.get_list_of_variants()

    # get current phi, psi, omega, and omega2
    old_phi = get_glycosidic_torsion( phi_dihedral, pose, seqpos )
    old_psi = get_glycosidic_torsion( psi_dihedral, pose, seqpos )
    old_omega = get_glycosidic_torsion( omega_dihedral, pose, seqpos )
    old_omega2 = get_glycosidic_torsion( omega2_dihedral, pose, seqpos )

    # get current phi, psi, omega, and omega2 of the following residue, if it exists
    # ie. check if this is not a tail residue by checking if it is not an UPPER_TERMINUS_VARIANT
    if not "UPPER_TERMINUS_VARIANT" in variant_types:
        following_seqpos = seqpos + 1
        old_phi_following = get_glycosidic_torsion( phi_dihedral, pose, following_seqpos )
        old_psi_following = get_glycosidic_torsion( psi_dihedral, pose, following_seqpos )
        old_omega_following = get_glycosidic_torsion( omega_dihedral, pose, following_seqpos )
        old_omega2_following = get_glycosidic_torsion( omega2_dihedral, pose, following_seqpos )

    # get the parent of this residue
    parent_seqpos = find_seqpos_of_saccharides_parent_residue( residue )
    parent = pose.residue( parent_seqpos )

    # get current phi, psi, omega, and omega2 of the parent residue
    old_phi_parent = get_glycosidic_torsion( phi_dihedral, pose, parent_seqpos )
    old_psi_parent = get_glycosidic_torsion( psi_dihedral, pose, parent_seqpos )
    old_omega_parent = get_glycosidic_torsion( omega_dihedral, pose, parent_seqpos )
    old_omega2_parent = get_glycosidic_torsion( omega2_dihedral, pose, parent_seqpos )

    # create the shear_delta to be used
    # this specific format is pulled from rosetta.protocols.simple_moves:ShearMover::make_move
    shear_delta = small_angle - rg().uniform() * big_angle

    # determine which torsions are available for sampling
    # to ensure generality, this takes the long way and stores each torsion this residue has by checking the reference atoms for that torsion
    torsions_available = []
    # for each torsion name possible ( phi, psi, omega 1, 2, and 3 ), check to see if this residue has atoms that define that torsion
    for torsion_name in MainchainTorsionType.names.keys():
        # skipping omega3 because it seems to not be implemented in Rosetta yet ( makes get_reference_atoms crash )
        if torsion_name != "omega3_dihedral":
            # get the corresponding TorsionType object
            torsion_type = MainchainTorsionType.names[ torsion_name ]
            # determine if there are any atoms that define this torsion
            # should always be either 4 atoms ( True, has torsion ) or 0 atoms ( False, does not have torsion ), therefore, can use bool( len )
            if bool( len( get_reference_atoms( torsion_type, pose, seqpos ) ) ):
                torsions_available.append( torsion_type )

    # determine which torsions are available in the parent residue
    # to ensure generality, this takes the long way and stores each torsion this residue has by checking the reference atoms for that torsion
    parent_torsions_available = []
    # for each torsion name possible ( phi, psi, omega 1, 2, and 3 ), check to see if this residue has atoms that define that torsion
    for torsion_name in MainchainTorsionType.names.keys():
        # skipping omega3 because it seems to not be implemented in Rosetta yet ( makes get_reference_atoms crash )
        if torsion_name != "omega3_dihedral":
            # get the corresponding TorsionType object
            torsion_type = MainchainTorsionType.names[ torsion_name ]
            # determine if there are any atoms that define this torsion
            # should always be either 4 atoms ( True, has torsion ) or 0 atoms ( False, does not have torsion ), therefore, can use bool( len )
            if bool( len( get_reference_atoms( torsion_type, pose, parent_seqpos ) ) ):
                parent_torsions_available.append( torsion_type )

    # randomly pick the torsion that will get sampled for this residue
    # using a Rosetta random module for consistency; random_range includes the ending value, so need to do -1
    torsion = torsions_available[ random_range( 0, len( torsions_available ) - 1 ) ]

    # check to see if the residue has an omega2 torsion, such as a core GlcNAc connected to a peptide
    if omega2_dihedral in torsions_available:
        # if the residue's phi torsion is changing, then change its own omega torsion by the inverse amount
        # or similarly, if the residue's omega is changing, then change its own phi torsion by the inverse amount
        if torsion == phi_dihedral or torsion == omega_dihedral:
            # make and set the new phi value of the residue
            new_phi = periodic_range( old_phi - shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
            # make and set the new omega negating value of the same residue
            new_omega = periodic_range( old_omega + shear_delta, 360.0 )
            set_glycosidic_torsion( omega_dihedral, pose, seqpos, new_omega )
        # if the residue's psi torsion is changing, then change its own omega2 torsion by the inverse amount
        # or similarly, if the residue's omega2 is changing, then change its own psi torsion by the inverse amount
        elif torsion == psi_dihedral or torsion == omega2_dihedral:
            # make and set the new psi value of the residue
            new_psi = periodic_range( old_psi - shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
            # make and set the new omega2 negating value of the same residue
            new_omega2 = periodic_range( old_omega2 + shear_delta, 360.0 )
            set_glycosidic_torsion( omega2_dihedral, pose, seqpos, new_omega2 )
        # otherwise, I don't know about this exception case yet, so exit
        else:
            print "\nThis residue pair presents me with an exception case! There's something I didn't think about... Exiting on case 1."
            sys.exit()
            
    # check to see if the residue only has an omega torsion, such as a residue with an exocyclic connection to a branch residue
    # since this is an elif statement checking for omega after checking for omega2, if the residue has both an omega and a omega2, 
    # it would not do both the omega2 move and the omega move, just the omega2 move
    elif omega_dihedral in torsions_available:
        # if the residue's phi torsion is changing, then change its own omega torsion by the inverse amount
        # or similarly, if the residue's omega is changing, then change its own phi torsion by the inverse amount
        if torsion == phi_dihedral or torsion == omega_dihedral:
            # make and set the new phi value of the residue
            new_phi = periodic_range( old_phi - shear_delta, 360.0 )
            set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
            # make and set the new omega negating value of the same residue
            new_omega = periodic_range( old_omega + shear_delta, 360.0 )
            set_glycosidic_torsion( omega_dihedral, pose, seqpos, new_omega )
        # if the residue's psi is changing, then change the following residues psi torsion by the inverse amount
        elif torsion == psi_dihedral:
            # make and set the new psi value of the residue
            new_psi = periodic_range( old_psi - shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
            # make and set the new psi negating value of the following residue
            new_psi_following = periodic_range( old_psi_following + shear_delta, 360.0 )
            set_glycosidic_torsion( psi_dihedral, pose, following_seqpos, new_psi_following )
        # otherwise, I don't know about this exception case yet, so exit
        else:
            print "\nThis residue pair presents me with an exception case! There's something I didn't think about... Exiting on case 2."
            sys.exit()
                
    # if this residue doesn't have an omega, then change this residue's phi and the following residue's ( n+1 ) psi
    # if no residue exists upstream of this residue ( ie. this residue is a tail residue ), then no negating change is needed
    else:
        # check to see if this is a UPPER_TERMINUS_VARIANT, ie. a tail residue
        if "UPPER_TERMINUS_VARIANT" in variant_types:
            # this is a tail residue, so just change its phi, psi, or omega torsion ( there is no following residue, so no negating change )
            if torsion == phi_dihedral:
                # make and set the new phi value
                new_phi = periodic_range( old_phi - shear_delta, 360.0 )
                set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
            elif torsion == psi_dihedral:
                # make and set the new psi value
                new_psi = periodic_range( old_psi - shear_delta, 360.0 )
                set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
            elif torsion == omega_dihedral:
                # make and set the new omega value
                new_omega = periodic_range( old_omega - shear_delta, 360.0 )
                set_glycosidic_torsion( omega_dihedral, pose, seqpos, new_omega )
            elif torsion == omega2_dihedral:
                # make and set the new omega2 value
                new_omega2 = periodic_range( old_omega2 - shear_delta, 360.0 )
                set_glycosidic_torsion( omega2_dihedral, pose, seqpos, new_omega2 )
            # otherwise, I don't know about this exception case yet, so exit
            else:
                print "\nThis residue pair presents me with an exception case! Likely it has an omega3 I didn't think about. Exiting on case 3."
                sys.exit()
                    
        # otherwise, this is not a tail residue, so change the torsion of the residue and the negating torsion of the following residue
        # or, if this residue is a branch point, change the torsion of the residue and the negating torsion of the PARENT residue
        else:
            # if this non-tail residue is a branch point
            if residue.is_branch_point():
                # check to see if the PARENT residue has an omega torsion, such as a residue with an exocyclic connection
                if omega_dihedral in parent_torsions_available:
                    # if we're changing phi, change omega of the PARENT residue by the opposite amount
                    if torsion == phi_dihedral:
                        # make and set the new phi value of the residue
                        new_phi = periodic_range( old_phi - shear_delta, 360.0 )
                        set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
                        # make and set the new omega negating value of the PARENT residue
                        new_omega_parent = periodic_range( old_omega_parent + shear_delta, 360.0 )
                        set_glycosidic_torsion( omega_dihedral, pose, parent_seqpos, new_omega_parent )
                    # if we're changing psi, change the psi of the PARENT residue by the opposite amount
                    elif torsion == psi_dihedral:
                        # make and set the new psi value of the residue
                        new_psi = periodic_range( old_psi - shear_delta, 360.0 )
                        set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
                        # make and set the new psi negating value of the PARENT residue
                        new_psi_parent = periodic_range( old_psi_parent + shear_delta, 360.0 )
                        set_glycosidic_torsion( psi_dihedral, pose, parent_seqpos, new_psi_parent )
                    # otherwise, I don't know about this exception case yet, so exit
                    else:
                        print "\nThis residue pair presents me with an exception case! There's something I didn't think about... Exiting on case 4."
                        sys.exit()
                # otherwise, there is no omega torsion in this branch residue
                # so change one the phi/psi of the residue and the psi/phi of the parent residue by the opposite amount
                else:
                    # if we're changing phi, change psi of the PARENT residue by the opposite amount
                    if torsion == phi_dihedral:
                        # make and set the new phi value of the residue
                        new_phi = periodic_range( old_phi - shear_delta, 360.0 )
                        set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
                        # make and set the new psi negating value of the PARENT residue
                        new_psi_parent = periodic_range( old_psi_parent + shear_delta, 360.0 )
                        set_glycosidic_torsion( psi_dihedral, pose, parent_seqpos, new_psi_parent )
                    # similarly, if we're changing psi, change phi of the PARENT residue by the opposite amount
                    elif torsion == psi_dihedral:
                        # make and set the new psi value of the residue
                        new_psi = periodic_range( old_psi - shear_delta, 360.0 )
                        set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
                        # make and set the new phi negating value of the PARENT residue
                        new_phi_parent = periodic_range( old_phi_parent + shear_delta, 360.0 )
                        set_glycosidic_torsion( phi_dihedral, pose, parent_seqpos, new_phi_parent )
                    # otherwise, I don't know about this exception case yet, so exit
                    else:
                        print "\nThis residue pair presents me with an exception case! There's something I didn't think about... Exiting on case 5."
                        sys.exit()

            # otherwise, if this residue is not a branch point, not a tail, and there is no omega
            else:
                # if we're changing phi, change psi of the following residue by the opposite amount
                if torsion == phi_dihedral:
                    # make and set the new phi value of the residue
                    new_phi = periodic_range( old_phi - shear_delta, 360.0 )
                    set_glycosidic_torsion( phi_dihedral, pose, seqpos, new_phi )
                    # make and set the new psi negating value of the following residue
                    new_psi_following = periodic_range( old_psi_following + shear_delta, 360.0 )
                    set_glycosidic_torsion( psi_dihedral, pose, following_seqpos, new_psi_following )
                # if we're changing psi, change phi of the following residue by the opposite amount
                elif torsion == psi_dihedral:
                    # make and set the new psi of the residue
                    new_psi = periodic_range( old_psi - shear_delta, 360.0 )
                    set_glycosidic_torsion( psi_dihedral, pose, seqpos, new_psi )
                    # make and set the new phi negating value of the following residue
                    new_phi_following = periodic_range( old_phi_following + shear_delta, 360.0 )
                    set_glycosidic_torsion( phi_dihedral, pose, following_seqpos, new_phi_following )
                # otherwise, I don't know about this exception case yet, so exit
                else:
                    print "\nThis residue pair presents me with an exception case! Likely it has an omega I didn't think about. Exiting on case 6."
                    sys.exit()

    return pose



def get_atom_nums_within_radius( seq_pos, atom_num, input_pose, radius ):
    """
    Get a dict of residue numbers : atom nums in <input_pose> that are within <radius> Ang around atom <atom_num> in residue <seq_pos>. Uses heavy, non-virtual atoms for determination.
    :param seq_pos: int( Pose number of residue of interest )
    :param atom_num: int( atom index number of atom of interest )
    :param input_pose: Pose
    :param radius: float( or int( residues around <seq_pos> within this cutoff distance ))
    :return: dict( keys = resiudes numbers, values = corresponding atom numbers within <radius> cutoff )
    """
    # container for the CA or the C1 xyz of each residue in <input_pose>
    CA_C1_xyz_of_pose = {}

    for res_num in range( 1, input_pose.n_residue() + 1 ):
        # skip the <seq_pos> if the user specified to do so
        if res_num == seq_pos and include_seq_pos == False:
            pass
        # add all xyz coordinates of each CA or C1 atom in <input_pose> residues
        else:
            if input_pose.residue( res_num ).is_carbohydrate():
                atom_num = input_pose.residue( res_num ).atom_index( "C1" )
            else:
                atom_num = input_pose.residue( res_num ).atom_index( "CA" )

            # store the coordinates of this atom in the dict
            CA_C1_xyz_of_pose[ res_num ] = input_pose.residue( res_num ).atom( atom_num ).xyz()

    # CA or C1 xyz of residue of interest
    if input_pose.residue( seq_pos ).is_carbohydrate():
        atom_num = input_pose.residue( seq_pos ).atom_index( "C1" )
    else:
        atom_num = input_pose.residue( seq_pos ).atom_index( "CA" )
    seq_pos_CA_C1_xyz = input_pose.residue( seq_pos ).atom( atom_num ).xyz()

    # container for residues within specified 3 * <radius> to <seq_pos>
    res_nums_inside_triple_radius = []

    # first collect all residue numbers of those whose CA/C1 distance to <seq_pos> CA/C1 is less than 3 * <radius>
    for res_num in CA_C1_xyz_of_pose.keys():
        CA_C1_res_xyz = CA_C1_xyz_of_pose[ res_num ]
        if CA_C1_res_xyz.distance( seq_pos_CA_C1_xyz ) <= ( 3*radius ):
            res_nums_inside_triple_radius.append( res_num )

    # collect the xyz of each heavy atom in <seq_pos>
    seq_pos_heavy_atom_xyzs = {}
    for seq_pos_atom_num in range( 1, len( input_pose.residue( seq_pos ).atoms() ) + 1 ):
        #if not input_pose.residue( seq_pos ).atom_is_hydrogen( seq_pos_atom_num ): this returns HO1 as a non-H atom
        if not input_pose.residue( seq_pos ).atom_is_hydrogen( seq_pos_atom_num ) and not input_pose.residue( seq_pos ).is_virtual( seq_pos_atom_num ):
            seq_pos_heavy_atom_xyzs[ seq_pos_atom_num ] = ( input_pose.residue( seq_pos ).atom( seq_pos_atom_num ).xyz() )

    # now iter through each heavy atom of residues within 3 * <radius> to <seq_pos>
    # and see if it's within <radius> of a heavy atom on <seq_pos>
    res_nums_inside_radius = []
    for res_num in res_nums_inside_triple_radius:
        # for each heavy atom in this residue in the <input_pose>
        for atom_num in range( 1, len( input_pose.residue( res_num ).atoms() ) + 1 ):

            # move on if this residue already has a heavy atom within the <radius>
            if res_num not in res_nums_inside_radius:

                # get the xyz coordinates for each heavy atom
                if not input_pose.residue( res_num ).atom_is_hydrogen( atom_num ) and not input_pose.residue( res_num ).is_virtual( atom_num ):
                    heavy_atom_xyz = input_pose.residue( res_num ).atom( atom_num ).xyz()

                    # check this heavy atom against every heavy atom in the <seq_pos>
                    for seq_pos_heavy_atom_num in seq_pos_heavy_atom_xyzs.keys():
                        seq_pos_heavy_atom_xyz = seq_pos_heavy_atom_xyzs[ seq_pos_heavy_atom_num ]
                        if seq_pos_heavy_atom_xyz.distance( heavy_atom_xyz ) <= radius:
                            res_nums_inside_radius.append( res_num )
                            break

    return res_nums_inside_radius



def make_pack_rotamers_mover( sf, input_pose, pack_branch_points = True, residue_range = None, use_pack_radius = False, pack_radius = PACK_RADIUS, verbose = False ):
    """
    Returns a standard pack_rotamers_mover restricted to repacking and allows for sampling of current residue conformations for the <pose>
    IMPORTANT: DON'T USE for a Pose you JUST mutated  --  it's not setup to handle mutations
    If you want a packer task for a specific set of residues, set <residue_range> to a list of the residue sequence positions of interest
    If you have one or more residues of interest where you want to pack within a radius around each residue, set <use_pack_radius> to True, give a <pack_radius> (or use default), and set <residue_range> to the residue sequence positions of interest
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param pack_branch_points: bool( allow packing at branch points? ). Default = True
    :param residue_range: list( of int( valid residue sequence positions ) ) or a single int()
    :param use_pack_radius: bool( Do you want to pack residues within a certain <pack_radius> around residues specified in <residue_range>? ). Default = False
    :param pack_radius: int( or float( the radius around which you want to pack additional residues around the residues from <residue_range>. MUST have <set use_pack_radius> to True to do this ). Default = PACK_RADIUS = 20.0 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: a pack_rotamers_mover object
    """
    # imports
    from rosetta import Pose, standard_packer_task, RotamerTrialsMover


    # copy a fresh pose
    pose = Pose()
    pose.assign( input_pose )

    # check to make sure <residue_range> is a list of valid residue numbers
    if residue_range is not None:
        if isinstance( residue_range, list ):
            # quick check to make sure the numbers are valid residue sequence positions
            for res_num in residue_range:
                if not 0 < res_num <= pose.n_residue():
                    print "\nIt seems that", res_num, "is not a valid sequence position in the Pose you gave me. Exiting"
                    sys.exit()
        # if they gave a single residue for <residue_range>, then that's okay as long as it is a valid residue sequence position
        elif isinstance( residue_range, int ):
            # if this isn't a valid residue number, print message and exit
            if residue_range not in range( 1, pose.n_residue() + 1 ):
                print "\nYou gave me a single residue number, but it is not a valid sequence position in the Pose you gave me. Check if", residue_range, "is what you intended to give me. Exiting"
                sys.exit()
            # otherwise, put the integer into a list by itself
            residue_range = [ residue_range ]
        else:
            # not sure what type of thing they gave me, print an error and exit
            print "\nI'm not sure what you gave me. <residue_range> is of type", type( residue_range ), "and I needed a list of integers or at least a single integer. Exiting"
            sys.exit()

    # ensure that <residue_range> is not empty if <use_pack_radius> is set to True
    if residue_range is None and use_pack_radius is True:
        print "\nYou want me to use a pack radius, but you didn't give me a residue range! Be sure to give me at least one residue if you intended to set <use_pack_radius> to True. Exiting"
        sys.exit()


    # argument checks passed, now speak to the user
    if verbose:
        print "Making a pack rotamers mover"

    # get a list of the residues in the protein
    protein_residues = range( 1, pose.n_residue() + 1 )
        
    # make the beginning packer task
    # this results in a packer task where each residue is allowed to be packed ( except for the first and last resiudes? )
    task = standard_packer_task( pose )
    task.or_include_current( True )
    task.restrict_to_repacking()

    # if <use_pack_radius> is True
    if use_pack_radius is True:
        # get all residues within the <pack_radius> including the residues within the list
        res_nums_inside_pack_radius = get_res_nums_within_radius_of_residue_list( residue_range, pose, 
                                                                                  radius = pack_radius, 
                                                                                  include_res_nums = True )

        # get all residues outside the <pack_radius> using a set difference
        res_nums_outside_pack_radius = list( set( protein_residues ) - set( res_nums_inside_pack_radius ) )

        # turn off repacking for the residues not in the <pack_radius>
        [ task.nonconst_residue_task( res_num ).prevent_repacking() for res_num in res_nums_outside_pack_radius ]
            
    # alternatively, if no <pack_radius> given, but a <residue_range> was given, prevent repacking for the other residues
    if use_pack_radius is False and residue_range is not None:
        # prevent repacking for every residue NOT given by the user through <residue_range>
        [ task.nonconst_residue_task( res_num ).prevent_repacking() for res_num in protein_residues if res_num not in residue_range ]

    # otherwise, all residues should be packed, which is the default behavior of instantiating a packer task
    # so don't need to prevent packing for any residues if residue_range was None

    # if specified, turn off packing for each branch point residue
    if pack_branch_points is False:
        if verbose:
            print "  Turning off packing for branch points"
        # prevent repacking for each residue that is a branch point
        [ task.nonconst_residue_task( res_num ).prevent_repacking() for res_num in protein_residues if pose.residue( res_num ).is_branch_point() ]
                
    # update the user on the status of the task
    if verbose:
        # task.being_packed returns True if that residue is being packed in the task
        num_packed = [ task.being_packed( res_num ) for res_num in protein_residues ].count( True )
        num_not_packed = pose.n_residue() - num_packed
        print "  Preventing", num_not_packed, "residues from repacking in this packer task"
        print "  Allowing", num_packed, "residues to be packed"

    # make the pack rotamers mover and return it
    pack_rotamers_mover = RotamerTrialsMover( sf, task )

    return pack_rotamers_mover



def make_min_mover( sf, input_pose, jumps = None, allow_sugar_chi = False, minimization_type = "dfpmin_strong_wolfe", verbose = False ):
    """
    Returns a min_mover object suitable for glycosylated proteins (turns off carbohydrate chi for now)
    IMPORTANT: If using a specific type of minimization, <minimization_type>, it DOES NOT check beforehand if the type you gave is valid, so any string actually will work. It will break when used
    See 'https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview' for minimization type options
    Prints error message and exits if there was a problem
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param jumps: list( Jump numbers of Jump(s) to be minimized ). Default = None (ie. all Jumps) (Give empty list for no Jumps)
    :param allow_sugar_chi: bool( allow the chi angles of sugars to be minimized ). Default = False
    :param minimization_type: str( the type of minimization you want to use ). Default = "dfpmin_strong_wolfe"
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: min_mover object
    """
    from rosetta import Pose, MoveMap, MinMover


    if verbose:
        print "Making a min mover"
        
    # copy a fresh pose
    pose = Pose()
    pose.assign( input_pose )

    # instantiate a MoveMap
    mm = MoveMap()

    # set all Jumps to True if <jumps> is None
    if jumps is None:
        if verbose:
            print "  Setting all Jumps to be minimized"
        for ii in range( 1, pose.num_jump() + 1 ):
            mm.set_jump( ii, True )

    # otherwise, set the user-specified Jumps to True
    else:
        # make sure <jumps> is of type list
        if isinstance( jumps, list ):
            for ii in jumps:
                if verbose:
                    print "  Setting only", jumps, "to be minimized"
                # ensure <jumps> only contains valid Jumps for the <pose>
                if 0 < ii <= pose.num_jump():
                    mm.set_jump( ii, True )
                else:
                    print "\n%s doesn't appear to be a valid Jump for your pose, exiting" %ii
                    sys.exit()
        else:
            print "\nYou didn't give me a list for the jumps param - I need a list of valid jump numbers. Exiting"
            sys.exit()

    # set backbone angles to be minimized for both protein and sugar residues
    mm.set_bb( True )

    # turn on chi angle minimization for all amino acids
    # turn on or off chi angle minimization for sugars based on user input
    if allow_sugar_chi:
        if verbose:
            print "  Turning on chi angle mobility for all residues"
        [ mm.set_chi( residue.seqpos(), True ) for residue in pose ]
    else:
        if verbose:
            print "  Turning on chi angle mobility only for protein residues"
        [ mm.set_chi( residue.seqpos(), True ) for residue in pose if not residue.is_carbohydrate() ]

    # create a MinMover with options
    min_mover = MinMover( mm, sf, minimization_type, 0.01, True )

    return min_mover



def make_movemap_for_loop( loop, allow_bb_movement = True, allow_chi_movement = True, verbose = False ):
    """
    Given a loop, creates and returns a MoveMap allowing for bb and chi angle freedom, with a +/- 1 buffer
    :param loop: Loop( loop of interest )
    :param allow_bb_movement: bool( Do you want to allow BackBone movement for the loop? ). Default = True
    :param allow_chi_movement: bool( Do you want to allow chi angle movement for the loop? ). Default = True
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for Loop
    """
    from rosetta import MoveMap


    if verbose:
        print "Making a MoveMap for a Loop"
    
    # instantiate a MoveMap
    mm = MoveMap()

    # set only the Loop area to true
    # if BB is True
    if allow_bb_movement:
        mm.set_bb_true_range( loop.start() + 1, loop.stop() - 1 )

    # if chi is True
    if allow_chi_movement:
        mm.set_chi_true_range( loop.start() + 1, loop.stop() - 1 )

    return mm



def make_movemap_for_jumps( JUMP_NUMbers, verbose = False ):
    """
    Given a jump number, creates and returns a MoveMap allowing for movement of only the Jump(s) given in the list
    :param JUMP_NUM: list( int( valid Jump number ) )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for Jump(s)
    """
    from rosetta import MoveMap


    # instantiate a MoveMap
    mm = MoveMap()

    # print different verbose statements and make the MoveMap differently depending on if it's a single number, or a list
    if isinstance( JUMP_NUMbers, int ):
        # this is only one Jump
        if verbose:
            print "Making a MoveMap for Jump number", JUMP_NUMbers

        # set only the given Jump to True
        mm.set_jump( JUMP_NUMbers, True )

    elif isinstance( JUMP_NUMbers, list ):
        # this is a list of Jump numbers
        if verbose:
            print "Making a MoveMap for Jump numbers", JUMP_NUMbers

        # set only the given Jumps to True
        for JUMP_NUM in JUMP_NUMbers:
            mm.set_jump( JUMP_NUM, True )

    else:
        # I don't know what they gave me  -  I need an integer or a list
        print
        print "I'm not sure what ", JUMP_NUMbers, "is. I need an integer or a list of integers. Exiting"
        sys.exit()

    return mm



def make_movemap_for_sugars( pose, allow_bb_movement = True, allow_chi_movement = False, verbose = False ):
    """
    Given a <pose> object, return a MoveMap allowing for bb and/or chi movement of only sugar residues.
    :param pose: Pose( a Pose object )
    :param allow_bb_movement: bool( Do you want to allow sugarback bone movement? ). Default = True
    :param allow_chi_movement: bool( Do you want to allow sugarchi angle movement? ). Default = False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for all sugar residues in <pose>
    """
    from rosetta import MoveMap


    if verbose:
        print "Making a MoveMap for a Loop"
    
    # instantiate a MoveMap
    mm = MoveMap()

    ## set movement for only sugar residues
    # if BB is True
    if allow_bb_movement:
        for res in pose:
            if res.is_carbohydrate():
                mm.set_bb( res.seqpos(), True )

    # if chi is True
    if allow_chi_movement:
        for res in pose:
            if res.is_carbohydrate():
                mm.set_chi( res.seqpos(), True )
    
    # return the sugar MoveMap
    return mm



def make_movemap_for_range( seqpos_list, allow_bb_movement = True, allow_chi_movement = True, verbose = False ):
    """
    Given a list of <seqpos_list>, return a MoveMap allowing for bb and/or chi movement of only the passed residues.
    :param seqpos_list: list( int( a list of sequence positions for your pose ) )
    :param allow_bb_movement: bool( Do you want to allow sugarback bone movement? ). Default = True
    :param allow_chi_movement: bool( Do you want to allow sugarchi angle movement? ). Default = False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: MoveMap for all sugar residues in <pose>
    """
    from rosetta import MoveMap


    # check to make sure <seqpos_list> is a list
    if not isinstance( seqpos_list, list ):
        print "You didn't pass me a list for your <seqpos_list> argument. That's what I need to make your MoveMap. Exiting."
        sys.exit()
        
    if verbose:
        print "Making a MoveMap for the following residues:", seqpos_list
    
    # instantiate a MoveMap
    mm = MoveMap()

    # if BB is True
    if allow_bb_movement:
        for num in seqpos_list:
            mm.set_bb( num, True )

    # if chi is True
    if allow_chi_movement:
        for num in seqpos_list:
            mm.set_chi( num, True )
    
    # return the MoveMap
    return mm



# TODO: Do I actually need this function?
def do_min_with_this_mm( mm, sf, pose, minimization_type = "dfpmin_strong_wolfe", verbose = False ):
    """
    Minimizes a given Pose using dfpmin_strong_wolfe and the user-supplied ScoreFunction <sf> and MoveMap <mm>
    IMPORTANT: If using a specific type of minimization, <minimization_type>, it DOES NOT check beforehand if the type you gave is valid, so any string actually will work. It will break when used
    See 'https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview' for minimization type options
    :param mm: finished MoveMap
    :param sf: ScoreFunction
    :param pose:  Pose
    :param minimization_type: str( the type of minimization you want to use ). Default = "dfpmin_strong_wolfe"
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: minimized Pose
    """
    from rosetta import MinMover


    # create a MinMover with options
    min_mover = MinMover( mm, sf, minimization_type, 0.01, True )

    if verbose:
        print "Minimizing the Pose with the following min_mover"
        print min_mover

    # apply the min_mover
    min_mover.apply( pose )

    return pose



def do_pack_min( sf, input_pose, residue_range = None, jumps = None, pack_branch_points = True, use_pack_radius = False, pack_radius = PACK_RADIUS, allow_sugar_chi = False, minimization_type = "dfpmin_strong_wolfe", verbose = False, pmm = None ):
    """
    Makes and applies a packer task and basic min mover to <pose> using the supplied ScoreFunction <sf>
    IMPORTANT: If using a specific type of minimization, <minimization_type>, it DOES NOT check beforehand if the type you gave is valid, so any string actually will work. It will break when used
    See 'https://www.rosettacommons.org/docs/latest/rosetta_basics/structural_concepts/minimization-overview' for minimization type options
    min_mover will set SUGAR BB freedom only - not CHI
    Will minimize over all Jumps, unless <jumps> is set to a list of Jump numbers to minimize instead
    If you want a packer task for a specific set of residues, set <residue_range> to a list of the residue sequence positions of interest
    If you have one or more residues of interest where you want to pack within a certain radius around each residue, set <use_pack_radius> to True, give a <pack_radius> (or use default value), and set <residue_range> to the residue sequence positions of interest
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param residue_range: list( of int( valid residue sequence positions ) )
    :param jumps = list( valid Jump numbers to be minimized if not all Jumps should be minimized). Default = None (ie. all jumps minimized) (Give empty list for no Jumps)
    :param pack_branch_points: bool( allow packing at branch points? ). Default = True
    :param use_pack_radius: bool( Do you want to pack residues within a certain <pack_radius> around residues specified in <residue_range>? ). Default = False
    :param pack_radius: int( or float( the radius around which you want to pack additional residues around the residues from <residue_range>. MUST have <set use_pack_radius> to True to do this ). Default = PACK_RADIUS = 20.0 Angstroms
    :param allow_sugar_chi: bool( allow chi minimization for sugars? ). Default = False
    :param minimization_type: str( the type of minimization you want to use ). Default = "dfpmin_strong_wolfe"
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: packed and minimized Pose
    """
    from rosetta import Pose


    if verbose:
        print "Packing and minimizing the Pose"
        
    # get a fresh copy of the input_pose
    pose = Pose()
    pose.assign( input_pose )
        
    # apply the PyMOL_Mover if passed
    if pmm is not None:
        # used in the function elsewhere to ensure the program doesn't use a broken mover
        pmm_worked = False
        try:
            pmm.apply( pose )
            pmm_worked = True
        except:
            print "Something was wrong with your PyMOL_Mover -- continuing without watching"
            pass

    # make and apply the pack_rotamers_mover
    pack_rotamers_mover = make_pack_rotamers_mover( sf, 
                                                    pose, 
                                                    pack_branch_points = pack_branch_points, 
                                                    residue_range = residue_range, 
                                                    use_pack_radius = use_pack_radius, 
                                                    pack_radius = pack_radius, 
                                                    verbose = verbose )
    pack_rotamers_mover.apply( pose )
    if pmm is not None and pmm_worked:
        pmm.apply( pose )

    # make and apply the min_mover
    min_mover = make_min_mover( sf, 
                                pose, 
                                jumps = jumps, 
                                allow_sugar_chi = allow_sugar_chi,
                                minimization_type = minimization_type, 
                                verbose = verbose )
    min_mover.apply( pose )
    if pmm is not None and pmm_worked:
        pmm.apply( pose )

    return pose



# used in functions within Protocols code block
def ramp_score_weight( sf, score_type_str, target_weight, current_step, total_steps, verbose = False ):
    """
    Ramps up or down a given ScoreType toward its target score. Uses percentage of completion to calculate amount of steps
    Current and Total steps -- Say you're doing 100 rounds, if you're on round 57, current_step = 57, total_steps = 100
    :param sf: ScoreFunction
    :param score_type_str: str( valid ScoreType )
    :param target_weight: int( or float( value of your target ScoreType weight )
    :param current_step: int( current step of how many rounds you're doing )
    :param total_steps: int( number of rounds you're doing )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: the ScoreFunction with the adjusted weight
    """
    from rosetta.core.scoring import score_type_from_name


    # get the ScoreType object given the str value passed in
    try:
        score_type = score_type_from_name( score_type_str )
    except RuntimeError:
        print
        print
        print score_type_str, "is not a valid score type, exiting"
        sys.exit()
    except:
        print
        print
        print "I'm not sure what happened. Check your input to ramp_score_weight. Exiting"
        sys.exit()

    # need to be able to do this a minimum of 10 times, otherwise the math will break
    if total_steps < 10:
        print
        print "You need to run a loop at least more than 10 times - otherwise the math in here won't work"
        sys.exit()

    # adjust the total_steps so that it's actually 10% less than what it actually is
    # this is so that the final 10% of the simulation will run with the target score weights
    total_steps = int( total_steps * 0.9 )

    # get the current weight of the score type
    current_weight = sf.get_weight( score_type )

    # make appropriate weight adjustments
    if current_weight > target_weight:  # we're decreasing the weight
        amount_left = target_weight - current_weight
        moves_left = float( total_steps - current_step )
        add_weight = amount_left / moves_left  # should be negative

        if verbose:
            print "Changing", score_type_str, "from", current_weight, "to", current_weight + add_weight

        # set the new weight
        sf.set_weight( score_type, current_weight + add_weight )

    elif target_weight > current_weight:  # we're increasing the weight
        amount_left = target_weight - current_weight
        moves_left = float( total_steps - current_step )
        add_weight = amount_left / moves_left  # should be positive

        if verbose:
            print "Changing", score_type_str, "from", current_weight, "to", current_weight + add_weight

        # set the new weight
        sf.set_weight( score_type, current_weight + add_weight )

    else:
        # we're at the target weight, so keep it there
        if verbose:
            print score_type_str, "is currently at the target weight of", target_weight

        # set to target weight
        sf.set_weight( score_type, target_weight )

    return sf



def CCD_loop_closure( loop, pose ):
    """
    Closes the given <loop> in <pose> using CCD
    :param loop: Loop( loop of interest )
    :param pose: Pose
    :return: Pose( <pose> with closed <loop> )
    """
    from rosetta.protocols.loops.loop_closure.ccd import CCDLoopClosureMover


    # make a MoveMap allowing only the loop to be flexible
    mm = make_movemap_for_loop( loop )

    # close the loop using CCD
    ccd = CCDLoopClosureMover( loop, mm )
    ccd.apply( pose )

    return pose



def calc_avg_of_list( in_list ):
    """
    Returns the average of a <in_list> of values, or None if the list is empty
    :param in_list: list( a list of integers or floats )
    :return: float( average of the list )
    """
    try:
        if len( in_list ) != 0:
            avg = float( sum( in_list ) / len( in_list ) )
        else:
            avg = None
    except:
        avg = None

    return avg



def calc_distance( vec1, vec2 ):
    """
    Calculates the distances between two points in 3D space
    :param vec1: list( x, y, z coordinates ) of point 1
    :param vec2: list( x, y, z coordinates ) of point 2
    :return: float( distance )
    """
    from math import sqrt, pow

    # assert that two lists were passed for coordinates
    if not isinstance( vec1, list ):
        print
        print "vec1 isn't a list - give me a list of xyz coordinates. Exiting"
        sys.exit()
    if not isinstance( vec2, list ):
        print
        print "vec2 isn't a list - give me a list of xyz coordinates. Exiting"
        sys.exit()

    # takes in two lists of xyz coordinates
    x1 = vec1[0]
    y1 = vec1[1]
    z1 = vec1[2]
    x2 = vec2[0]
    y2 = vec2[1]
    z2 = vec2[2]

    # calculate the distance
    dist = sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) + pow( z2 - z1, 2 ) )

    return dist



def get_res_nums_within_radius( res_num_in, input_pose, radius, include_res_num = False ):
    """
    Use the nbr_atom_xyz to find residue numbers within <radius> of <pose_num> in <pose>
    The nbr_atom seems to be C4 on carbohydrates
    :param res_num_in: int( Pose residue number )
    :param input_pose: Pose
    :param radius: int or float( radius around <pose_num> to use to select resiudes )
    :param include_res_num: bool( do you want to include <res_num> in the return list? ) Default = False
    :return: list( Pose residue numbers within <radius> of <pose_num>
    """
    # clone the <input_pose>
    pose = input_pose.clone()

    # container for the centers of each residue in pose
    centers_of_res = []

    # fill up the centers container
    for res_num in range( 1, pose.n_residue() + 1 ):
        center = pose.residue( res_num ).nbr_atom_xyz()
        centers_of_res.append( center )

    # container for residues inside the <radius>
    res_nums_in_radius = []

    # nbr_xyz of the residue of interest
    res_num_xyz = pose.residue( res_num_in ).nbr_atom_xyz()

    for res_num in range( 1, pose.n_residue() + 1 ):
        # this will get the xyz of the residue of interest, but it will be removed from the final list if desired
        # (since it will be added as 0 will always be less than <radius>)
        # get the center of the residue
        center = pose.residue( res_num ).nbr_atom_xyz()

        # keep the residue number if the nbr_atom_xyz is less than <radius>
        if center.distance( res_num_xyz ) <= radius:
            res_nums_in_radius.append( res_num )

    # if the user didn't want the residue of interest in the return list, remove it
    if not include_res_num:
        res_nums_in_radius.remove( res_num_in )

    return res_nums_in_radius



def get_res_nums_within_radius_of_residue_list( residues, input_pose, radius, include_res_nums = False ):
    """
    Find all residue numbers around the list of <residues> given in <input_pose> within <radius> Angstroms.
    Set <include_residues> if you want to include the list of passed <residues> in the return list of residue numbers.
    Uses the nbr_atom to calculate distance. The nbr_atom seems to be C4 on carbohydrates
    :param residues: list( Pose residue numbers )
    :param input_pose: Pose
    :param radius: int() or float( radius in Angstroms )
    :param include_res_nums: bool( do you want to include the passed <residues> in the return list of resiude numbers? ) Default = False
    :return: list( residues around passed <residues> list within <radius> Angstroms
    """
    # argument check: ensure passed <residues> argument is a list
    if type( residues ) != list:
        print "\nArgument error. You're supposed to past me a list of residue numbers for the <residues> argument. Returning None."
        return None

    # use get_res_nums_within_radius to get all residue numbers
    residues_within_radius = []
    for res_num in residues:
        residues_within_radius.extend( get_res_nums_within_radius( res_num, input_pose, radius, include_res_num = include_res_nums ) )


    # get the set of the list and sort the residue numbers
    set_of_residues_within_radius = [ res for res in set( residues_within_radius ) ]

    # it is possible that there are still residues from <residues> in the list, so remove them one by one if desired
    if not include_res_nums:
        for res in residues:
            try:
                set_of_residues_within_radius.remove( res )
            except ValueError:
                pass

    # sort
    set_of_residues_within_radius.sort()

    return set_of_residues_within_radius



def compare_single_residue_energies_between_poses( sf, res_num, native, decoy ):
    """
    Shows the breakdown of the <pose1>'s score by printing the score for the <res_num> of each nonzero weighted ScoreType in <sf> in <pose1> and <pose2>
    :param sf: ScoreFunction
    :param res_num: int( Pose residue number )
    :param native: Pose
    :param decoy: Pose
    """
    # imports
    from rosetta import score_type_from_name


    # get the string version of the residue's energies from Pose
    # native
    sf( native )
    res_energies_obj1 = native.energies().residue_total_energies( res_num )
    res_energies1 = res_energies_obj1.show_nonzero().strip().split()
    score_terms_in_native = []
    for term in res_energies1:
        try:
            # if this can't be converted into a float, it is a string. Therefore it is a score term
            float( term )
        except ValueError:
            score_terms_in_native.append( term.replace( ':', '' ) )

    # decoy
    sf( decoy )
    res_energies_obj2 = decoy.energies().residue_total_energies( res_num )
    res_energies2 = res_energies_obj2.show_nonzero().strip().split()
    score_terms_in_decoy = []
    for term in res_energies2:
        try:
            float( term )
        except ValueError:
            score_terms_in_decoy.append( term.replace( ':', '' ) )

    # get the score terms that are in common between the two poses
    score_terms = list( set( score_terms_in_native ) & set( score_terms_in_decoy ) )

    # print out each score of the residue from native vs decoy
    print "{0:<13} {1:^7} {2:^7} {3:>7}".format( "ScoreType", "Native", "Decoy", "Delta" )
    for score_term in score_terms:
        e1 = res_energies_obj1.get( score_type_from_name( score_term ) )
        e2 = res_energies_obj2.get( score_type_from_name( score_term ) )
        if score_term == "total_score":
            print "{0:13} {1:^7.3f} {2:7.3f} {3:7.3f}*".format( score_term, e1, e2, e2 - e1 )
        else:
            print "{0:13} {1:^7.3f} {2:7.3f} {3:7.3f}".format( score_term, e1, e2, e2 - e1 )



def compare_residue_energies_between_poses_with_cutoff( sf, pose1, pose2, diff_cutoff = 1.0, detailed_analysis = False, compare_using_this_scoretype = None, weight = 1.0, verbose = False ):
    """
    Uses the ScoreFunction <sf> to compare the individual energy per residue of <pose1> and returns the count of residues that have a difference of greater than +/- 1 of <diff_cutoff> compared to the corresponding residue in <pose2>
    If you are collecting this information for use in a detailed analysis, consider setting <detailed_analysis> to True and this function will return the residue objects of <pose2> that satisfied the energy <diff_cutoff>
    If <compare_using_this_scoretype> is not already a ScoreType that is part of the <sf>, this function will add it to <sf> with a weight of <weight>
    :param sf: ScoreFunction
    :param pose1: Pose( original pose )
    :param pose2: Pose( new pose )
    :param diff_cutoff: int( or float( energy cutoff in REUs ). Positive number. Default = 1.0 REU
    :param detailed_analysis: bool( changes return values if you want the residues that satisfied the energy <diff_cutoff>). Default = False
    :param compare_using_this_scoretype: str( what specific ScoreType do you want to use for comparison? ). Default = None = total_energy
    :param weight: int( or float( weight of the desired <str_scoretype> if it's not already a part of <sf> ). Default = 1.0
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: detailed_analysis = False, int( count of residues that satisfied the energy <diff_cutoff> )
    :return: detailed_analysis = True, int( count of residues that satisfied the energy <diff_cutoff> ), list( residues that satisfied the energy <diff_cutoff> of type Residue() )
    """
    # ensure that a number was given as a cutoff value
    try:
        diff_cutoff = abs( diff_cutoff )
        neg_diff_cutoff = diff_cutoff * -1
    except TypeError:
        print
        print
        print "Did you mean to give me a number for the diff_cutoff? Exiting"
        sys.exit()
    except:
        print
        print
        print "Something I don't understand happened. Raising the error."
        raise

    # ensure that the poses have the same number of amino acids (sequence identity doesn't matter)
    if pose1.total_residue() != pose2.total_residue():
        print "Your poses don't have the same number of residues, this comparison function won't work for this case. Exiting"
        print pose1
        print "vs"
        print pose2
        sys.exit()

    # otherwise continue on
    else:
        if verbose:
            if compare_using_this_scoretype is not None:
                print "Comparing residue v. residue energy using the specific ScoreType", compare_using_this_scoretype, "with a difference in energy cutoff value of", diff_cutoff
            else:
                print "Comparing residue v. residue total energy with a difference in energy cutoff value of", diff_cutoff

        # getting score only to allow access to the Pose's energies object
        sf( pose1 )
        sf( pose2 )

        # instantiate list that will hold Residue objects of residues that satisfied the energy <diff_cutoff>
        pose2_residues_with_greater_than_cutoff = []

        # number of residues in <pose2> that satisfied the energy <diff_cutoff>
        count = 0
        for ii in range( 1, pose1.total_residue() + 1 ):
            # if a ScoreType was given to use for comparison, calculate that energy
            if compare_using_this_scoretype is not None:
                # get ScoreType energy of the two residues
                pose1_res_E = get_residue_score_by_scoretype( sf, compare_using_this_scoretype, ii, pose1, weight = weight )
                pose2_res_E = get_residue_score_by_scoretype( sf, compare_using_this_scoretype, ii, pose2, weight = weight )

            # else, calculate total_energy and use that for comparison
            else:
                # get the total energy of each residue
                pose1_res_E = pose1.energies().residue_total_energy( ii )
                pose2_res_E = pose2.energies().residue_total_energy( ii )

            # get the energy difference, whether from total_energy or from a specific ScoreType
            E_diff = pose2_res_E - pose1_res_E

            # if E_diff is outside the +/- diff_cutoff
            if not neg_diff_cutoff < E_diff < diff_cutoff:
                pose2_residues_with_greater_than_cutoff.append( pose2.residue( ii ) )
                count += 1
                ##print "Pose1", pose1_res_E, "vs Pose2", pose2_res_E, "at position", ii

        if detailed_analysis:
            return count, pose2_residues_with_greater_than_cutoff
        else:
            return count



def compare_these_poses_by_score( sf, pose1, pose2, compare_using_this_scoretype = None, verbose = False ):
    """
    Calculates the energy of <pose1> and compares it to the energy of <pose2>, returning the Pose with the lower score
    :param sf: ScoreFunction
    :param pose1: Pose( pose1 )
    :param pose2: Pose( pose2 )
    :param compare_using_this_scoretype: str( what specific ScoreType do you want to use for comparison? ). Default = None = total_energy
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: the Pose with the lower energy of the two inputs
    """
    # return the better scored structure based on total_energy
    if compare_using_this_scoretype is None:
        if verbose:
            print "Returning the better scored Pose based on total energy"

        # return whichever Pose has the lower total energy score
        if sf( pose1 ) < sf( pose2 ):
            return pose1
        else:
            return pose2

    # otherwise, score by the given ScoreType and return the better scored Pose
    else:
        if verbose:
            print "Returning the better scored Pose based on", compare_using_this_scoretype, "energy"

        # get the specific ScoreType energy for the Poses
        pose1_E = get_score_by_scoretype( sf, compare_using_this_scoretype, pose1 )
        pose2_E = get_score_by_scoretype( sf, compare_using_this_scoretype, pose2 )

        # check to make sure neither energy is None (which means there was an error). If it is, exit this function
        if pose1_E is None or pose1_E is None:
            # don't need to print an error message because get_score_by_scoretype already does that
            sys.exit()

        # else return the better scored pose based on ScoreType
        else:
            if pose1_E < pose2_E:
                return pose1
            else:
                return pose2



# TODO-see how many inner and outer trials are actually necessary to hit a decent minimum
def get_best_structure_based_on_score( sf, pose, outer_trials = 3, inner_trials = 3, compare_using_this_scoretype = None, dump_best_pose = False, dump_pose_name = None, dump_dir = None, verbose = False, pmm = None ):
    """
    Packs and minimizes a <pose> <inner_trials> times, then uses the best Pose based on total score to pack and minimize again <outer_trials> times using the ScoreFunction <sf>
    :param sf: ScoreFunction
    :param pose: Pose
    :param outer_trials: int( number of times to run <inner_trials> ). Default = 3
    :param inner_trials: int( number of times to pack and minimize before calling that the temporary "best" structure ). Default = 3
    :param compare_using_this_scoretype: str( what specific ScoreType do you want to use for comparison? ). Default = None = total_energy
    :param dump_best_pose: bool( after finding the lowest energy Pose, dump the structure into the current directory (or to <dump_dir>). Default = False
    :param dump_pose_name: bool( filename of the Pose to be dumped ). Default = None (ie. "Best_" + current Pose name)
    :param dump_dir: str( path/to/dir/where/pose/will/be/dumped ). Default = None (ie. current working directory)
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: Pose( the Pose with the lowest total score after the trials of packing and minimization )
    """
    from rosetta import Pose


    # inform user of options depending on how much they want to hear
    if verbose:
        print "\nUsing", outer_trials, "outer loops and", inner_trials, "inner loops to find a low-energy structure based on total score for pose", pose.pdb_info().name().split( '/' )[-1]
    else:
        print "Getting a base structure for your pose using packing and minimization..."
        
    # apply the PyMOL_Mover if passed
    if pmm is not None:
        # used in the function elsewhere to ensure the program doesn't use a broken mover
        pmm_worked = False
        try:
            if verbose:
                print "Watching your pose now"
            pmm.apply( pose )
            pmm_worked = True
        except:
            print "Something was wrong with your PyMOL_Mover -- continuing without watching"
            pass

    # make a dummy best pose - will be replaced as a new best pose is found
    best_pose = Pose()
    best_pose.assign( pose )

    # pack and minimize a temp pose
    for ii in range( outer_trials ):
        pose.assign( best_pose )
        
        # print current outer trial
        if verbose:
            print "Outer trial", ii + 1
            
        # watch pose in PyMOL
        if pmm is not None and pmm_worked:
            pmm.apply( pose )

        # now cycle through the inner pack_min loops
        for jj in range( inner_trials ):
            # print current inner trial
            if verbose:
                print "  Inner trial", jj + 1
                
            # this makes a temporary Pose which eventually gets replaced with an updated best Pose, make the minimization trajectory a little more directed
            temp_pose = Pose()
            temp_pose.assign( pose )
            if pmm is not None and pmm_worked:
                pmm.apply( temp_pose )
                
            # pack and minimize
            temp_pose.assign( do_pack_min( sf, temp_pose ) )
            if pmm is not None and pmm_worked:
                pmm.apply( temp_pose )
            best_pose.assign( compare_these_poses_by_score( sf, best_pose, temp_pose, compare_using_this_scoretype, verbose = verbose ) )


    # dump the Pose, if desired
    if dump_best_pose:
        # set up the correct Pose filename if none given
        if dump_pose_name is None:
            # default Pose name
            dump_pose_name = "Best_" + pose.pdb_info().name()

        # set up the correct dumping directory
        if dump_dir is None:
            # dump directory is the current working directory since it was not specified at input
            dump_dir = os.getcwd()

            if verbose:
                print "Dumping the Pose now named", dump_pose_name, "into", dump_dir

        else:
            # check if given directory is a valid path, else, dump to current working directory instead and inform user
            if not os.path.isdir( dump_dir ):
                print dump_dir, "is not a valid directory. Dumping the Pose now named", dump_pose_name, "to current working directory (", dump_dir, ") instead"
                dump_dir = os.getcwd()

            # the specified dump directory is valid, so the Pose will be dumped there
            else:
                if verbose:
                    print "Dumping the Pose now named", dump_pose_name, "into the directory", dump_dir

        # dump the Pose by navigating to dump directory
        cur_dir = os.getcwd()
        os.chdir( dump_dir )

        # dump the pose and return to the current working directory
        best_pose.dump_pdb( dump_pose_name )
        os.chdir( cur_dir )

    return best_pose



####################
#### LOOPS CODE ####
####################


def make_loop( start, stop, cut = None ):
    """
    Makes a loop starting at <start>, ending at <end>, and with a cutpoint at either <cut>, or the midpoint between <start> and <stop>
    :param start: int( start of loop )
    :param stop: int( end of loop )
    :param cut: int( cutpoint ). Default = None (ie. the midpoint between <start> and <stop>)
    :return: a Loop object
    """
    from rosetta import Loop


    # if a cut was specified on input
    if cut is not None:
        loop = Loop( start, stop, cut )

    # otherwise, make the cut the midpoint between start and stop
    else:
        loop = Loop( start, stop, int( ( start + stop ) / 2 ) )

    return loop



def make_loops( starts, stops, cuts = None ):
    """
    Makes a Loops object out of a list of <starts> and <stops> with an optional list of <cuts>
    :param starts: list( int( start of loops ) )
    :param stops: list( int( end of loops ) )
    :param cuts: list( int( cutpoints ) ). Default = None (ie. the midpoint between <start> and <stop>)
    :return: a Loops object
    """
    from rosetta import Loops


    # check that the arguments passed are lists - need lists for code to work
    if not isinstance( starts, list ):
        print
        print starts, "isn't a list. I need a list of int( start points ). Exiting"
        sys.exit()
    if not isinstance( stops, list ):
        print
        print stops, "isn't a list. I need a list of int( stop points ). Exiting"
        sys.exit()
    if cuts is not None:
        if not isinstance( cuts, list ):
            print
            print cuts, "isn't a list. I need a list of int( cutpoints ). Exiting"
            sys.exit()

    # check to make sure there are as many starts as stops and as cuts (if passed)
    if cuts is None:
        if len( starts ) != len( stops ):
            print
            print "You didn't give me the right amount of starts and ends. Exiting"
            sys.exit()
    else:
        if len( starts ) != len( stops ) != len( cuts ):
            print
            print "You didn't give me the right amount of starts/ends/cuts. Exiting"
            sys.exit()

    # now make the Loops object
    loops = Loops()
    for ii in range( len( starts ) ):
        # if cut points were specified at input, use those
        if cuts is not None:
            loop = make_loop( starts[ii], stops[ii], cuts[ii] )

        # otherwise, make the cutpoints the midpoint between the starts and stops
        else:
            loop = make_loop( starts[ii], stops[ii] )
        loops.add_loop( loop )

    return loops



def make_loops_from_file( loops_file, pose, PDB_numbering = False, anchor_loops = False, verbose = False ):
    """
    Reads the <loops_file> and creates the corresponding <loops>. Returns None if there was an error.
    IMPORTANT: can go from PDB --> Pose numbering ONLY IF in your LOOPs file you include the chain LETTER after the stop or cut value. Set PDB_numbering to True for this feature
    :param loops_file: str( /path/to/LOOPs/file )
    :param pose: Pose
    :param PDB_numbering: bool( is this file numbered using PDB numbering? ) Default = False
    :param anchor_loops: bool( do you want to add +2-residue anchors to the ends of your loops? ). Default = True
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: a Loops object
    """
    from rosetta import Loop, Loops


    # makes the Loops object that will hold all the Loops created
    loops = Loops()

    try:
        # make sure the <loops_file> is a valid file and can be opened
        with open( loops_file, 'rb' ) as fh:
            lines = fh.readlines()
    except IOError:
        # if it's not a valid path, print an error message and exit
        print
        print
        print loops_file, "does not appear to be a valid path. Exiting"
        sys.exit()
    except:
        # if something else happened, print an error message and exit
        print
        print
        print "Something I don't understand happened. Raising the error."
        raise

    if verbose:
        if PDB_numbering:
            print "Making", len( lines ), "Loops from the given LOOPs file and changing from PDB numbering to Pose numbering\n"
        else:
            print "Making", len( lines ), "Loops from the given LOOPs file\n"

    # create the Loops from the file
    for line in lines:
        if line != '':
            # strip out new line characters and white space  -  used to make string parsing more efficient for Loop creation
            line = line.rstrip( '\n' )
            white_space_columns = line.split( ' ' )

            columns = []
            # splits up the line by columns, ignoring anything that is a Blank Space (Taylor Swift approves!)
            for ii in white_space_columns:
                if ii != '':
                    columns.append( ii )

            # the start and stop point should always be at [1] and [2]  -  "EDGE" is always at [0]
            start = int( columns[1] )
            stop = int( columns[2] )
            
            # if anchor_loops is True, add subtract 2 from start and add 2 to stop
            if anchor_loops:
                start -= 2
                stop += 2

            # the first try and the following three excepts are relevant to finding or making the cutpoint
            # in the try case, we are trying to see if an integer cutpoint is specified at columns[3]
            try:
                # an integer cutpoint should be here
                cut = int( columns[3] )

                # check to make sure cut is between the given start and stop, otherwise print an error message and exit
                if not start < cut < stop:
##                if not ( start < cut < stop or stop < cut < start ):
                    print
                    print "The cut you gave me,", cut, ",is not between your start and stop points:", start, "and", stop, " -  Check your file. Exiting"
                    sys.exit()

            # in this case, the cutpoint doesn't exist, but the chain ID does  -  so make cutpoint instead
            # trying to make an integer out of a string results in a ValueError, so this will work for this case, but not for making a string out of an integer
            except ValueError:
                # no cut given, so make one by getting the midpoint of the start and stop points
                cut = int( ( start + stop ) / 2 )

                # ensure a chain ID does indeed exist here by checking if it's a letter
                chain = str( columns[3] )
                if chain not in all_letters_list:
                    print
                    print
                    print "I'm not sure what this is:", "'%s'" %chain, " -  I was expecting a letter from A through Z. Check your file. sys.exit()"
                    sys.exit()

            # in this case, the cutpoint doesn't exist, and there is no chain ID there either  -  so make the cutpoint, or exit if there was supposed to be a chain ID (if <PDB_numbering> is True)
            except IndexError:
                # no cut given, so make one by getting the midpoint of the start and stop points
                cut = int( ( start + stop ) / 2 )

                # chain ID doesn't exist in this case, exit if user set <PDB_numbering> to True. They need to have added a chain here
                if PDB_numbering:
                    print
                    print
                    print "You set PDB_numbering to True, but did not add a chain ID after the stop residue. Exiting"
                    sys.exit()

            # otherwise, something weird happened
            except:
                print
                print
                print "Something I don't understand happened. Raising the error."
                raise

############

            # now, look for a chain ID, if <PDB_numbering> is True (otherwise, it wouldn't be there anyway)
            if PDB_numbering:
                # try to see if anything exists at columns[4]  -  don't need to check columns[3] because that was done above
                try:
                    chain = columns[4]

                # if there is nothing there, exit because <PDB_numbering> is True, thus we need a chain ID
                except IndexError:
                    print
                    print
                    print "You set PDB_numbering to True, but did not add a chain ID after the cut residue. Exiting"
                    sys.exit()

                # if something weird happened, exit
                except:
                    print
                    print
                    print "Something I don't understand happened. Raising the error."
                    raise

                # try except won't work for there being a ValueError because it is allowed to make a string out of an integer, so we need another way to ensure that the chain ID given is a valid letter
                # check to make sure that columns[4] is a string  -  exit if it is an integer or anything else
                if isinstance( chain, str ):
                    # check to make sure that the chain is actually only a letter
                    if chain not in all_letters_list:
                        print
                        print "I'm not sure what this is:", "'%s'" %chain, " -  I was expecting a letter from A through Z. Check your file. Exiting"
                        sys.exit()

                # if it's an integer, exit since we are expecting a string
                elif isinstance( chain, int ):
                    print
                    print "I was not expecting an integer here, I need a string for the chain ID. Exiting"
                    sys.exit()

                # otherwise, I'm not sure what is in columns[4], so exit
                else:
                    print
                    print "I'm not sure what this is", columns[4], "I need a string for the chain ID. Exiting"
                    sys.exit()


            # if <PDB_numbering> is True, go from PDB to Pose numbering using the positions and chain IDs from file
            if PDB_numbering:
                start = pose.pdb_info().pdb2pose( chain, start )
                stop = pose.pdb_info().pdb2pose( chain, stop )
                cut = pose.pdb_info().pdb2pose( chain, cut )

            # add each Loop to the Loops object
            loop = Loop( start, stop, cut )
            loops.add_loop( loop )

    # print loop to user if verbose is True
    if verbose:
        print "Returning the following Loops object\n\t", loops
        
    return loops



def add_cutpoint_variants( loops, pose ):
    """
    Adds cutpoint_variants to the <pose> for each loop given in <loops>
    :param loops: a Loops object (can be one loop)
    :param pose: Pose
    :return: Pose( <pose> with <loops> that are cutpoint_variants )
    """
    from rosetta import add_single_cutpoint_variant


    # for every Loop in the Loops object, add it as a cutpoint variant to the <pose>
    for loop_num in range( 1, loops.num_loop() + 1 ):
        add_single_cutpoint_variant( pose, loops[ loop_num ] )

    return pose



########################
#### FOLD TREE CODE ####
########################

def foldtree_to_string_and_stripped( pose, verbose = False ):
    """
    Turns the FoldTree of <pose> into a list of strings only containing the relevant numbers to each EDGE
    Removes the words FOLD_TREE and EDGE and strips whitespace and splits the rest into chunks
    Example output:
    :param pose: Pose
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: list( relevant information from each EDGE)
    """
    if verbose:
        print "Turning the Pose's FoldTree into a string by splitting on 'EDGE' and stripping off whitespace"

    # get string version of FoldTree
    ft = pose.fold_tree().to_string()
    ft = ft.rstrip()

    # split FoldTree on "EDGE"
    ft_split_on_edge = ft.split( "EDGE" )

    # split on white space to get access to numbers and remove the FOLD_TREE tag at the beginning
    ft_split_on_white_space = []
    for edge in ft_split_on_edge:
        # strips off white space from the ends
        edge = edge.strip()
        # ignores the FOLD_TREE chunk
        if edge != "FOLD_TREE":
            # ignores any blank spaces
            if edge != '':
                ft_split_on_white_space.append( edge )

    return ft_split_on_white_space



def get_JUMP_NUM_from_seq_pos( seq_pos, pose, downstream = False, verbose = False ):
    """
    Return the associated up- or downstream jump number given a specific <seq_pos> of a <pose>
    Returns the upstream jump number by default, set <downstream> to True to retrieve upstream jump number of <seq_pos>
    :param seq_pos: int( sequence position of residue in question )
    :param pose: Pose
    :param downstream: bool( returns upstream residue by default, set <downstream> to True to return downstream jump num )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: int( up- or downstream jump number )
    """
    # first check to see that the <seq_pos> entered is a valid position
    if seq_pos <= 0 or seq_pos > pose.n_residue():
        print
        print "Invalid sequence position given:", seq_pos, "Exiting"
        sys.exit()

    if verbose:
        if downstream:
            print "Given the residue sequence number", seq_pos, "this will find the closest downstream Jump number"
        else:
            print "Given the residue sequence number", seq_pos, "this will find the closest upstream Jump number"

    # retrieve both upstream and downstream jump numbers relative to <seq_pos>
    upstream_jump = 0
    downstream_jump = 0
    last_jump_seen = 0

    # get string version without whitespace of FoldTree EDGE's
    ft_split_on_white_space = foldtree_to_string_and_stripped( pose )

    # split up and use each EDGE from the FoldTree
    for edge in ft_split_on_white_space:
        # split on whitespace to get only numbers
        edge_split = edge.split( ' ' )
        edge_split_on_white = []

        for e in edge_split:
            # ignore resulting white space and blank spaces
            if e != ' ':
                if e != '':
                    # append only the relevant numbers for this EDGE
                    edge_split_on_white.append( e )

        # get relevant EDGE numbers  -  not using start because the end residue is more relevant to use when looking for the closest Jump number
        end_res = int( edge_split_on_white[1] )
        edge_type = int( edge_split_on_white[2] )

        # if this is a jump EDGE  (if edge_type is a positive integer)
        if edge_type > 0:
            # store it as the last_jump_seen to be used in case the residue in question is beyond the last Jump
            last_jump_seen = edge_type

            # if the <seq_pos> is within this Jump EDGE  (or at least less than the end of it)
            if seq_pos <= end_res:
                # store both the upstream and downstream jump
                upstream_jump = int( edge_type )
                downstream_jump = upstream_jump - 1

                # return upstream or downstream jump depending on user input when it is found  -  no need to go through all of the FoldTree
                if downstream:
                    return downstream_jump
                else:
                    return upstream_jump

    # if the last_jump_seen != 0 but upstream_jump == 0, it means there are jumps, but the seq_pos is beyond the last j=Jump. So set the Jump to last_jump_seen
    if last_jump_seen > 0 and upstream_jump == 0:
        upstream_jump = last_jump_seen
        downstream_jump = last_jump_seen - 1

    # return upstream or downstream jump depending on user input
    if downstream:
        return downstream_jump
    else:
        return upstream_jump



def get_seq_pos_from_JUMP_NUM( JUMP_NUM, pose, start = False, verbose = False):
    """
    Return the associated start or stop sequence position corresponding to the <JUMP_NUM> in a <pose>
    By default, returns the stop position. Set <start> to True if you want the start position
    :param JUMP_NUM: int( Jump number )
    :param pose: Pose
    :param start: bool( do you want the start sequence position? Please set to True). Default = False (ie. returns stop position )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: int( sequence position of the start or stop residue )
    """
    # check to make sure a valid <JUMP_NUM> was passed
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print "Invalid jump number given:", JUMP_NUM, "Exiting"
        sys.exit()

    if verbose:
        if start:
            print "Getting the start residue closest to the Jump number", JUMP_NUM
        else:
            print "Getting the stop residue closest to the Jump number", JUMP_NUM

    # get string version without whitespace of FoldTree EDGE's
    ft_split_on_white_space = foldtree_to_string_and_stripped( pose )

    # split up and use each EDGE from the FoldTree
    for edge in ft_split_on_white_space:
        # split on whitespace to get only numbers
        edge_split = edge.split( ' ' )
        edge_split_on_white = []

        for e in edge_split:
            # ignore resulting white space and blank spaces
            if e != ' ':
                if e != '':
                    # append only the relevant numbers for this EDGE
                    edge_split_on_white.append( e )


        # find the start or stop residue of the given <JUMP_NUM>
        start_pos = int( edge_split_on_white[0] )
        stop_pos = int( edge_split_on_white[1] )
        edge_type = int( edge_split_on_white[2] )

        # if the edge_type is the <JUMP_NUM>
        if edge_type == JUMP_NUM:
            # return start or stop position, given user input
            if start:
                return start_pos
            else:
                return stop_pos



# LOOPs file --> appropriate fold tree
def hard_code_fold_tree( orig_loops, pose, verbose = False ):
    from rosetta import FoldTree, Loops


    # instantiate an empty FoldTree
    ft = FoldTree()
    
    # make a copy of the original loops as to not alter them for later use
    loops = Loops()
    for ii in range( 1, orig_loops.num_loop() + 1 ):
        loops.add_loop( orig_loops[ii].start(), orig_loops[ii].stop(), orig_loops[ii].cut() )
    if verbose:
        print "Hard coding FoldTree with the following Loops object:\n\t", loops
    
    # attach the Loops object to the <pose> to be used when tearing down the FoldTree with restore_original_fold_tree
    pose.loops = Loops()
    for ii in range( 1, loops.num_loop() + 1 ):
        pose.loops.add_loop( loops[ii] )
    
    # print original FoldTree, if verbose
    if verbose:
        print "Original FoldTree:\n\t", pose.fold_tree()
    
    # parse through the loops first
    starts_from_loops = []
    for ii in range( 1, loops.num_loop() + 1 ):
        starts_from_loops.append( int( loops[ii].start() ) )  # used for indexing
        
    hold_ft = []
    break_out = False
    JUMP_NUM = 1
    
    orig_ft = pose.fold_tree().to_string().split( "EDGE" )
    hold_ft = []
    for edge in orig_ft:
        edge = edge.strip()
        if edge != "FOLD_TREE":
            if edge != '':
                hold_ft.append( edge )
                        
    # outer loop should be the loop_starts in case there are two loops within one EDGE
    for edge in hold_ft:
        # split on whitespace
        temp_columns = edge.split( ' ' )
        columns = []
        for ii in temp_columns:
            if ii != '':
                columns.append( ii )
                    
        # check to see that this is a -1 EDGE ( this is the only EDGE I'd be altering anyway )
        if int( columns[2] ) == -1:
            loops_to_add = []
            for loop_start in starts_from_loops:
                # second number will be the end of the edge or jump, check to see that the new loop start isn't smaller than this sequence position
                if loop_start < int( columns[1] ):
                    loops_to_add.append( loop_start )
                        
            if len( loops_to_add ) != 0:
                if len( loops_to_add ) == 1:
                    for loop_start in loops_to_add:
                        # split this EDGE up to make room for the loop
                        # get the loop's index number within the loops object to get it's end and cut info
                        loop_index_number = starts_from_loops.index( loop_start ) + 1  # Loops are indexed starting at 1
                        loop = loops[ loop_index_number ]
                        
                        # add the edge from the original start to this loop's start point
                        ft.add_edge( int( columns[0] ), loop.start(), -1 )
                        
                        # add the edge from the start of the loop to its cutpoint
                        ft.add_edge( loop.start(), loop.cut(), -1 )
                        
                        # add the jump from the start to the end of the loop
                        ft.add_edge( loop.start(), loop.stop(), JUMP_NUM )
                        JUMP_NUM += 1
                        
                        # add the edge from the end of the loop to its upper cutpoint
                        ft.add_edge( loop.stop(), loop.cut() + 1, -1 )
                        
                        # finally, go from end of loop to the original end point
                        ft.add_edge( loop.stop(), int( columns[1] ), -1 )
                        
                        # remove the loop from the loops list
                        starts_from_loops.remove( loop_start )
                        loops.delete_loop( loop.start(), loop.stop() )
                else:
                    ii = 0
                    jj = 1
                    while ( ii + 1 ) != len( loops_to_add ):
                        loop_index_number = starts_from_loops.index( loops_to_add[ii] ) + 1  # Loops are indexed starting at 1
                        loop1 = loops[ loop_index_number ]
                        
                        loop_index_number = starts_from_loops.index( loops_to_add[jj] ) + 1  # Loops are indexed starting at 1
                        loop2 = loops[ loop_index_number ]
                        
                        if ii == 0:
                            # the start is from the current edge start to the new loop start
                            ft.add_edge( int( columns[0] ), loop1.start(), -1 )
                            
                            # add the edge from the start of the loop to its cutpoint
                            ft.add_edge( loop1.start(), loop1.cut(), -1 )
                            
                            # add the jump from the start to the end of the loop
                            ft.add_edge( loop1.start(), loop1.stop(), JUMP_NUM )
                            JUMP_NUM += 1
                            
                            # add the edge from the end of the loop to its upper cutpoint
                            ft.add_edge( loop1.stop(), loop1.cut() + 1, -1 )
                            
                            # finally, go from end of loop to the start of the next loop
                            ft.add_edge( loop1.stop(), loop2.start(), -1 )
                            
                            # remove the loop from the loops list
                            if ii != ( len( loops_to_add ) - 1 ):
                                starts_from_loops.remove( loop1.start() )
                                loops.delete_loop( loop1.start(), loop1.stop() )
                            
                            # up the loop indices
                            ii += 1
                            jj += 1
                                
                        else:
                            # add the edge from the start of the loop to its cutpoint
                            ft.add_edge( loop1.start(), loop1.cut(), -1 )
                            
                            # add the jump from the start to the end of the loop
                            ft.add_edge( loop1.start(), loop1.stop(), JUMP_NUM )
                            JUMP_NUM += 1
                            
                            # add the edge from the end of the loop to its upper cutpoint
                            ft.add_edge( loop1.stop(), loop1.cut() + 1, -1 )
                            
                            # finally, go from end of loop to the start of the next loop
                            ft.add_edge( loop1.stop(), loop2.start(), -1 )
                            
                            # remove the loop from the loops list
                            if ii != ( len( loops_to_add ) - 2 ):
                                starts_from_loops.remove( loop2.start() )
                                loops.delete_loop( loop2.start(), loop2.stop() )
                                
                            # up the loop indices
                            ii += 1
                            jj += 1
                            
                    # add the final loop
                    loop_start = loops_to_add[-1]
                    loop_index_number = starts_from_loops.index( loop_start ) + 1  # Loops are indexed starting at 1
                    loop2 = loops[ loop_index_number ]
                    
                    # add the edge from the start of the loop to its cutpoint
                    ft.add_edge( loop2.start(), loop2.cut(), -1 )
                    
                    # add the jump from the start to the end of the loop
                    ft.add_edge( loop2.start(), loop2.stop(), JUMP_NUM )
                    JUMP_NUM += 1
                    
                    # add the edge from the end of the loop to its upper cutpoint
                    ft.add_edge( loop2.stop(), loop2.cut() + 1, -1 )
                    
                    # finally, go from end of loop to the original end point
                    ft.add_edge( loop2.stop(), int( columns[1] ), -1 )
                    
                    # remove the loop from the loops list
                    starts_from_loops.remove( loop2.start() )
                    loops.delete_loop( loop2.start(), loop2.stop() )
                    
            else:
                ft.add_edge( int( columns[0] ), int( columns[1] ), -1 )
                    
        if int( columns[2] ) == -2:
            # chemical edge ( sugar )
            ft.add_edge( int( columns[0] ), int( columns[1] ), columns[3], columns[4] )
            
        if int( columns[2] ) > 0:
            # a jump
            ft.add_edge( int( columns[0] ), int( columns[1] ), JUMP_NUM )
            JUMP_NUM += 1
            
    if ft.check_fold_tree():
        if verbose:
            print "New FoldTree:\n\t", ft
        return ft
    else:
        print "There was a problem hard coding your fold tree"
        sys.exit()



def setup_new_fold_tree( loops_file, pose, PDB_numbering = False, anchor_loops = True, add_cutpoints = True, verbose = False ):
    """
    Creates a FoldTree given a <loops_file> and the <pose>'s current FoldTree and then gives it to <pose>
    If the Loops result in a invalid FoldTree, the program will result in an error and exit
    :param loops_file: str( /path/to/valid/LOOPs/file )
    :param pose: Pose
    :param PDB_numbering: bool( is this file numbered using PDB numbering? ) Default = False
    :param anchor_loops: bool( do you want to add +2-residue anchors to the ends of your loops? ). Default = True
    :param add_cutpoints bool( do you want to add cutpoints in the FoldTree for each of the loop's cutpoints? ). Default = True
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: Pose with the new FoldTree
    """
    from rosetta import Loops


    # attach the Loops object to the <pose> to be used when tearing down the FoldTree with restore_original_fold_tree
    pose.loops = Loops()
    
    # check to see if user passed a Loop(s) object or the file name of a Loop(s) object
    if isinstance( loops_file, str ):
        # add the name of the loops file to the pose to be used when returning it to its normal FoldTree
        pose.loops_file = loops_file
        
        # make the Loops from LOOPs file
        loops = make_loops_from_file( loops_file, pose, PDB_numbering = PDB_numbering, anchor_loops = anchor_loops, verbose = verbose )
        for ii in range( 1, loops.num_loop() + 1 ):
            pose.loops.add_loop( loops[ii] )

    else:
        print
        print "I'm not quite sure what you passed to me. I need the name of the file containing your Loops, or an actual Loop or Loops object. Exiting."
        sys.exit()
    
    # apply the new FoldTree by altering the current one to include the new Loops
    pose.fold_tree( hard_code_fold_tree( loops, pose, verbose = verbose ) )

    # add cutpoint variants for each new Loop
    if add_cutpoints:
        if verbose:
            print "Pose's original set of cutpoints:\n\t", pose.fold_tree().cutpoints()
        pose.assign( add_cutpoint_variants( loops, pose ) )
        if verbose:
            print "Pose's new set of cutpoints:\n\t", pose.fold_tree().cutpoints()
    
    return pose



def restore_original_fold_tree( pose, verbose = False ):
    """
    Returns the Pose's FoldTree back to its original version and removes the added cutpoint variants
    The Pose's original FoldTree was attached to itself during the function load_pose
    :param pose: Pose
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: Pose with original FoldTree
    """
    from rosetta.core.chemical import VariantType
    from rosetta.core.pose import remove_variant_type_from_pose_residue


    # check to make sure there's actually a reason to restore the original FoldTree
    # if pose.loops_file is still None, that means it hasn't been changed  -  therefore nothing to do for this <pose>
    if pose.loops_file is None:
        print
        print "I see you haven't actually changed the fold tree. Is this true?"
        return pose

    # otherwise, restore it
    else:
        if verbose:
            print "Returning the Pose's FoldTree back to its original one by removing the cutpoint variants added and then applying the original FoldTree stored from from the load_pose function"

        ## IMPORTANT: have to get the loops again because the fold tree thinks the sugars are cutpoints - will mess it up
        loops = pose.loops
        if verbose:
            print "  Tearing down the following Loops:\n\t", loops

        # remove cutpoint variants by reading the LOOPs file attached to the <pose>
        for ii in range( 1, loops.size() + 1 ):
            # make the Loop from the Loops object
            loop = loops[ii]

            # get its lower and upper cutpoints
            lower_cut = loop.cut()
            upper_cut = lower_cut + 1

            # remove the variant types for the Loop from the LOOPs file
            remove_variant_type_from_pose_residue( pose, VariantType.CUTPOINT_LOWER, lower_cut )
            remove_variant_type_from_pose_residue( pose, VariantType.CUTPOINT_UPPER, upper_cut )

        # apply the original FoldTree back to the <pose>
        pose.fold_tree( pose.orig_fold_tree )
        
        if verbose:
            print "Restored Pose to following FoldTree:\n\t", pose.fold_tree()
            print "Pose now has the following cutpoints:\n\t", pose.fold_tree().cutpoints()

        # return the LOOPs file object in the <pose> back to None since it's now back to original FoldTree
        pose.loops_file = None
        pose.loops = None

        return pose





#####################################
#### MUTATIONAL WORKER FUNCTIONS ####
#####################################

def mutate_residue( pose_num, new_res_name, input_pose, sf, pdb_num = False, pdb_chain = None, pack_radius = 5, do_pack = True, do_min = True, do_full_pack_min = False ):
    """
    Mutate residue at position <pose_num> to <new_res_name>
    <new_res_name> can be a single-letter or three-letter residue code
    If you are giving a pdb number, set <pdb_num> to True AND give me a <pdb_chain> letter id
    :param pose_num: int( Pose number for residue )
    :param new_res_name: str( one- or three-letter code for the new amino acid. Example 'A' or "THR" )
    :param input_pose: Pose
    :param sf: ScoreFunction ( used for packing )
    :param pdb_num: bool( did you give me a PDB number instead? Set to True if so. Give me a <pdb_chain> too then ) Default = False (Pose number)
    :param pdb_chain: str( PDB chain id such as 'A' or 'X'. Must have set <pdb_num> to True as well
    :param pack_radius: int or float( how far out in Angstroms do you want to pack around the mutation site? ) Default = 5
    :param do_pack: bool( do you want to pack around the mutation? ) Default = True
    :param do_min: bool( do you want to minimize around the mutation? ) Default = True
    :param do_full_pack_min: bool( instead of doing a local pack and/or min, do you want to do it to the whole Pose instead? ) Default = False
    :return: mutated Pose
    """
    # imports
    from rosetta import Pose, pose_from_sequence, ResidueFactory, MoveMap, MinMover


    # copy over the input pose
    pose = input_pose.clone()

    # check if <pdb_chain> was given if <pdb_num> is True
    if pdb_num == True:
        if pdb_chain is None:
            print "\nYou told me you gave me a PDB number, but you did not give me a PDB chain id. Set <pdb_chain> to the appropriate chain id. Returning the original pose."
            return pose

    # check the logic of the input arguments
    # ensure the <new_res_name> is a valid residue
    if len( new_res_name ) != 1 and len( new_res_name ) != 3:
        print "\nYou did not give me a single- or three-letter amino acid code. '%s' did not work. Returning the original pose." %new_res_name
        return pose
    # if it is a single-letter code
    if len( new_res_name ) == 1 and new_res_name.upper() not in AA_name1_list:
        print "\nIt appears that '%s' is not a valid single-letter amino acid code. Returning the original pose." %new_res_name
        return pose
    # if it is a three-letter code
    elif len( new_res_name ) == 3 and new_res_name.upper() not in AA_name3_list:
        print "\nIt appears that '%s' is not a valid three-letter amino acid code. Returning the original pose." %new_res_name
        return pose
   # otherwise, use the <new_res_name> argument to get the appropriate three-letter amino acid code
    if len( new_res_name ) == 1:
        single_new_res_name = new_res_name.upper()
        new_res_name = AA_name1_to_name3[ single_new_res_name ]
    else:
        new_res_name = new_res_name.upper()
        single_new_res_name = AA_name3_to_name1[ new_res_name ]

    # ensure <pose_num> (and <pdb_chain>) exists in the pose
    if not pdb_num:
        if not 1 <= pose_num <= pose.n_residue():
            print "\nYou appear to have given me an invalid Pose residue number. Ensure residue number %s exists in your Pose. Returning the original pose." %pose_num
            return pose
    # if it's a PDB number, check it exists as well using the <pdb_chain> too
    else:
        # get the actual pose number
        pose_num = pose.pdb_info().pdb2pose( pdb_chain, pose_num )
        if pose_num == 0:
            print "\nYour PDB number and chain ( %s chain %s ) don't seem to exist in the pose. Check your input. Returning the original pose." %( pose_num, pdb_chain )
            return pose

    # move on to the mutation
    # instantiate a ResidueFactory
    res_factory = ResidueFactory()

    # create a three-mer of the <new_res_name> desired
    # want a three-mer because it's easier to deal with a new amino acid that does not have a special end VariantType
    threemer = pose_from_sequence( single_new_res_name * 3 )

    # get the ResidueType from the middle <new_res_name> in the threemer
    res_type = threemer.conformation().residue_type( 2 )

    # build the new residue and preserve the CB information from the original pose
    new_residue = res_factory.create_residue( res_type, 
                                              current_rsd = pose.residue( pose_num ), 
                                              conformation = pose.conformation(), 
                                              preserve_c_beta = True )

    # replace the old residue in the pose
    pose.replace_residue( pose_num, new_residue, orient_backbone = True )

    # if a pack and/or min is happening
    if do_pack or do_min:
        # if they don't want to do a full Pose pack and/or min
        if not do_full_pack_min:
            # get residue numbers (including mutation site) to be packed and minimized
            res_nums_around_mutation_site = get_res_nums_within_radius( pose_num, pose, pack_radius, include_res_num = True )
        # if a full pack and/or min is desired, make the residue range the full size of the pose
        else:
            res_nums_around_mutation_site = range( 1, pose.n_residue() + 1 )

    # pack around mutation, if desired
    if do_pack:
        pack_rotamers_mover = make_pack_rotamers_mover( sf, pose,
                                                        pack_branch_points = True,
                                                        residue_range = res_nums_around_mutation_site )
        pack_rotamers_mover.apply( pose )

    # minimize around mutation, if desired
    if do_min:
        min_mm = MoveMap()
        # turn on packing for the bb and chi
        for res_num in res_nums_around_mutation_site:
            min_mm.set_bb( res_num, True )
            min_mm.set_chi( res_num, True )
        # create and apply the MinMover
        min_mover = MinMover( movemap_in = min_mm,
                              scorefxn_in = sf,
                              min_type_in = "dfpmin_strong_wolfe",
                              tolerance_in = 0.01,
                              use_nb_list_in = True )
        min_mover.apply( pose )

    return pose



def make_mutation_packer_task( amino_acid, seq_pos, sf, pose, pack_radius = PACK_RADIUS ):
    """
    OUTDATED AS OF 17 AUGUST 2016 BUT KEEPING FOR REFERENCE
    Returns a packer task that can handle a <pose> with a SINGLE mutation at <seq_pos>
    :param amino_acid: str( ONE LETTER amino acid code for the SINGLE mutated residue )
    :param seq_pos: int( sequence position of the mutated residue )
    :param sf: ScoreFunction
    :param pose: Pose
    :param pack_radius: int( or float( the distance in Angstroms around the mutated residue you want to be packed ). Default = PACK_RADIUS = 20
    :return: packer_task ready to pack a Pose with a SINGLE mutation
    """
    from rosetta import standard_packer_task, RotamerTrialsMover, \
        aa_from_oneletter_code
    from rosetta.utility import vector1_bool


    # create a packer task handling a SINGLE mutated residue
    # tell packer task that only the current <amino_acid> should be allowed to be modified
    aa_bool = vector1_bool()
    mutant_amino_acid = aa_from_oneletter_code( amino_acid )
    for ii in range( 1, 20 + 1 ):
        aa_bool.append( ii == mutant_amino_acid )

    # get center of mutated residue to be used in distance calculation for pack radius
    center = pose.residue( seq_pos ).nbr_atom_xyz()

    # make the packer task
    task = standard_packer_task( pose )
    task.or_include_current( True )
    task.nonconst_residue_task( seq_pos ).restrict_absent_canonical_aas( aa_bool )
    task.restrict_to_repacking()

    # turn off repacking for residues not in the pack radius
    for ii in range( 1, pose.total_residue() + 1 ):
        if center.distance( pose.residue( ii ).nbr_atom_xyz() ) > pack_radius:
            task.nonconst_residue_task( ii ).prevent_repacking()
    pack_rotamers_mover = RotamerTrialsMover( sf, task )

    return pack_rotamers_mover



def do_mutation_pack( seq_pos, amino_acid, sf, mutated_pose, pack_radius = PACK_RADIUS ):
    """
    OUTDATED AS OF 17 AUGUST 2016 BUT KEEPING FOR REFERENCE
    Returns a packed and minimized <pose> that has a SINGLE mutation of <amino_acid> at <seq_pos>
    :param seq_pos: int( sequence position of the mutated residue _
    :param amino_acid: str( ONE LETTER amino acid code for the SINGLE mutated residue )
    :param sf: ScoreFunction
    :param mutated_pose: Pose with the SINGLE mutation
    :param pack_radius: int( or float( the distance in Angstroms around the mutated residue you want to be packed ). Default = PACK_RADIUS = 20
    :return: Pose packed around the SINGLE mutation
    """
    from rosetta import Pose


    pose = Pose()
    pose.assign( mutated_pose )
    
    # pack the <mutated_pose>
    pack_rotamers_mover = make_mutation_packer_task( amino_acid, seq_pos, sf, pose, pack_radius )
    pack_rotamers_mover.apply( pose )

    return pose



###################################
#### MAIN MUTANT POSE CREATION ####
###################################

def get_best_3ay4_mutant_of_20_symmetrical( pose_num, sf, input_pose, pdb_num = False, pack_radius = 5, verbose = False, pmm = None ):
    """
    For a single position <pos_num> (both chain A and chain B), mutate to all 20 amino acids and return the best mutant pose
    You can give a Pose residue position or a PDB residue position. Just specify if it is a PDB position by setting <pdb_num> to True
    :param pose_num: int( Pose number for residue )
    :param sf: ScoreFunction ( used for packing )
    :param input_pose: Pose
    :param pdb_num: bool( did you give me a PDB number instead? Set to True if so ) Default = False
    :param pack_radius: int or float( how far out in Angstroms do you want to pack around the mutation site? ) Default = 5
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: the mutated Pose of the lowest total score out of the twenty mutations
    """
    # make a copy of the pose
    pose = input_pose.clone()

    # if desired, send the native Pose to PyMOL
    if pmm is not None:
        pose.pdb_info().name( "native" )
        pmm.apply( pose )

    # get the original amino acid at this position
    # if the PDB number was given, get the Pose number
    if pdb_num is True:
        pose_num1 = pose.pdb_info().pdb2pose( 'A', pose_num )
        pose_num2 = pose.pdb_info().pdb2pose( 'B', pose_num )

        # make sure the amino acid at the PDB position specified is the same on chain A as it is on chain B
        if pose.residue( pose_num1 ).name1() != pose.residue( pose_num2 ).name1():
            print "There is something wrong with your symmetrical Pose. The residues at position", pose_num, "on chain A and B are not the same. Exiting."
            sys.exit()

        # get the original amino acid
        orig_amino_acid = pose.residue( pose_num1 ).name1()
    # else, if the Pose number was given, get the original amino acid directly
    else:
        orig_amino_acid = pose.residue( pose_num ).name1()

    # get the starting score
    start_score = sf( pose )

    # talk to user, if desired
    if verbose:
        if pdb_num is True:
            print "Mutating", orig_amino_acid, "at PDB position", pose_num, "on chain A and B"
        else:
            print "Mutating", orig_amino_acid, "at Pose position", pose_num

    # create data holders
    # mutants: key = new amino acid, value = mutant Pose
    mutants = {}

    # mutant_scores will be used to keep track of the lowest E mut
    # key = new amino acid, value = mutant Pose score
    mutant_scores = {}
    
    # mutate this position to all 20 amino acids
    for amino_acid in AA_name1_list:
        # make a copy of the pose
        mutant_pose = pose.clone()

        if verbose:
            if orig_amino_acid == amino_acid:
                print "Now mutating", orig_amino_acid, "back to", amino_acid, "..."
            else:
                print "Now mutating", orig_amino_acid, "to", amino_acid, "..."

        # mutate to a new residue
        mutant_pose = mutate_residue( pose_num, amino_acid, mutant_pose, sf, pdb_num = pdb_num, pdb_chain = 'A', pack_radius = pack_radius )
        mutant_pose = mutate_residue( pose_num, amino_acid, mutant_pose, sf, pdb_num = pdb_num, pdb_chain = 'B', pack_radius = pack_radius )

        # if desired, send the mutant Pose to PyMOL
        if pmm is not None:
            mutant_pose.pdb_info().name( orig_amino_acid + "to" + amino_acid )
            pmm.apply( mutant_pose )

        # add the mutant_pose to the dict with the key being the mutation amino acid
        mutants[ amino_acid ] = mutant_pose

        # collect the energy information of this mutant
        mutant_scores[ amino_acid ] = sf( mutant_pose )

    # find the best mutant Pose by using mutant_scores dict information
    lowest_score = None
    best_mutation = None
    for amino_acid in mutant_scores.keys():
        # if this is the first mutation checked, update lowest_score and amino_acid with this mutation's information
        if lowest_score is None:
            lowest_score = mutant_scores[ amino_acid ]
            best_mutation = amino_acid
        # otherwise it is not empty, so compare the other mutations to this one. Keep the lowest score
        else:
            if mutant_scores[ amino_acid ] < lowest_score:
                lowest_score = mutant_scores[ amino_acid ]
                best_mutation = amino_acid

    if verbose:
        ddG = lowest_score - start_score
        print "The best mutation was from", orig_amino_acid, "to", best_mutation, "with a total_score ddG of", ddG

    # return the best mutant from the mutants dictionary using the score information
    return mutants[ best_mutation ]



def get_best_3ay4_mutant_of_20_asymmetrical( pose_num, sf, input_pose, pdb_num = False, pdb_chain = None, pack_radius = 5, verbose = False, pmm = None ):
    """
    For a single position <pos_num>, mutate to all 20 amino acids and return the best mutant pose
    You can give a Pose residue position or a PDB residue position. Just specify if it is a PDB position and give the pdb_chain
    :param pose_num: int( Pose number for residue )
    :param sf: ScoreFunction ( used for packing )
    :param input_pose: Pose
    :param pdb_num: bool( did you give me a PDB number instead? Set to True if so. Give me a <pdb_chain> too then ) Default = False (Pose number)
    :param pdb_chain: str( PDB chain id such as 'A' or 'X'. Must have set <pdb_num> to True as well
    :param pack_radius: int or float( how far out in Angstroms do you want to pack around the mutation site? ) Default = 5
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: the mutated Pose of the lowest total score out of the twenty mutations
    """
    # make a copy of the pose
    pose = input_pose.clone()

    # argument check
    # check if <pdb_chain> was given if <pdb_num> is True
    if pdb_num == True:
        if pdb_chain is None:
            print "\nYou told me you gave me a PDB number, but you did not give me a PDB chain id to create an asymmetrical antibody. Set <pdb_chain> to the appropriate chain id. Returning the original pose."
            return pose

    # if desired, send the native Pose to PyMOL
    if pmm is not None:
        pose.pdb_info().name( "native" )
        pmm.apply( pose )

    # get the original amino acid at this position
    # if the PDB number was given, get the Pose number
    if pdb_num is True:
        pose_num1 = pose.pdb_info().pdb2pose( pdb_chain, pose_num )

        # get the original amino acid
        orig_amino_acid = pose.residue( pose_num1 ).name1()
    # else, if the Pose number was given, get the original amino acid directly
    else:
        orig_amino_acid = pose.residue( pose_num ).name1()

    # get the starting score
    start_score = sf( pose )

    # talk to user, if desired
    if verbose:
        if pdb_num is True:
            print "Mutating", orig_amino_acid, "at PDB position", pose_num, "on chain", pdb_chain
        else:
            print "Mutating", orig_amino_acid, "at Pose position", pose_num

    # create data holders
    # mutants: key = new amino acid, value = mutant Pose
    mutants = {}

    # mutant_scores will be used to keep track of the lowest E mut
    # key = new amino acid, value = mutant Pose score
    mutant_scores = {}
    
    # mutate this position to all 20 amino acids
    for amino_acid in AA_name1_list:
        # make a copy of the pose
        mutant_pose = pose.clone()

        if verbose:
            if orig_amino_acid == amino_acid:
                print "Now mutating", orig_amino_acid, "back to", amino_acid, "..."
            else:
                print "Now mutating", orig_amino_acid, "to", amino_acid, "..."

        # mutate to a new residue
        mutant_pose = mutate_residue( pose_num, amino_acid, mutant_pose, sf, pdb_num = pdb_num, pdb_chain = pdb_chain, pack_radius = pack_radius )

        # if desired, send the mutant Pose to PyMOL
        if pmm is not None:
            mutant_pose.pdb_info().name( orig_amino_acid + "to" + amino_acid + "_chain" + pdb_chain )
            pmm.apply( mutant_pose )

        # add the mutant_pose to the dict with the key being the mutation amino acid
        mutants[ amino_acid ] = mutant_pose

        # collect the energy information of this mutant
        mutant_scores[ amino_acid ] = sf( mutant_pose )

    # find the best mutant Pose by using mutant_scores dict information
    lowest_score = None
    best_mutation = None
    for amino_acid in mutant_scores.keys():
        # if this is the first mutation checked, update lowest_score and amino_acid with this mutation's information
        if lowest_score is None:
            lowest_score = mutant_scores[ amino_acid ]
            best_mutation = amino_acid
        # otherwise it is not empty, so compare the other mutations to this one. Keep the lowest score
        else:
            if mutant_scores[ amino_acid ] < lowest_score:
                lowest_score = mutant_scores[ amino_acid ]
                best_mutation = amino_acid

    if verbose:
        ddG = lowest_score - start_score
        print "The best mutation was from", orig_amino_acid, "to", best_mutation, "with a total_score ddG of", ddG

    # return the best mutant from the mutants dictionary using the score information
    return mutants[ best_mutation ]



def read_mutation_file( mutation_filepath ):
    """
    Return a list of mutations to be made by reading a mutation file designated by the <mutation_filepath>
    Assumes PDB numbering!!!
    Chain designations are separated by '_', multiple mutations are separated by '-'. See below for examples
    Second column can be a ratio of binding constants, if known (i.e. the Shields paper). Some kind of ordered ranking system, basically
    Five possible formats for a mutation string
    1) Single point mutation on both sides (symmetrical)
       A123T (Ala at position 123 on chain A and B to Thr)
    2) Single point mutation on a specific side (asymmetrical)
       A123T_B (Ala at position 123 on chain B to Thr)
    3) Multiple point mutations on both sides (symmetrical)
       A123T-S456Y (Ala at position 123 to Thr and Ser at position 456 to Tyr both on chain A and B)
    4) Multiple point mutations on a specific chain (asymmetrical)
       A123T_A-S456Y_B (Ala at position 123 to Thr on chain A and Ser at position 456 to Tyr on chain B)
    5) Multiple point mutations with a combination of symmetrical and asymmetrical
       A123T-S456Y_B (Ala at position 123 to Thr on chain A and B and Ser at position 456 to Try on chain B)
    Chain designation is optional. If the mutation should be made on both sides, then don't add a chain designation
    This is designed to read mutations for 3ay4 at the moment as it assumes chain A and B
    Commented lines in the file are ignored are ignored
    :param mutation_filepath: str( /path/to/file with mutation strings desired )
    :return: list( mutations to be made )
    """
    # try to open the file
    try:
        f = open( mutation_filepath, "rb" )
        lines = f.readlines()
    except:
        print "Something is wrong with your <mutation_filepath> ( %s ). Please check your input." %mutation_filepath
        raise


    # using an object as to hold mutation names and ratios, if the ratios are in the file
    data_holder = DataHolder()

    # for each line specifying a mutation
    all_mutations = []
    mutation_to_ratio = {}
    mutation_to_normalized_ratio = {}
    for line in lines:
        # strip off the carriage return
        line = line.strip()

        # skip empty lines and commented out lines
        if line != '' and not line.startswith( '#' ):
            # the second column can contain information
            mutation = line.split( ' ' )[0]
            try:
                ratio = line.split( ' ' )[1]
            except IndexError:
                ratio = None
                pass

            # add the mutation to the list
            all_mutations.append( mutation )
            # add the ratio, if any, to the mutation_to_ratio dict
            if ratio is not None:
                mutation_to_ratio[ mutation ] = float( ratio )
            else:
                mutation_to_ratio[ mutation ] = ratio

    '''
    # normalize the ratios given, if they were given
    # the min, max, and mean should ignore None's
    no_none_ratios = [ val for val in mutation_to_ratio.values() if val is not None ]
    if len( no_none_ratios ) == 0 or len( no_none_ratios ) == 1:
        mutat
    min_ratio = min( no_none_ratios )
    max_ratio = max( no_none_ratios )
    mean_ratio = sum( no_none_ratios ) / len( no_none_ratios )
    for mut_name, ratio in mutation_to_ratio.items():
        if ratio is not None:
            mutation_to_normalized_ratio[ mut_name ] = round( ( float( ratio ) - min_ratio ) / ( max_ratio - min_ratio ), 4 )
        else:
            mutation_to_normalized_ratio[ mut_name ] = None
    '''
    
    # attach the data to the DataHolder object
    data_holder.all_mutations = all_mutations
    data_holder.mutation_to_ratio = mutation_to_ratio
    #data_holder.mutation_to_normalized_ratio = mutation_to_normalized_ratio

    return data_holder





########################
#### DATA FUNCTIONS ####
########################

def get_full_contact_map( pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Returns a dictionary of each residue in <pose> that has a contact with another residue in the <pose> less than the given <cutoff> distance
    :param pose: Pose
    :param cutoff: int( or float( distance cutoff for what defines a contact). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: dict( key = residue num in pose, value = list of residue nums within <cutoff> distance )
    """
    if verbose:
        print "Getting the residue contact map for the Pose"

    # holds the resulting contact map
    contact_map = {}

    for seq_pos in range( 1, pose.n_residue() + 1 ):
        # holder for the residue numbers of residues within the given <cutoff> distance
        contacts = []

        # get the center for residue one
        center_1 = list( pose.residue( seq_pos ).nbr_atom_xyz() )

        # for every residue in the <pose>
        for residue in pose:
            seq_pos_2 = residue.seqpos()

            # if it's not the same residue (obviously they contact)
            if seq_pos != seq_pos_2:
                # get the center of the second residue
                center_2 = list( residue.nbr_atom_xyz() )

                # calculate the distance and store the residue number if it's below the given <cutoff>
                if calc_distance( center_1, center_2 ) < cutoff:
                    contacts.append( seq_pos_2 )

        # fill the contact_map with the contacts to that residue, only if there any
        if len( contacts ) != 0:
            contact_map[ seq_pos ] = contacts

    return contact_map



def get_contact_map_between_range1_range2( range1, range2, pose, cutoff = CUTOFF_DISTANCE, return_more_info = False, verbose = False ):
    """
    Returns a dictionary of each residue in <pose> that has a contact with another residue in the <pose> less than the given <cutoff> distance
    :param range1: list( residue pose numbers in range 1 )
    :param range2: list( residue pose numbers in range 2 )
    :param pose: Pose
    :param cutoff: int( or float( distance cutoff for what defines a contact). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param return_more_info: bool( do you want to return the number of contacts made as well? ) Default = False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: dict( key = residue num in pose, value = list of residue nums within <cutoff> distance )
    """
    if verbose:
        print "Getting the residue contact map for the Pose"

    # if possible, sort the two passed ranges
    try:
        range1.sort()
        range2.sort()
    except:
        pass

    # holds the resulting contact map
    contact_map = {}

    for seqpos_1 in range1:
        # holder for the residue numbers of residues within the given <cutoff> distance
        contacts = []

        # get the center for residue one
        center_1 = list( pose.residue( seqpos_1 ).nbr_atom_xyz() )

        # for every residue in the <pose>
        for seqpos_2 in range2:
            # get the center of the second residue
            center_2 = list( pose.residue( seqpos_2 ).nbr_atom_xyz() )

            # calculate the distance and store the residue number if it's below the given <cutoff>
            if calc_distance( center_1, center_2 ) < cutoff:
                contacts.append( seqpos_2 )

        # fill the contact_map with the contacts to that residue, only if there any
        if len( contacts ) != 0:
            contact_map[ seqpos_1 ] = contacts

    if return_more_info:
        num_contacts = sum( [ len( contacts ) for contacts in contact_map.values() ] )
        return contact_map, num_contacts
    else:
        return contact_map



def get_contact_map_with_JUMP_NUM( JUMP_NUM, pose, cutoff = CUTOFF_DISTANCE, return_more_info = False, verbose = False ):
    """
    Returns a dictionary of each residue in <pose> between the interface defined by <JUMP_NUM> that has a contact with another residue in the <pose> less than the given <cutoff> distance
    :param JUMP_NUM: int( JUMP number that defines the interface of interest )
    :param pose: Pose
    :param cutoff: int( or float( distance cutoff for what defines a contact). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param return_more_info: bool( do you want to return the number of contacts made as well? ) Default = False
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: dict( key = residue num in pose, value = list of residue nums within <cutoff> distance )
    """
    # get a list of all residue numbers on either side of the decoy's interface
    pose_r1 = []
    pose_r2 = []

    # find all of the residue numbers that correspond to side 1 and side 2
    for ii in range( 1, pose.total_residue() + 1 ):
        if ii < pose.fold_tree().downstream_jump_residue( JUMP_NUM ):
            pose_r1.append( ii )
        else:
            pose_r2.append( ii )

    # check to see that neither of the lists are empty
    if len( pose_r1 ) == 0 or len( pose_r2 ) == 0:
        print
        print "It appears that the jump number", JUMP_NUM, "does not actually define an interface. Exiting."
        sys.exit()

    # get the contact maps for the decoy and native poses
    if return_more_info:
        contact_map, num_contacts = get_contact_map_between_range1_range2( pose_r1, pose_r2, pose, cutoff = cutoff, return_more_info = return_more_info, verbose = verbose )
        return contact_map, num_contacts
    else:
        contact_map = get_contact_map_between_range1_range2( pose_r1, pose_r2, pose, cutoff = cutoff, return_more_info = return_more_info, verbose = verbose )
        return contact_map



def analyze_contact_map( contact_map, pose ):
    """
    Return the fraction of protein-protein, protein-carbohydrate, and carbohydrate_carbohydrate contacts found in the <pose>'s <contact_map> 
    :param contact_map: dict( key = residue num in pose, value = list of residue nums within a cutoff distance )
    :param pose: Pose( pose responsible for contact map )
    :return: data object that holds float( fraction pro-pro ), float( fraction pro-carb ), float( fraction carb-carb )
    """
    # instantiate contact type data counts
    pro_pro_contacts = 0
    pro_carb_contacts = 0
    carb_carb_contacts = 0

    # instantiate contact distances data list
    #contact_distances = []

    # instantiate carbohydrate to residue contact polarity data counts
    carb_to_polar_contacts = 0
    carb_to_nonpolar_contacts = 0
    carb_to_aromatic_contacts = 0
    
    # for each residue in contact map
    for resnum in contact_map.keys():
        # for each residue it makes contact with in the pose
        for contact_resnum in contact_map[ resnum ]:
            ## record the type of each contact made
            # protein to protein
            if pose.residue( resnum ).is_protein() and pose.residue( contact_resnum ).is_protein():
                pro_pro_contacts += 1
            # carbhydrate to carbohydrate
            elif pose.residue( resnum ).is_carbohydrate() and pose.residue( contact_resnum ).is_carbohydrate():
                carb_carb_contacts += 1
            # protein to carbohydrate
            else:
                pro_carb_contacts += 1

                # since this is a protein to carbohydrate contact, determine the polarity of the amino acid
                if pose.residue( resnum ).is_protein():
                    pro_resnum = resnum
                else:
                    pro_resnum = contact_resnum
                if pose.residue( pro_resnum ).is_polar():
                    carb_to_polar_contacts += 1
                if pose.residue( pro_resnum ).is_apolar():
                    carb_to_nonpolar_contacts += 1
                if pose.residue( pro_resnum ).is_aromatic():
                    carb_to_aromatic_contacts += 1

            ## record the contact distance using nbr atoms
            # get the atoms used as the center for these residues
            #resnum_nbr_xyz = list( pose.residue( resnum ).nbr_atom_xyz() )
            #contact_resnum_nbr_xyz = list( pose.residue( contact_resnum ).nbr_atom_xyz() )
            #contact_dist = calc_distance( resnum_nbr_xyz, contact_resnum_nbr_xyz )
            #contact_distances.append( contact_dist )

    # calculate fraction of contact types
    num_contacts = sum( [ len( contacts ) for contacts in contact_map.values() ] )
    try:
        pro_pro_fraction = round( float( pro_pro_contacts ) / float( num_contacts ), 2 )
    except ZeroDivisionError:
        pro_pro_fraction = "None"
    try:
        pro_carb_fraction = round( float( pro_carb_contacts ) / float( num_contacts ), 2 )
    except ZeroDivisionError:
        pro_carb_fraction = "None"
    try:
        carb_carb_fraction = round( float( carb_carb_contacts ) / float( num_contacts ), 2 )
    except ZeroDivisionError:
        carb_carb_fraction = "None"

    # calculate the fraction of carbohydrate to protein polarity contacts
    try:
        carb_to_polar_fraction = round( float( carb_to_polar_contacts ) / float( pro_carb_contacts ), 2 )
    except ZeroDivisionError:
        carb_to_polar_fraction = "None"
    try:
        carb_to_nonpolar_fraction = round( float( carb_to_nonpolar_contacts ) / float( pro_carb_contacts ), 2 )
    except ZeroDivisionError:
        carb_to_nonpolar_fraction = "None"
    try:
        carb_to_aromatic_fraction = round( float( carb_to_aromatic_contacts ) / float( pro_carb_contacts ), 2 )
    except ZeroDivisionError:
        carb_to_aromatic_fraction = "None"

    # add all data to data holder
    data_holder = DataHolder()
    data_holder.pro_pro_contacts = pro_pro_contacts
    data_holder.pro_carb_contacts = pro_carb_contacts
    data_holder.carb_carb_contacts = carb_carb_contacts
    data_holder.pro_pro_fraction = pro_pro_fraction
    data_holder.pro_carb_fraction = pro_carb_fraction
    data_holder.carb_carb_fraction = carb_carb_fraction
    data_holder.carb_to_polar_contacts = carb_to_polar_contacts
    data_holder.carb_to_nonpolar_contacts = carb_to_nonpolar_contacts
    data_holder.carb_to_aromatic_contacts = carb_to_aromatic_contacts
    data_holder.carb_to_polar_fraction = carb_to_polar_fraction
    data_holder.carb_to_nonpolar_fraction = carb_to_nonpolar_fraction
    data_holder.carb_to_aromatic_fraction = carb_to_aromatic_fraction
    #data_holder.contact_distances = contact_distances
    #data_holder.contact_distance_avg = calc_avg_of_list( contact_distances )
    #data_holder.contact_distance_max = max( contact_distances )
    #data_holder.contact_distance_min = min( contact_distances )

    return data_holder



def calc_Fnats_with_contact_maps( decoy_contact_map, decoy, native_contact_map, native, decoy_to_native_res_map = None ):
    """
    Return Fnat as calculated between two contact maps
    :param decoy_contact_map: dict( key = residue num in decoy, value = list of residue nums within a cutoff distance )
    :param decoy: Pose( decoy Pose responsible for decoy_contact_map )
    :param native_contact_map: dict( key = residue num in native, value = list of residue nums within a cutoff distance )
    :param native: Pose( native Pose responsible for native_contact_map )
    :param decoy_to_native_res_map: dict( if you're comparing a decoy to a native that has a different numbering scheme, pass in a dictionary that maps each resnum from the decoy to its corresponding resnum in the native )
    :return: DataHolder( object that holds floats of Fnat recovered, Fnat protein-protein recovered, Fnat protein-carbohydrate recovered, Fnat carbohydrate-carbohydrate recovered, Fnat carbohydrate to polar recovered, Fnat carbohydrate to nonpolar recovered, and Fnat carbohydrate to aromatic recovered for both protein-to-protein and protein-to-cabrohydate contacts )
    """
    # instantiate a DataHolder object
    data_holder = DataHolder()

    # get a DataHolder object with the data for the native
    native_data_holder = analyze_contact_map( native_contact_map, native )

    # instantiate contact type data counts
    pro_pro_contacts_recovered = 0
    pro_carb_contacts_recovered = 0
    carb_carb_contacts_recovered = 0

    # instantiate carbohydrate to residue contact polarity data counts
    carb_to_polar_contacts_recovered = 0
    carb_to_nonpolar_contacts_recovered = 0
    carb_to_aromatic_contacts_recovered = 0

    # get the number of contacts total found in the native
    num_native_contacts = sum( [ len( contacts ) for contacts in native_contact_map.values() ] )

    # determine the difference in contacts between decoy and native
    num_decoy_contacts_recovered = 0
    for decoy_resnum in decoy_contact_map.keys():
        # get the corresponding native residue number, if needed
        if decoy_to_native_res_map is not None:
            corresponding_native_resnum = decoy_to_native_res_map[ decoy_resnum ]
        # else the numbers are the same between the decoy and the native
        else:
            corresponding_native_resnum = decoy_resnum

        # for each contact this residue makes within decoy
        for decoy_contact_resnum in decoy_contact_map[ decoy_resnum ]:
            # get the corresponding native residue number for the contact, if needed
            if decoy_to_native_res_map is not None:
                corresponding_native_contact_resnum = decoy_to_native_res_map[ decoy_contact_resnum ]
            # else the numbers are the same between the decoy and the native
            else:
                corresponding_native_contact_resnum = decoy_contact_resnum

            # check to see if this is a contact made in the native
            try:
                native_contacts_at_this_decoy_resnum = native_contact_map[ corresponding_native_resnum ]

                # if the contact made in the decoy is the same as the contact made in the native, count it
                if corresponding_native_contact_resnum in native_contacts_at_this_decoy_resnum:
                    num_decoy_contacts_recovered += 1

                    # since this decoy contact is found in the native, analyze the contact for residue and polarity type
                    # protein to protein
                    if decoy.residue( decoy_resnum ).is_protein() and decoy.residue( decoy_contact_resnum ).is_protein():
                        pro_pro_contacts_recovered += 1
                    # carbhydrate to carbohydrate
                    elif decoy.residue( decoy_resnum ).is_carbohydrate() and decoy.residue( decoy_contact_resnum ).is_carbohydrate():
                        carb_carb_contacts_recovered += 1
                    # protein to carbohydrate
                    else:
                        pro_carb_contacts_recovered += 1

                        # since this is a protein to carbohydrate contact, determine the polarity of the amino acid
                        if decoy.residue( decoy_resnum ).is_protein():
                            decoy_pro_resnum = decoy_resnum
                        else:
                            decoy_pro_resnum = decoy_contact_resnum
                        if decoy.residue( decoy_pro_resnum ).is_polar():
                            carb_to_polar_contacts_recovered += 1
                        if decoy.residue( decoy_pro_resnum ).is_apolar():
                            carb_to_nonpolar_contacts_recovered += 1
                        if decoy.residue( decoy_pro_resnum ).is_aromatic():
                            carb_to_aromatic_contacts_recovered += 1

            # if this residue doesn't make contacts in the native pose, skip it
            except KeyError:
                pass

    # calculate Fnats
    try:
        Fnat_tot_contacts_recovered = round( ( float( num_decoy_contacts_recovered ) / float( num_native_contacts ) ) * 100, 2 )
    except ZeroDivisionError:
        Fnat_tot_contacts_recovered = "None"
    try:
        Fnat_pro_pro_contacts_recovered = round( ( float( pro_pro_contacts_recovered ) / float( native_data_holder.pro_pro_contacts ) ) * 100, 2 )
    except ZeroDivisionError:
        Fnat_pro_pro_contacts_recovered = "None"
    try:
        Fnat_pro_carb_contacts_recovered = round( ( float( pro_carb_contacts_recovered ) / float( native_data_holder.pro_carb_contacts ) ) * 100, 2 )
    except ZeroDivisionError:
        Fnat_pro_carb_contacts_recovered = "None"
    try:
        Fnat_carb_carb_contacts_recovered = round( ( float( carb_carb_contacts_recovered ) / float( native_data_holder.carb_carb_contacts ) ) * 100, 2 )
    except ZeroDivisionError:
        Fnat_carb_carb_contacts_recovered = "None"
    try:
        Fnat_carb_to_polar_contacts_recovered = round( ( float( carb_to_polar_contacts_recovered ) / float( native_data_holder.carb_to_polar_contacts ) ) * 100, 2 )
    except ZeroDivisionError:
        Fnat_carb_to_polar_contacts_recovered = "None"
    try:
        Fnat_carb_to_nonpolar_contacts_recovered = round( ( float( carb_to_nonpolar_contacts_recovered ) / float( native_data_holder.carb_to_nonpolar_contacts ) ) * 100, 2 )
    except ZeroDivisionError:
        Fnat_carb_to_nonpolar_contacts_recovered = "None"
    try:
        Fnat_carb_to_aromatic_contacts_recovered = round( ( float( carb_to_aromatic_contacts_recovered ) / float( native_data_holder.carb_to_aromatic_contacts ) ) * 100, 2 )
    except ZeroDivisionError:
        Fnat_carb_to_aromatic_contacts_recovered = "None"

    # attach data to DataHolder object
    data_holder.Fnat_tot_contacts_recovered = Fnat_tot_contacts_recovered
    data_holder.Fnat_pro_pro_contacts_recovered = Fnat_pro_pro_contacts_recovered
    data_holder.Fnat_pro_carb_contacts_recovered = Fnat_pro_carb_contacts_recovered
    data_holder.Fnat_carb_carb_contacts_recovered = Fnat_carb_carb_contacts_recovered
    data_holder.Fnat_carb_to_polar_contacts_recovered = Fnat_carb_to_polar_contacts_recovered
    data_holder.Fnat_carb_to_nonpolar_contacts_recovered = Fnat_carb_to_nonpolar_contacts_recovered
    data_holder.Fnat_carb_to_aromatic_contacts_recovered = Fnat_carb_to_aromatic_contacts_recovered

    return data_holder



def find_missing_native_contacts( decoy_contact_map, native_contact_map, native_to_decoy_res_map = None ):
    """
    Determine which contacts in the native are not found in the decoy
    :param decoy_contact_map: dict( key = residue num in decoy, value = list of residue nums within a cutoff distance )
    :param native_contact_map: dict( key = residue num in native, value = list of residue nums within a cutoff distance )
    :param decoy_to_native_res_map: dict( decoy numbering of contacts that are not made in decoy but found in native )
    """
    # determine the difference in contacts between decoy and native
    num_native_contacts_lost = 0
    native_contacts_lost_res_map = {}

    for native_resnum in native_contact_map.keys():
        # get the corresponding decoy residue number, if needed
        if native_to_decoy_res_map is not None:
            corresponding_decoy_resnum = native_to_decoy_res_map[ native_resnum ]
        # else the numbers are the same between the decoy and the native
        else:
            corresponding_decoy_resnum = native_resnum

        # for each contact this residue makes within native
        for native_contact_resnum in native_contact_map[ native_resnum ]:
            # get the corresponding decoy residue number for the contact, if needed
            if native_to_decoy_res_map is not None:
                corresponding_decoy_contact_resnum = native_to_decoy_res_map[ native_contact_resnum ]
            # else the numbers are the same between the decoy and the native
            else:
                corresponding_decoy_contact_resnum = native_contact_resnum

            # check to see if this is a contact made in the decoy
            try:
                decoy_contacts_at_this_native_resnum = decoy_contact_map[ corresponding_decoy_resnum ]

                # if the contact made in the native is lost in the decoy, count it
                if corresponding_decoy_contact_resnum not in decoy_contacts_at_this_native_resnum:
                    num_native_contacts_lost += 1

                    # add the lost contacts to the residue map
                    if corresponding_decoy_resnum in native_contacts_lost_res_map.keys():
                        if corresponding_decoy_contact_resnum not in native_contacts_lost_res_map[ corresponding_decoy_resnum ]:
                            native_contacts_lost_res_map[ corresponding_decoy_resnum ].append( corresponding_decoy_contact_resnum )
                    else:
                        native_contacts_lost_res_map[ corresponding_decoy_resnum ] = []
                        native_contacts_lost_res_map[ corresponding_decoy_resnum ].append( corresponding_decoy_contact_resnum )
            except:
                pass

    return native_contacts_lost_res_map, num_native_contacts_lost



def count_atomic_contacts_between_range1_range2( range1, range2, pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Counts the number of non-hydrogen atomic contacts in a <pose> between residues of <range1> and <range2> given the <cutoff> distance
    :param range1: list of ints( Pose residue numbers of side 1 )
    :param range2: list of ints( Pose residue numbers of side 2 )
    :param pose: Pose
    :param cutoff: int( or float( cutoff distance in Angstroms). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( number of atomic contacts at interface ), list( str( atm1 res name + atm name + atm1 chain + '_' + atm2 res name + atm2 atm name + atm2 chain
    """
    if verbose:
        print "Finding atomic contacts (except hydrogens) between", len( range1 ), "residues on side 1 and", len( range2 ), "residues on side 2"

    # instantiate counter and the list to hold the unique contact names
    contacts = 0
    contact_list = []
    
    # loop through each atom of each residue of side 1 specified and see if its within <cutoff> distance of any atom from side 2
    for res_num_1 in range1:
        # get the number of atoms in residue 1
        n_atms_1 = pose.residue( res_num_1 ).natoms()

        # loop over each atom
        for atm_1 in range( 1, n_atms_1 + 1 ):
            # if it's not a hydrogen atom, get its...
            if not pose.residue( res_num_1 ).atom_is_hydrogen( atm_1 ):
                # xyz coordinates
                atm1_xyz = list( pose.residue( res_num_1 ).atom( atm_1 ).xyz() )
                # residue name
                atm1_res_name = pose.residue( res_num_1 ).name3()
                # atom name
                atm1_atm_name = pose.residue( res_num_1 ).atom_name( atm_1 ).replace( ' ', '' )
                # and which chain it's on
                atm1_chain = pose.pdb_info().chain( res_num_1 )

                # do the same for residue 2
                for res_num_2 in range2:
                    # get the number of atoms in residue 2
                    n_atms_2 = pose.residue( res_num_2 ).natoms()

                    for atm_2 in range( 1, n_atms_2 + 1 ):
                        # if it's not a hydrogen atom, get its...
                        if not pose.residue( res_num_2 ).atom_is_hydrogen( atm_2 ):
                            # xyz coordinates
                            atm2_xyz = list( pose.residue( res_num_2 ).atom( atm_2 ).xyz() )
                            # residue name
                            atm2_res_name = pose.residue( res_num_2 ).name3()
                            # atom name
                            atm2_atm_name = pose.residue( res_num_2 ).atom_name( atm_2 ).replace( ' ', '' )
                            # and which chain it's on
                            atm2_chain = pose.pdb_info().chain( res_num_2 )

                            # if the two atoms are within the appropriate distance given by <cutoff>
                            atm1_atm2_dist = calc_distance( atm1_xyz, atm2_xyz )
                            if atm1_atm2_dist <= cutoff:
                                # makes unique names for each contact for ease of analysis
                                uniq_name_atm1 = atm1_res_name + '_' + atm1_atm_name + '_' + atm1_chain + '_' + str( res_num_1 )
                                uniq_name_atm2 = atm2_res_name + '_' + atm2_atm_name + '_' + atm2_chain + '_' + str( res_num_2 )

                                # unique name is the combination of the atom information from both atom1 and atom2
                                unique_name = uniq_name_atm1 + '.' + uniq_name_atm2

                                # if this specific contact has not already been found, add it to the list and up the contacts counter
                                if unique_name not in contact_list:
                                    contact_list.append( unique_name )
                                    contacts += 1

    return contacts, contact_list



def count_interface_atomic_contacts( JUMP_NUM, pose, cutoff = CUTOFF_DISTANCE, verbose = False ):
    """
    Counts the atom-to-atom (discluding hydrogens) contacts given a cutoff of <cutoff> Angstroms at the interface given the jump number <JUMP_NUM>
    Returns the number of contacts as a float and a list of the unique contacts made between the atoms of side 1 to side 2
    :param JUMP_NUM: int( Jump number that defines the interface )
    :param pose: Pose
    :param cutoff: int( or float( cutoff distance in Angstroms). Default = CUTOFF_DISTANCE = 5 Angstroms
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :return: float( number of atomic contacts at interface ), list( str( atm1 res name + atm name + atm1 chain + '_' + atm2 res name + atm2 atm name + atm2 chain
    """
    # check to see that the jump number is valid
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print "You gave me an invalid jump number, try again"
        sys.exit()

    if verbose:
        print "Counting interface contacts across jump number", JUMP_NUM, "discluding hydrogen atoms..."

    # get a list of all residue numbers from side 1 to 2
    side_1_list = []
    side_2_list = []

    # find all of the residue numbers that correspond to side 1 and side 2
    for ii in range( 1, pose.total_residue() + 1 ):
        if ii < pose.fold_tree().downstream_jump_residue( JUMP_NUM ):
            side_1_list.append( ii )
        else:
            side_2_list.append( ii )

    # check to see that neither of the lists are empty
    if len( side_1_list ) == 0 or len( side_2_list ) == 0:
        print
        print "It appears that the jump number", JUMP_NUM, "does not actually define an interface"
        sys.exit()

    # find all residues from side 1 and side 2 that are within 3 * CUTOFF_DISTANCE angstroms of each other
    # speeds up the counting process a bit
    side_1_list_close = []
    side_2_list_close = []

    # triple the given <cutoff> distance
    TRIPLE_CUTOFF_DISTANCE = cutoff * 3

    if verbose:
        print "Finding the residues within 3 times the given cutoff distance to speed up counting a bit"

    # loop over each residue from side 1
    for res_num_1 in side_1_list:
        # find the center of the residue in side 1
        center = pose.residue( res_num_1 ).nbr_atom_xyz()

        # now loop over the residues in side 2 for comparison
        for res_num_2 in side_2_list:
            # compare the distance to the other residue in side 2  -  using triple the distance at first to theoretically speed up the counting process
            if center.distance( pose.residue( res_num_2 ).nbr_atom_xyz() ) <= TRIPLE_CUTOFF_DISTANCE:
                # if the residues are within triple the <cutoff>, and if they haven't already been added, add them to the close list
                if res_num_1 not in side_1_list_close:
                    side_1_list_close.append( res_num_1 )

                if res_num_2 not in side_2_list_close:
                    side_2_list_close.append( res_num_2 )

    # loops over side 1 list and see which atoms are within 5 Ang of side 2 for native pose
    num_interface_contacts, interface_contact_list = count_atomic_contacts_between_range1_range2( side_1_list_close, side_2_list_close, pose, cutoff )

    return float( num_interface_contacts ), interface_contact_list



def calc_interface_sasa( pose, JUMP_NUM ):
    """
    Use rosetta.calc_total_sasa to compute the SASA of the total pose - SASA of the split-apart pose at <JUMP_NUM>
    :param pose: Pose
    :param JUMP_NUM: int( valid Jump number defining interface )
    :return: float( interface_SASA value )
    """
    # imports
    from rosetta import Pose, calc_total_sasa
    from rosetta.numeric import xyzVector_Real

    
    # make sure a valid JUMP_NUM was passed in
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print JUMP_NUM, "is an invalid Jump number. Exiting"
        sys.exit()
    
    # make and split the temporary pose
    temp_pose = Pose()
    temp_pose.assign( pose )
    jump = temp_pose.jump( JUMP_NUM )

    #TODO-get current xyz location and multiply by 500 or something instead
    vec = xyzVector_Real( 1000, 1000, 1000 )
    jump.set_translation( vec )
    temp_pose.set_jump( JUMP_NUM, jump )
    
    # calculate the SASA values
    total_sasa = calc_total_sasa( pose, PROBE_RADIUS )
    split_sasa = calc_total_sasa( temp_pose, PROBE_RADIUS )
    delta_interface_sasa = total_sasa - split_sasa
    
    return delta_interface_sasa



def analyze_interface( pose, JUMP_NUM, pack_separated = True, verbose = False ):
    """
    Use rosetta.protocols.analysis.Interface Analyzer to compute various interface metrics
    :param pose: Pose
    :param JUMP_NUM: int( valid Jump number defining interface )
    :param pack_separated: bool( Do you want to pack the protein after you split them apart? ). Default = True
    :param verbose: bool( if you want the function to print out the high-energy residues, set to True ). Default = False
    :return: InterfaceAnalyzerMover with data ready to be extracted from it
    """
    from rosetta.protocols.analysis import InterfaceAnalyzerMover as IAM


    # make sure a valid JUMP_NUM was passed in
    if JUMP_NUM <= 0 or JUMP_NUM > pose.num_jump():
        print
        print JUMP_NUM, "is an invalid Jump number. Exiting"
        sys.exit()

    # instantiate an InterfaceAnalyzer mover
    IAmover = IAM()

    # set the interface jump to the jump number passed in
    IAmover.set_interface_jump( JUMP_NUM )

    # TODO-see what's worth calculating and then actually return the value
    # set what you want to calculate
    IAmover.set_compute_interface_delta_hbond_unsat( True )
    IAmover.set_compute_interface_sc( True )
    IAmover.set_compute_separated_sasa( True )

    # set the relevant options
    IAmover.set_input_pose( pose )
    IAmover.set_pack_separated( pack_separated )

    if verbose:
        print "Analyzing interface..."
    IAmover.reset_status()
    IAmover.apply( pose )

    # retrieve data
    # TODO-retrieve new, relevant data
    #interface_dSASA = IAmover.get_interface_delta_sasa()
    ##unsat_hbond =  IAmover.get_interface_delta_hbond_unsat()
    #interface_dG = IAmover.get_interface_dG()  # I already calculate this with the fxn I wrote
    ##num_interface_residues = IAmover.get_num_interface_residues()

    return IAmover



def get_interface_score( JUMP_NUM, sf, pose, watch = False ):
    """
    Given a jump number that defines the interface, calculates Rosetta's ddG interface
    Splits apart the two domains defined by the <JUMP_NUM>, scores it, then subtracts that from the total score of the <pose>  -  result is the interface score
    :param JUMP_NUM: int( valid Jump number of the interface )
    :param sf: ScoreFunction
    :param pose: Pose
    :param watch: bool( do you want to debug this in PyMol? ) Default = False
    :return: float( ddG interface score )
    """
    # imports
    from rosetta import Pose
    from rosetta.numeric import xyzVector_Real
    if watch:
        from rosetta import PyMOL_Mover
        pmm = PyMOL_Mover()


    # get start score
    start_score = sf( pose )

    # make and split the temporary pose
    temp_pose = Pose()
    temp_pose.assign( pose )
    jump = temp_pose.jump( JUMP_NUM )
    if watch:
        temp_pose.pdb_info().name( "temp_pose" )
        pmm.apply( temp_pose )

    # multiply the jump's translation vector by 500
    vec = jump.get_translation() * 500
    jump.set_translation( vec )
    temp_pose.set_jump( JUMP_NUM, jump )
    if watch:
        temp_pose.pdb_info().name( "split_pose" )
        pmm.apply( temp_pose )

    # get and return interface score
    split_apart_score = sf( temp_pose )
    interface_score = start_score - split_apart_score

    return interface_score



def get_phi_psi_omega_of_res( pose, seqpos ):
    """
    Returns the phi, psi, and omega values of the passed residue <seqpos> in <pose>
    :param pose: Pose
    :param seqpos: int( the sequence position of the residue of interest )
    :return: float( phi ), float( psi ), float( omega )
    """
    try:
        return pose.phi( seqpos ), pose.psi( seqpos ), pose.omega( seqpos )
    except:
        print "%s does not seem to be a residue in the passed Pose. Check your input" %seqpos
        return None
        
    

def determine_amino_acid_composition( pose ):
    """
    Prints the amino acid composition of the <pose> using Pandas DataFrame
    Ignores sugars ( single letter code 'Z' ) in composition analysis
    :param pose: Pose
    :return: a DataFrame or dictionary of the data, depending on if Pandas can be imported
    """
    try:
        import pandas as pd
        pandas_on = True
    except ImportError:
        pandas_on = False
        print "Skipping Pandas import - consider downloading it! Who doesn't love Pandas??"
        pass


    # initialize a dictionary for the data
    # won't be used if Pandas import worked
    data_dict = {}

    # build lists for pandas DataFrame
    res_count = []
    percentage = []

    # get number of amino acid residues
    res_total = ( pose.total_residue() - pose.sequence().count( 'Z' ) )

    for res_type in AA_name1_list:
        count = float( pose.sequence().count( res_type ) )
        res_count.append( count )
        percentage.append( round( ( count / res_total * 100 ), 2 ) )

    # if global Pandas import was successful, create a DataFrame
    if pandas_on:
        df_AA_composition = pd.DataFrame()
        df_AA_composition["Res_Type"] = AA_name1_list
        df_AA_composition["Count"] = res_count
        df_AA_composition["Percentage"] = percentage
        
    # else create a dictionary
    else:
        data_dict["Res_Type"] = AA_name1_list
        data_dict["Count"] = res_count
        data_dict["Percentage"] = percentage
        
    # print out relevant data
    print
    print "This protein has..."
    if pandas_on:
        print df_AA_composition.sort( "Percentage", ascending = False )
        return df_AA_composition
    else:
        print data_dict
        return data_dict



def compare_all_rotamers( pose1, pose2, E_diff = 5 ):
    """
    Uses is_similar_rotamer on each residue, only if they are the same amino acids
    Will exit if the two Poses don't have the same number of amino acids - can't guarantee the validity of the comparison otherwise
    Skips over sugars and glycines
    :param pose1: Pose 1
    :param pose2: Pose 2
    :param E_diff: int( or float( the energy difference required to be considered different ). Default = +/- 5
    :return: a DataFrame or dictionary of the data, depending on if Pandas can be imported (does not print)
    """
    # ensure that the two poses have the same number of residues
    if not pose1.total_residue() == pose2.total_residue():
        print
        print "Check your poses, they don't have the same number of residues. Exiting"
        sys.exit()

    # get the absolute value of the <E_diff>, and also get its negation
    E_diff = abs( E_diff )
    neg_E_diff = -1 * E_diff

    # make the lists to hold data
    resi_num_pose = []
    resi_chain = []
    resi_name = []

    # compare each residue rotamer
    total_residues = pose1.total_residue()
    for ii in range( 1, total_residues + 1 ):
        res1_AA = pose1.residue( ii ).name1()
        res2_AA = pose2.residue( ii ).name1()

        # only compare if the two residues are the same amino acid, and are not sugars (I'm assuming sugars will change enough)
        if res1_AA == res2_AA:
            res1 = pose1.residue( ii )
            res2 = pose2.residue( ii )

            # if they aren't sugars
            if not res1.is_carbohydrate():
                if not res1.name1() != 'G':
                    if not res2.is_carbohydrate():
                        if not res2.name1() != 'G':
                            for chi in range( 1, len( res1.chi() ) + 1 ):
                                diff = res2.chi()[ chi ] - res1.chi()[ chi ]
                                if diff >= E_diff or diff <= neg_E_diff:
                                    resi_num_pose.append( ii )
                                    # get the chain of the residue by going from pose --> PDB numbering and getting it from there
                                    resi_chain.append( pose1.pdb_info().pose2pdb( ii )[ -2:len( pose1.pdb_info().pose2pdb( ii ) ) ] )
                                    resi_name.append( res1.name1() )

    # make and return Pandas DataFrame if Pandas was imported
    if pandas_on:
        df = pd.DataFrame()
        df["Pose Num"] = resi_num_pose
        df["Chain"] = resi_chain
        df["AA"] = resi_name
        
        return df
    
    # if Pandas wasn't imported, make and return a dictionary of the data
    else:
        # initialize a dictionary for the data
        data_dict = {}
        
        # append data
        data_dict["Pose Num"] = resi_num_pose
        data_dict["Chain"] = resi_chain
        data_dict["AA"] = resi_name
        
        return data_dict



def check_E_per_residue( sf, pose, energy_cutoff = 1.5, verbose = False ):
    """
    Returns a list of the residues with a total energy value greater than <energy_cutoff>
    :param sf: ScoreFunction
    :param pose: Pose
    :param energy_cutoff: int( or float( residue energy cutoff where residues of interest have a greater energy value ) ). Default = 1.5
    :param verbose: bool( if you want the function to print out the high-energy residues, set to True ). Default = False
    :return: list( Residue objects of residues with total energy above <energy_cutoff> )
    """
    # get the absolute value of the passed energy_cutoff
    abs_cutoff = abs( energy_cutoff )
    
    # score the pose to get access to its energy data
    sf( pose )
    
    # iterate through each residue in pose looking at its energy
    high_E_residues = []
    for residue in pose:
        energy = pose.energies().residue_total_energy( residue.seqpos() )
        if energy > abs_cutoff:
            if verbose:
                print residue.name(), '\t', pose.pdb_info().pose2pdb( residue.seqpos() ), '\t', energy
            high_E_residues.append( residue )
            
    return high_E_residues


def compare_native_vs_decoy_per_res( sf, native, decoy ):
    """
    Uses Pandas to construct a DataFrame of residues between the two poses that show a difference in scores
    :param sf: ScoreFunction
    :param native: Pose
    :param decoy: Pose
    :return: DataFrame
    """
    # imports
    import pandas as pd
    from rosetta import score_type_from_name
    from rosetta.core.scoring.sasa import per_res_sc_sasa


    # score the poses to get access to their energy objects
    sf( native )
    sf( decoy )

    # absolute per residue sidechain SASA
    native_per_res_sc_sasa = per_res_sc_sasa( native )
    decoy_per_res_sc_sasa = per_res_sc_sasa( decoy )

    # instantiate the DataFrame object
    df = pd.DataFrame()
    residue_numbers = range( 1, native.n_residue() + 1 )
    df[ "res_num" ] = residue_numbers

    # get delta scores between native and the two passed decoys
    # decoy
    df[ "dfa_atr" ] = [ decoy.energies().residue_total_energies( res ).get( score_type_from_name( "fa_atr" ) ) - 
                              native.energies().residue_total_energies( res ).get(score_type_from_name( "fa_atr" ) ) for res in residue_numbers ]
    df[ "dfa_rep" ] = [ decoy.energies().residue_total_energies( res ).get( score_type_from_name( "fa_rep" ) ) - 
                              native.energies().residue_total_energies( res ).get(score_type_from_name( "fa_rep" ) ) for res in residue_numbers ]
    df[ "dfa_sol" ] = [ decoy.energies().residue_total_energies( res ).get( score_type_from_name( "fa_sol" ) ) - 
                              native.energies().residue_total_energies( res ).get(score_type_from_name( "fa_sol" ) ) for res in residue_numbers ]
    df[ "per_res_dsasa" ] = [ decoy_per_res_sc_sasa[ ii ] - native_per_res_sc_sasa[ ii ] for ii in residue_numbers ]

    # set the residue numbers as the index
    # there has to be a smarter way of doing this
    df = df.set_index( "res_num" )

    # return a DataFrame that excludes rows that are all zero
    nonzero_df = df[(df.T != 0).any()]

    return nonzero_df



def get_pymol_res_selection( residues, input_pose ):
    """
    Using the list of Pose <residues> numbers, turn them into a string selection using pdb_info() from <input_pose>
    :param residues: list( Pose residue numbers )
    :param input_pose: Pose
    :return: str( pymol selection of residues )
    """
    pymol_selection = "select resi "
    for ii in range( len( residues ) ):
        # for all residues except the last one
        if ii != len( residues ) - 1:
            pymol_selection += input_pose.pdb_info().pose2pdb( residues[ii] ).strip().split( ' ' )[0] + " and chain " + input_pose.pdb_info().pose2pdb( residues[ii] ).split( ' ' )[1] + " + resi "
        else:
            pymol_selection += input_pose.pdb_info().pose2pdb( residues[ii] ).strip().split( ' ' )[0] + " and chain " + input_pose.pdb_info().pose2pdb( residues[ii] ).split( ' ' )[1]

    return pymol_selection



def get_rank_order_of_list( input ):
    """
    Function retreived from http://codereview.stackexchange.com/questions/65031/creating-a-list-containing-the-rank-of-the-elements-in-the-original-list
    """
    indices = list(range(len(input)))
    indices.sort(key=lambda x: input[x])
    output = [0] * len(indices)
    for i, x in enumerate(indices):
        # i'm adding a plus one because this ranks from 0 - n-1, I want 1 - n
        # also added the float() part
        output[x] = float( i + 1 )
    return output



def make_RotamerTrialsMover( moveable_residues, sf, input_pose, pack_radius = None ):
    """
    Given a list of <moveable_residues>, get all additional residues within <pack_radius> Angstroms around them in <input_pose> and return a RotamerTrialsMover that will pack these residues
    If no <pack_radius> is given, then only residues in <moveable_residues> will be allowed to pack
    :param moveable_residues: list( Pose numbers )
    :param sf: ScoreFunction
    :param input_pose: Pose
    :param pack_radius: int( or float( radius in Angstroms to pack around the <moveable_residues>. Uses nbr_atom to determine residues in proximity ) ) Default = None which means that only residues in <moveable_residues> get packed
    :return: RotamerTrialsMover
    """
    # imports
    from rosetta import standard_packer_task, RotamerTrialsMover


    # copy over the input_pose
    pose = input_pose.clone()

    # make the PackRotamersMover from the passed MoveMap
    # default of standard_packer_task is to set packing for residues to True
    task = standard_packer_task( pose )
    task.or_include_current( True )
    task.restrict_to_repacking()

    # if a pack_radius was not given, then everything gets packed. So the task does not need to be adjusted as the default option is packing True for all
    # otherwise, if a pack_radius was given, turn off repacking for residues outside the pack_radius
    if pack_radius is not None:
        # get all the protein residues within pack_radius of the moveable_residues
        # inclue_passed_res_nums means that the function will return a list of residues that includes all numbers in moveable_residues
        # I am adding them in myself for clarity, so this setting is set to off
        nearby_protein_residues = get_res_nums_within_radius( moveable_residues, pose, 
                                                              radius = pack_radius, 
                                                              include_passed_res_nums = False )

        # create a list of residue numbers that can be packed
        # meaning, the moveable carbohydrate residues and the residues around them
        packable_residues = [ res_num for res_num in moveable_residues ]
        packable_residues.extend( nearby_protein_residues )
        packable_residues = list( set( packable_residues ) )

        # turn off packing for all residues that are NOT packable
        # i.e. for all residues in the pose, turn OFF packing if they are NOT in the packable_residues list
        [ task.nonconst_residue_task( res_num ).prevent_repacking() for res_num in range( 1, pose.n_residue() + 1 ) if res_num not in packable_residues ]

    # otherwise, only residues specified by moveable_residues should be allowed to be packed
    else:
        # turn off repacking for all residues in the pose that are NOT in moveable_residues
        # no pack_radius was given, so all residues not specified in moveable_residues should not be packed
        [ task.nonconst_residue_task( res_num ).prevent_repacking() for res_num in range( 1, pose.n_residue() + 1 ) if res_num not in moveable_residues ]

    # make the pack_rotamers_mover with the given ScoreFunction and created task
    pack_rotamers_mover = RotamerTrialsMover( sf, task )

    return pack_rotamers_mover





############################
#### INITIALIZE ROSETTA ####
############################

if __name__ == '__main__':
    # initialize rosetta with sugar flags
    from rosetta import init
    
    print "Initializing Rosetta with sugar flags"
    #init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records" )
    init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -write_pdb_link_records" )
    #init( extra_options="-include_sugars -override_rsd_type_limit -write_pdb_link_records -constant_seed" )
    #init( extra_options="-include_sugars -override_rsd_type_limit -write_pdb_link_records -constant_seed -out:level 400" )

############################
#### INITIALIZE ROSETTA ####
############################
