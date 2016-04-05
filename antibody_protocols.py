__author__ = 'morganlnance'


#####################################
# HOME OF USEFUL PROTOCOL FUNCTIONS #
#####################################

# import needed Rosetta functions
from rosetta import MonteCarlo, SmallMover, ShearMover
from rosetta.protocols.rigid import RigidBodyPerturbMover
import random

# TODO -- add "watch" possibility for each protocol (take in an arg for a pymol mover)


# all relevant imports and function calls should be made possible after importing from mutational_analysis.py
from antibody_functions import *


def make_base_pack_min_pose( sf, pose, trials = 2, dump_best_pose = False, dump_pose_name = None, dump_dir = None, verbose = False, pmm = None):
    """
    Take the given <pose>, pack and minimize it once to get a base, then get a low E structure using get_best_structure_based_on_pose
    :param sf: ScoreFunction                   
    :param pose: Pose                          
    :param trials: int( number of times to pack and minimize the pose ). Default = 2
    :param dump_best_pose: bool( after finding the lowest energy Pose, dump the structure into the current directory (or to <dump_dir>). Default = False
    :param dump_pose_name: bool( filename of the Pose to be dumped ). Default = None (ie. "Best_" + current Pose name)
    :param dump_dir: str( path/to/dir/where/pose/will/be/dumped ). Default = None (ie. current working directory)
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: Pose( the Pose with the lowest total score after the trials of packing and minimization )
    """
    working_pose = Pose()
    working_pose.assign( pose )
    
    for ii in range( trials ):
        # pack
        pack_rotamers_mover = make_pack_rotamers_mover( sf,
                                                        working_pose,
                                                        pack_branch_points = True,
                                                        residue_range = None,
                                                        use_pack_radius = False,
                                                        verbose = False )
        pack_rotamers_mover.apply( working_pose )
        
        # minimize
        min_mover = make_min_mover( sf,
                                    working_pose,
                                    # jumps = None actually means minimize all jumps
                                    jumps = None,
                                    allow_sugar_chi = False,
                                    verbose = False )

    return working_pose



def make_rigid_body_moves( sf, pose, loops_file = None, anchor_loops = True, temperature = kT, rotation = 4, translation = 0.6, trials = 100, ramp_elec_score = False, adjust_mc_weights = True, dont_perturb_these_jumps = None, verbose = False, pmm = None ):
    """
    Make rigid body movements in <pose> of the domains defined by the <loops_file> (if desired) using the ScoreFunction <sf> over <trials>
    Iterates through each jump, including jumps created by the <loops_file> and samples rigid body movements
    If you don't want the MonteCarlo translation and rotation weights to change during the trials, set <adjust_mc_weights> to False
    :param sf: ScoreFunction
    :param pose: Pose
    :param loops_file: str( a LOOPs file that defines the Pose's domains through the Loops ). Default = None
    :param anchor_loops: bool( do you want to add +2-residue anchors to the ends of your loops? ). Default = True
    :param temperature: int( or float( temperature used for MonteCarlo object ). Default = kT = 0.7
    :param rotation: int( or float( rotational freedom allowed in degrees ). Default = 4
    :param translation: int( or float( translational freedom allowed in Angstroms ). Default = 0.6
    :param trials: int( number of motion trials to be completed ). Default = 100 (minimum of 10)
    :param ramp_elec_score: bool( do you want to do a high-to-low ramping of fa_elec during the trials? ). Default = False
    :param adjust_mc_weights: bool( adjust the MonteCarlo's translation and rotation freedom during trial? ). Default = True
    :param dont_perturb_these_jumps: list( a list of valid int( Jump numbers ) of any jumps you don't want to perturb )
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: the <pose> after making the rigid body movements
    """
    # make a list of numbers of 1 to the total number of jumps in pose
    # will remove jumps specified in dont_perturb_these_jumps by using .remove()
    jump_numbers = range( 1, pose.num_jump() + 1 )

    # check the validity of the dont_perturb_these_jumps parameter, if not None
    if dont_perturb_these_jumps is not None:
        # first, check to make sure it's a list
        if not isinstance( dont_perturb_these_jumps, list):
            print "You didn't give me a list for the dont_perturb_these_jumps parameter:", dont_perturb_these_jumps, "I need a list of int( Jump numbers ) to make this work. Exiting"
            sys.exit()

        # second, check to make sure that the list contains valid jump numbers
        num_jumps = pose.num_jump()
        for jump in dont_perturb_these_jumps:
            if isinstance( jump, int ):
                # ensure the Jump number is within the range of the number of Jumps the pose has
                if jump <= 0 or jump > num_jumps:
                    print "One of the values in your list of Jumps for the dont_perturb_these_jumps parameter isn't a valid Jump number:", jump, "Please ensure the validity of the Jump numbers. Returning None"
                    return None
            # if the Jump number isn't an integer, return None and print an error message
            else:
                print "One of the values in your list of Jumps for the dont_perturb_these_jumps paramter isn't an integer:", jump, "I need an integer specifically. Returning None"
                return None

    # now that the validity of the jump numbers have been checked, remove the jump numbers
    if dont_perturb_these_jumps is not None:
        for jump_num in dont_perturb_these_jumps:
            jump_numbers.remove( jump_num )
            
    # get the Loops from the <loops_file>
    if loops_file is not None:
        loops = make_loops_from_file( loops_file, pose, anchor_loops = anchor_loops )
        
        # if there are 0 loops, inform user and return None
        if loops.num_loop() == 0:
            print "The LOOPs file passed resulting in 0 Loops - check your input/LOOPs file. Returning None."
            return None
    # otherwise make an empty Loops object
    else:
        loops = Loops()

    # change the fold tree to the necessary one allowing domain motion given the <loops_file>
    if loops_file is not None:
        pose = setup_new_fold_tree( loops_file, pose, anchor_loops = anchor_loops, add_cutpoints = False )

    # add chainbreak term to sf for use in loop scoring
    # rigid body motion allows for loops to break - adding chainbreak will capture the energetic cost of this
    sf.set_weight( score_type_from_name( "chainbreak" ), 1 )

    # apply sugar branch point constraints to sf
    print "Constraining distances and angles at sugar branch points"
    apply_sugar_constraints_to_sf( sf, pose )

    # if no problem with FoldTree creation, load the needed Rosetta modules
    print "Making rigid body moves..."
    from rosetta.protocols.rigid import RigidBodyPerturbMover

    # store the original fa_atr, fa_rep, and fa_elec weights
    FA_ATR_ORIG = sf.get_weight( score_type_from_name( "fa_atr" ) )
    FA_REP_ORIG = sf.get_weight( score_type_from_name( "fa_rep" ) )
    FA_ELEC_ORIG = sf.get_weight( score_type_from_name( "fa_elec" ) )

    # raise the fa_atr term and lower the fa_rep term in the ScoreFunction to be able to do a ramping of its contribution
    FA_ATR_NEW = FA_ATR_ORIG * 2
    sf.set_weight( score_type_from_name( "fa_atr" ), FA_ATR_NEW )

    FA_REP_NEW = FA_REP_ORIG * 0.2
    sf.set_weight( score_type_from_name( "fa_rep" ), FA_REP_NEW )
    
    # if desired, raise the fa_elec term to allow ramping
    if ramp_elec_score:
        FA_ELEC_NEW = FA_ELEC_ORIG * 2
        sf.set_weight( score_type_from_name( "fa_elec" ), FA_ELEC_NEW )
        

    # TODO-fix centroid mode for glycosylated proteins
    # switch the pose from full atom to centroid mode
    ##from rosetta import SwitchResidueTypeSetMover
    ##switch = SwitchResidueTypeSetMover( "centroid" )
    ##switch.apply( pose )

    # create the MonteCarlo object and attach the pose, ScoreFunction with the adjusted weights, and a temperature
    mc = MonteCarlo( pose, sf, temperature )
    num_accepts = 0
    num_rejects = 0

    # over a given number of trials, make rigid body moves and accept or reject them using the MonteCarlo object
    # start temperature, rotation, and translation are used to revert to original weights if too many rejects have happened
    working_temperature = temperature
    start_temperature = temperature
    working_rotation = rotation
    start_rotation = rotation
    working_translation = translation
    start_translation = translation

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
        
    # the jump number counter starts at 0 instead of 1 because we're indexing from a normal list
    cur_jump_index = 0
    
    # inform the user of options passed in
    if adjust_mc_weights:
        print "Allowing the adjustment of temperature, rotational freedom, and translation freedom during protocol for MonteCarlo use"
    print "Starting temperature:", start_temperature, "rotational freedom:", start_rotation, "translational freedom:", start_translation

    for ii in range( 1, trials + 1 ):
        # ramp up or down the appropriate scoring terms and get it back to the MonteCarlo object
        sf = ramp_score_weight( sf, "fa_atr", FA_ATR_ORIG, ii - 1, trials )
        sf = ramp_score_weight( sf, "fa_rep", FA_REP_ORIG, ii - 1, trials )
        sf = ramp_score_weight( sf, "fa_elec", FA_ELEC_ORIG, ii - 1, trials )
        mc.score_function( sf )

        # every 10 trials...
        if ii % 10 == 0 and ii != 0:
            # calculate and print out the acceptance rate
            acceptance_rate = float( num_accepts ) / float( mc.total_trials() )
            print "  Acceptance rate:", str( round( acceptance_rate * 100, 3 ) ) + '%', "after", ii, "trials of", trials

            # if the user allows for the weights to be adjusted during the rigid move sampling...
            if adjust_mc_weights:
                # if the acceptance rate is too low, adjust temp, rotation, or translation accordingly
                if acceptance_rate < .495:
                    working_temperature *= 1.15  # increase temperature and give it to the MonteCarlo object
                    mc.set_temperature( working_temperature )
                    working_rotation *= 0.85  # decrease rotation
                    working_translation *= 0.85  # decrease translation
                    print "    * Acceptance rate too low, changing weights"
                    print "    temp:", round( working_temperature, 3 ), "rotation:", round( working_rotation, 3 ), "translation:", round( working_translation, 3 )

                    # clear the number of rejected moves
                    num_rejects = 0

                # if the acceptance rate is pretty high, adjust temp, rotation, or translation accordingly
                if 0.505 < acceptance_rate <= 1.0:
                    working_temperature *= 0.85  # decrease temperature and give it to the MonteCarlo object
                    mc.set_temperature( working_temperature )
                    working_rotation *= 1.15  # increase rotation
                    working_translation *= 1.1  # slightly increase translation  -  this isn't very realistic, so keep it low
                    print "    * Acceptance rate is high, changing weights"
                    print "    temp:", round( working_temperature, 3 ), "rotation:", round( working_rotation, 3 ), "translation:", round( working_translation, 3 )

                    # clear the number of rejected moves
                    num_rejects = 0

        # if the past 15 trials ( 5 full rounds with 3 MonteCarlo's each round ) have been rejected, reset the original temperature, rotation, and translation
        if num_rejects >= 15:
            print "* A critical number of rejections has been hit in a row - resetting weights to starting values"
            working_temperature = start_temperature
            mc.set_temperature( working_temperature )
            working_rotation = start_rotation
            working_translation = start_translation
            num_rejects = 0

        # perturb the domains defined by the <loops_file> by looping over each Jump in the <pose>
        jump_num = jump_numbers[ cur_jump_index ]
        pert_mover = RigidBodyPerturbMover( jump_num, working_rotation, working_translation )
        pert_mover.apply( pose )
        if pmm is not None and pmm_worked:
            pmm.apply( pose )
            
        # minimize the now moved domain
        mm = make_movemap_for_jumps( jump_num )
        do_min_with_this_mm( mm, sf, pose )
        
        # get the stop residue number of the current Jump number
        stop_res_of_current_jump = get_seq_pos_from_jump_num( jump_num, pose )
        
        # iterate through the stop residue numbers of each Loop made by the <loops_file>
        for loop_num in range( 1, loops.num_loop() + 1 ):
            loop = loops[ loop_num ]
            stop = loop.stop()
            
            # if this Jump is from a Loop made from the <loops_file>, minimize that Loop
            if stop == stop_res_of_current_jump:
                mm = make_movemap_for_loop( loop )
                do_min_with_this_mm( mm, sf, pose )
                if pmm is not None and pmm_worked:
                    pmm.apply( pose )
                    
        # accept or reject movement and adjust counters accordingly
        if mc.boltzmann( pose ):
            num_accepts += 1
            num_rejects = 0
        else:
            num_rejects += 1
        if pmm is not None and pmm_worked:
            pmm.apply( pose )
            
        # update the jump_num index
        if cur_jump_index == len( jump_numbers ) - 1:
            cur_jump_index = 0
        else:
            cur_jump_index += 1
        
        
    if pmm is not None and pmm_worked:
        pmm.apply( pose )
    
    # switch the pose back to full atom mode
    ##switch = SwitchResidueTypeSetMover( "fa_standard" )
    ##switch.apply( pose )

    # wrap up with pack and minimization
    print "Running a pack and minimization"
    pose = do_pack_min( sf, pose )
    
    # change the FoldTree back to the original
    if loops_file is not None:
        restore_original_fold_tree( pose )

    return pose



def make_loop_perturbations( loops_file, sf, pose, temperature = kT, trials = 100, PDB_numbering = False, anchor_loops = True, model_loops = True, verbose = False, pmm = None ):
    """
    Adds Loops to <pose> by reading <loops_file> and then, over <trials>, makes Small and Shear moves on loop then scores them using ScoreFunction <sf>
    <model_loops> = True, means that the Loops will be broken then closed after perturbation - used for modeling loop movement
    <model_loops> = False, means that the Loops won't break after perturbation - used more for sampling linker movement
    IMPORTANT: Will not work if Loop includes sugars
    :param loops_file: str( /path/to/LOOPs/file )
    :param sf: ScoreFunction
    :param pose: Pose
    :param temperature: int( or float( temperature used for MonteCarlo object and Small and Shear movers ). Default = kT = 0.7
    :param trials: int( number of trials to be completed ). Default = 100. Minimum of 10 trials
    :param PDB_numbering: bool( is this file numbered using PDB numbering? ) Default = False
    :param anchor_loops: bool( do you want to add +2-residue anchors to the ends of your loops? ). Default = True
    :param model_loops: bool( apply a FoldTree with the Loops to the Pose to model the Loops? ). Default = True
    :param verbose: bool( if you want the function to print out statements about what its doing, set to True ). Default = False
    :param pmm: PyMOL_Mover( pass a PyMOL_Mover object if you want to watch the protocol ). Default = None
    :return: the <pose> after Small and Shear moves have been made
    """
    # get the Loops from the <loops_file>
    # this doesn't need to be verbose because this is just to make sure that the loops file actually makes loops
    loops = make_loops_from_file( loops_file, pose, PDB_numbering = PDB_numbering )

    # LoopModelerMover!!!!
    
    # ensure that valid Loops were made by checking the number of Loops
    if loops.num_loop() == 0:
        print "It seems the loops_file that you gave me did not make any Loops - check your input. Returning None"
        return None

    # make a FoldTree with the Loops and give it to the Pose, if <model_loops> is True
    if model_loops:
        # new FoldTree for Pose with cutpoints from <loops_file>
        pose = setup_new_fold_tree( loops_file, pose, PDB_numbering = PDB_numbering, anchor_loops = anchor_loops, verbose = verbose )

    # apply sugar branch point constraints to sf
    print "Constraining distances and angles at sugar branch points"
    sf = apply_sugar_constraints_to_sf( sf, pose )

    # make and prepare the MonteCarlo object
    mc = MonteCarlo( pose, sf, temperature )

    # instantiate some data and counters
    num_small_moves = 0
    num_shear_moves = 0
    num_accepts = 0

    # make a list of max angles allowed for the loop residues, to be randomly selected from
    max_angles = [ 5.0, 6.0, 7.0 ]  # 6 is default

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
        
    # get the number of loops -- remember they are indexed starting at 1
    num_loops = loops.num_loop()
    cur_loop = 1
 
    # TODO-fix centroid mode for glycosylated proteins
    # switch the pose from full atom to centroid mode
    ##from rosetta import SwitchResidueTypeSetMover
    ##switch = SwitchResidueTypeSetMover( "centroid" )
    ##switch.apply( pose )
    
    # for each trial, perturbate each loop using  a MoveMap with Small or Shear moves
    # if <model_loops> is True, close the loop. Otherwise, just minimize it (as it wouldn't have been broken)
    if verbose:
        print "Making Small and Shear moves now..."
    for ii in range( 1, trials + 1 ):
        # pull out a loop to perturbate
        loop = loops[ cur_loop ]
        loop_start = loop.start()
        loop_stop = loop.stop()

        # update MoveMap to allow backbone and chi angle perturbation
        mm = MoveMap()
        mm.set_bb_true_range( loop_start, loop_stop )
        mm.set_chi_true_range( loop_start, loop_stop )
        
        
        # print out acceptance rate at intervals of 10
        if ii % 10 == 0:
            # calculate and print out the acceptance rate
            acceptance_rate = float( num_accepts ) / float( mc.total_trials() )
            print "  Acceptance rate:", str( round( acceptance_rate * 100, 3 ) ) + '%', "after", ii, "trials of", trials

        # either do a small or shear move, but not both
        if random.randint( 1, 2 ) == 1:  # 50% chance of one or the other
            num_small_moves += 1
            small_mover = SmallMover( mm, temperature, 5 )
           
            # use small mover and set a random angle max for loop
            small_mover.angle_max( 'L', random.choice( max_angles ) )  # L for loop
            small_mover.apply( pose )
            if pmm is not None and pmm_worked:
                pmm.apply( pose )
                
            # minimize Loop with the same MoveMap
            do_min_with_this_mm( mm, sf, pose )
            
            # close the loop and then minimize it if <model_loops> is True
            if model_loops:
                pose = CCD_loop_closure( loop, pose )
                do_min_with_this_mm( mm, sf, pose )
                if pmm is not None and pmm_worked:
                    pmm.apply( pose )

            # accept or reject move using MonteCarlo object
            if mc.boltzmann( pose ):
                num_accepts += 1
                    
        else:
            num_shear_moves += 1
            shear_mover = ShearMover( mm, temperature, 5 )
            
            # use shear mover and set a random angle max for loop
            shear_mover.angle_max( 'L', random.choice( max_angles ) )  # L for loop
            shear_mover.apply( pose )
            if pmm is not None and pmm_worked:
                pmm.apply( pose )
                
            # minimize Loop with the same MoveMap
            do_min_with_this_mm( mm, sf, pose )
            
            # close the loop and then minimize it if <model_loops> is True
            if model_loops:
                pose = CCD_loop_closure( loop, pose )
                do_min_with_this_mm( mm, sf, pose )
                if pmm is not None and pmm_worked:
                    pmm.apply( pose )
                    
            # accept or reject move using MonteCarlo object
            if mc.boltzmann( pose ):
                num_accepts += 1
                    
        # update the current loop index number, reseting to one if it's larger than the number of loops in the object
        if cur_loop == num_loops:
            cur_loop = 1
        else:
            cur_loop += 1
            
    # switch the pose back to full atom mode
    ##switch = SwitchResidueTypeSetMover( "fa_standard" )
    ##switch.apply( pose )
    
    # wrap up with pack and minimization
    print "Running a pack and minimization"
    pose = do_pack_min( sf, pose )
    
    if pmm is not None and pmm_worked:
        pmm.apply( pose )

    # if loops were modeled, restore the original FoldTree
    if model_loops:
        pose = restore_original_fold_tree( pose, verbose = verbose )
    
    # print out the percentage of Small moves made versus Shear moves during protocol
    print str( round( ( float( num_small_moves ) / float( trials ) * 100 ), 2 ) ) + '%', "small moves vs", str( round( ( float( num_shear_moves ) / float( trials ) * 100 ), 2 ) ) + '%', "shear moves given a total number of", trials, "trials"

    return pose

