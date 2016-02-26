import argparse

# parse and store arguments
parser = argparse.ArgumentParser(description="Use Rosetta to make and score point mutations.")
parser.add_argument("pdb_filename", help="the filename of the PDB structure to be evaluated")
parser.add_argument("-get_best_decoy", action="store_true", help="Do you want to run the mutational analysis multiple times to get the most accurate answer?")
args = parser.parse_args()

from mutational_analysis import *
from mother_nature_mutator import *

import random
import pandas as pd

class Master():
    def __init__(self):
        self.mutation = Mutation(args)
        self.mutation.load_pose()
        self.mn_mutator = MN_Mutator()
        

    def go(self):
        # load up the pose
        self.mutation.create()
        
        # create list of best pdb structures
        self.best_pdb_structures = []
        
        # initiate lists for dataframe data, and define the magic number of iterations
        self.mutation_names = []
        self.num_mutations = []
        self.mutation_score = []
        self.mutation_ddG_score = []
        self.mutation_ddG_interface = []
        self.mutation_hbond_ratio = []
        self.mutation_interface_contact_ratio = []
        self.final_comparison_score = []
        ii = 1
        magic_num = 2
        
        # get sequence positions of all residues in chain C for later atomic contact analysis
        # get sequence positions of all residues in chain A and B, because that's what I want
        # also ignore any residues that are sugars, or that are making non-polymeric connections
        A_B_positions = []
        C_positions = []
        for res_num in range(1, self.mutation.pose.total_residue() + 1):
            # C residues
            if self.mutation.pose.pdb_info().chain(res_num) == 'C':
                if not self.mutation.pose.residue(res_num).is_carbohydrate():
                    C_positions.append(res_num)
            # A and B residues
            if self.mutation.pose.pdb_info().chain(res_num) == 'A' or self.mutation.pose.pdb_info().chain(res_num) == 'B':
                if not self.mutation.pose.residue(res_num).is_carbohydrate():
                    if self.mutation.pose.residue(res_num).n_non_polymeric_residue_connections() == 0:
                        A_B_positions.append(res_num)
        
        # loop over all residues given a magic number, ignoring ones which should not be mutated
        print "Now using the mother nature mutation protocol"
        while ii < magic_num + 1:
            # pick random number of point mutations to make, 1 through 6
            num_point_mutations = random.randint(1, 6)
            print "Making", num_point_mutations, "point mutation(s)"
            
            orig_AA_list = []
            seq_pos_list = []
            for jj in range(num_point_mutations):
                # pick a random residue from chain A or B 
                random_res_pos = random.choice(A_B_positions)
                orig_AA_list.append(self.mutation.pose.residue(random_res_pos).name1())
                seq_pos_list.append(random_res_pos)
                
            # pick new AA to mutate to
            new_AA_list = []
            for AA in orig_AA_list:
                # generates a random AA to mutate to based on most likely codon mutation
                new_AA = self.mn_mutator.decide_my_mutation(AA)
                new_AA_list.append(new_AA)
            
            # make the name of the full mutation (make a list of the names and join the list with a '_')
            name_list = []
            for jj in range(num_point_mutations):
                name = orig_AA_list[jj] + str(seq_pos_list[jj]) + new_AA_list[jj]
                name_list.append(name)
            mutation_name = "_".join(name_list)
            
            # check to see if the mutation has already been done
            if mutation_name in self.mutation_names:
                # skip this mutation! it has already been done
                # set ii to what it currently is (ie. don't up the count)
                ii = ii
            
            else:
                # mutation time!
                # append the current information to the lists for the dataframe
                self.mutation_names.append(mutation_name)
                self.num_mutations.append(num_point_mutations)
                print "Mutation", ii, "is", mutation_name
                
                pack_min_rounds = 5 + 1 # 5 rounds, with 1 extra making the initial packed and minimized structure
                print "Packing and minimizing..."
                for index in range(num_point_mutations):
                    mut_AA = new_AA_list[index]
                    seq_pos = seq_pos_list[index]
                    
                    # if this is the first (or only) mutation, change the orig pose to return a mutant pose
                    if index == 0:
                        mutant_pose = Pose()
                        mutant_pose.assign(self.mutation.make_single_point_mutation_no_packmin(self.mutation.pose, seq_pos, mut_AA))
                        
                        # pack and minimize over specified interval, saving the best one
                        # start with a dummy best pose, and create a Pose object to hold the packed and minimized pose
                        best_mutant_pose = Pose()
                        best_mutant_score = self.mutation.sf(mutant_pose)
                        packmin_pose_holder = Pose()
                        for kk in range(pack_min_rounds):
                            # assign a temp pose so pack rounds are separate
                            temp_pose = Pose()
                            temp_pose.assign(mutant_pose)
                            
                            # pack and minimize while handling the mutated residue region
                            packmin_pose_holder.assign(self.mutation.do_mutation_pack_min(mut_AA, seq_pos, temp_pose))
                            
                            # score and compare, saving best score pose
                            check_score = self.mutation.sf(packmin_pose_holder)
                            if check_score < best_mutant_score:
                                best_mutant_pose.assign(packmin_pose_holder)
                                best_mutant_score = check_score
                    
                    # if this is anything but the first mutation, change the already-mutated pose again
                    else:
                        mutant_pose.assign(self.mutation.make_single_point_mutation_no_packmin(mutant_pose, seq_pos, mut_AA))
                        
                        # pack and minimize over specified interval, saving the best one
                        # start with a dummy best pose, and create a Pose object to hold the packed and minimized pose
                        best_mutant_pose = Pose()
                        best_mutant_score = self.mutation.sf(mutant_pose)
                        packmin_pose_holder = Pose()
                        for kk in range(pack_min_rounds):
                            # assign a temp pose so pack rounds are separate
                            temp_pose = Pose()
                            temp_pose.assign(mutant_pose)
                            
                            # pack and minimize while handling the mutated residue region
                            packmin_pose_holder.assign(self.mutation.do_mutation_pack_min(mut_AA, seq_pos, temp_pose))
                            
                            # score and compare, saving best score pose
                            check_score = self.mutation.sf(packmin_pose_holder)
                            if check_score < best_mutant_score:
                                best_mutant_pose.assign(packmin_pose_holder)
                                best_mutant_score = check_score
                                
                # add best_mutant to the best pose list
                self.best_pdb_structures.append(best_mutant_pose)
                
                # collect scoring data
                best_E = self.mutation.sf(best_mutant_pose)
                self.mutation_score.append(best_E)
                mut_ddG_score = best_E - self.mutation.nat_post_min_score
                self.mutation_ddG_score.append(mut_ddG_score)
                
                # get interface score
                temp_pose = Pose()
                temp_pose.assign(best_mutant_pose)
                jump = temp_pose.jump(2)
                vec = xyzVector_Real(400, 400, 400)
                jump.set_translation(vec)
                temp_pose.set_jump(2, jump)
                mut_split_E = self.mutation.sf(temp_pose)
                mut_interface_E = best_E - mut_split_E
                mut_ddG_interface = mut_interface_E - self.mutation.nat_interface_score
                self.mutation_ddG_interface.append(mut_ddG_interface)
                    
                # collect and store hydrogen bond ratio
                new_nhbonds = float( get_hbonds(best_mutant_pose).nhbonds() )
                hbond_ratio = self.mutation.orig_nhbonds/new_nhbonds
                self.mutation_hbond_ratio.append(hbond_ratio)

                # count, compare, and store interface contact ratio
                new_ninterface_contacts = self.mutation.count_interface_contacts(best_mutant_pose)
                interface_contact_ratio = self.mutation.orig_ninterface_contacts/new_ninterface_contacts
                self.mutation_interface_contact_ratio.append(interface_contact_ratio)
                
                # store final comparison score data
                # subtract 1 to hbond and interface contact ratio because 0 means no change (standard)
                # comparison score is ddG score and interface, hbond and interface ratio -1 and times -1
                # set a weight of 5 for hbonds and contacts to better catch any differences
                weight = 5
                final_score = 0
                final_score += mut_ddG_score
                final_score += mut_ddG_interface
                final_score += (weight * (-1 * (hbond_ratio - 1) ) )
                final_score += (weight * (-1 * (interface_contact_ratio - 1) ) )
                self.final_comparison_score.append(final_score)
                
                # up the counter and do it all again!
                ii += 1

        # collect all the data and store it in a pandas dataframe
        self.mutation_df = pd.DataFrame()
        self.mutation_df["Mut Name"] = self.mutation_names
        self.mutation_df["Num Muts"] = self.num_mutations
        self.mutation_df["Mut Score"] = self.mutation_score
        self.mutation_df["ddG Score"] = self.mutation_ddG_score
        self.mutation_df["ddG Interface"] = self.mutation_ddG_interface
        self.mutation_df["Hbond Ratio"] = self.mutation_hbond_ratio
        self.mutation_df["Interface Contact Ratio"] = self.mutation_interface_contact_ratio
        self.mutation_df["Final Weighted Score"] = self.final_comparison_score
        
        
        # store dataframe as a csv in current directory
        data_filename = "Mother_Nature_Mutation.csv"
        self.mutation_df.to_csv(data_filename)
        print 'Data written to:', data_filename
        sorted_df = self.mutation_df.sort("Final Weighted Score", ascending=False)
        print sorted_df
