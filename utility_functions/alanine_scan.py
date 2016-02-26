# import rosetta-specific modules
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
    standard_packer_task, change_cys_state, \
    Pose, MoveMap, RotamerTrialsMover, MinMover, \
    PyMOL_Mover
from toolbox import mutate_residue

# imports for data analysis
import os.path
import sys
import pandas as pd

# create global pymol object for viewing
pmm = PyMOL_Mover()

# path for mutational analysis data
data_dir = "/Users/Research/pyrosetta_repo/mutational_data/"
data_filename = "alanine_scan.csv"

# initialize rosetta with sugar flags
init(extra_options="-include_sugars -override_rsd_type_limit -read_pdb_link_records")

# create pose object, assign info object, and send pose to pymol
pose = pose_from_pdb("Fc_FcgRIII.pdb")
info = pose.pdb_info()
#pmm.apply(pose)

# instantiate full atom score function and score pose
sf = get_fa_scorefxn()
start_score = sf(pose)

# create packer task
task = standard_packer_task(pose)
task.restrict_to_repacking()
task.or_include_current(True)

# run a RotamerTrialsMover
pack_rotamers_mover = RotamerTrialsMover(sf, task)
pack_rotamers_mover.apply(pose)
#pmm.apply(pose)

# insantiate a movemap
mm = MoveMap() 
mm.set_bb(True)
mm.set_chi(True)
mm.set_nu(True)

# turn off chi and nu for sugar residues
## this has been hardcoded!
sugar_res_pos = []
for x in range(216,223 + 1) + range(440,447 + 1) + range(608,618 + 1):
    sugar_res_pos.append(x)
for x in sugar_res_pos:
    mm.set_chi(x, False)
    mm.set_nu(x, False)
    
# minimize using MinMover
min_mover = MinMover(mm, sf, "dfpmin", 0.01, True)
min_mover.apply(pose)

# grab post-minimization score and send to pymol
post_min_score = sf(pose)
#pmm.apply(pose)



# check to see if alanine scan has already been done
# if yes, print data
if(os.path.exists(data_dir + data_filename)):
    print "\nAlanine scan has already been ran! Check mutational_data for csv file."
    df = pd.DataFrame.from_csv(data_dir + data_filename)
    print df
    sys.exit()
    
# if no, get data
else:
    # create some constant data
    AA = 'A'
    PACK_RADIUS = 10.0
    
    # insantiate lists for data (will be put into pandas dataframe)
    orig_AA = []
    position = []
    native_E = []
    mut_E = []
    ddG = []
    
    # for all residues, mutate to alanine
    for seq_pos in range(1, pose.total_residue() + 1):
        res = pose.residue(seq_pos)
        
        # make sure residue is not a disulfide, sugar, or a branch point
        if res.name() != "CYD":
            if not res.is_carbohydrate():
                if res.name1() != 'A':
                    if not res.is_branch_point():
                        # append original residue, position, and score to list
                        orig_AA.append(pose.residue(seq_pos).name1())
                        position.append(seq_pos)
                        native_E.append(sf(pose))
                        
                        # mutate!
                        print "Mutating position", seq_pos, "..."
                        mutant = mutate_residue(pose, seq_pos, AA, PACK_RADIUS, sf)
                    	#pmm.apply(mutant)
                        
                        # minimize mutation
                        min_mover.apply(mutant)
                    	#pmm.apply(mutant)
                        
                        # score and add to list for dataframe
                        new_E = sf(mutant)
                        mut_E.append(new_E)
                        ddG.append(new_E - sf(pose))
                    
                
# print data and add data to dataframe
df = pd.DataFrame()
df["orig_AA"] = orig_AA
df["seq_pos"] = position
df["native_E"] = native_E
df["mut_E"] = mut_E
df["ddG"] = ddG
print df

# output results to file
with open(data_dir + data_filename, 'a') as f:
    df.to_csv(data_dir + data_filename)
print 'Data written to:', data_dir + data_filename
