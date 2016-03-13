#!/usr/bin/python

import argparse

# TODO - for each structure, look through each decoy in directory and compare RMSD to score and pull out which one is better. When plotting on graph, have as many colors as there are number of decoys (so blue = lowest energy and to red = high energy, something like that

parser = argparse.ArgumentParser(description="Use Rosetta to score and compare mutated structures")
parser.add_argument("native_pdb_filename", help="the filename of the PDB structure to serve as the base structure")
parser.add_argument("interface_jump_num", type=int, help="what jump number defines the interface location? Used for calculating binding energy differences.")
parser.add_argument("structure_dir", type=str, help="where are the structures you wish to analyze against the native? Give me the directory")
parser.add_argument("resulting_filename", type=str, help="what do you want the resulting csv file to be called? This program will add the .csv extension for you")
parser.add_argument("--verbose", "-v", default=False, action="store_true", help="do you want the program to tell you what it's doing?")
input_args = parser.parse_args()

from antibody_functions import *
import os

'''
# define the residue numbers that comprise the CH2 and CH3 domains, loop regions, and the sugars ( pose numbering )
# only relevant to 3ay4 FcgRIIIa
CH2_A_domain = range( 1, 108 + 1 )  # from CYS at linker region to right before loop A
loop_A = range(109, 118 + 1 ) # connects CH2 to CH3 domain in chain A
CH3_A_domain = range( 119, 215 + 1 )  # from end of loop A to the rest of chain A
sugar_A = range( 216, 223 + 1 )  # the entire glycan on ASN 297 of chain A
CH2_B_domain = range( 224, 331 + 1 )  # from CYS at linker region to right before loop B
loop_B = range( 332, 341 + 1 ) # connects CH2 to CH3 domain in chain B
CH3_B_domain = range( 342, 439 + 1 )  # from end of loop B to the rest of chain B
sugar_B = range( 440, 447 + 1 )  # the entire glycan on ASN 297 of chain B
receptor_protein = range( 448, 607 + 1 )  # all amino acids within the receptor
receptor_sugar = range( 608, 615 + 1 )  # the glycan that docks with the Fc region
'''

# check validity of structure directory
working_dir = os.getcwd()
if os.path.isdir( input_args.structure_dir ):
    if input_args.structure_dir.endswith( '/' ):
        structure_dir = input_args.structure_dir
    else: 
        structure_dir = input_args.structure_dir + '/'           
else:
    print "It seems like you gave me an incorrect path, exiting"
    sys.exit()

# get all the structure names from the structure directory 
os.chdir( structure_dir )
structures = []
structure_names = []
for f in os.listdir( os.getcwd() ):
    if f.endswith( ".pdb" ):
        structures.append( os.path.abspath( f ) )
        # I don't want the full path - just the name of the file ie. the structure name
        structure_names.append( f.split( '/' )[-1] )        
os.chdir( working_dir )

# inform the user of the structure directory and files to be analyzed
print "From", structure_dir, "analyzing:"
print structure_names

        
class Analyze():
    def __init__(self):
        # for use within go()
        self.sf = get_fa_scorefxn()
        
        
    def go(self):
        # set up and store the input native structure as the base pdb
        print "Wild type pose:", input_args.native_pdb_filename
        self.WT_pose = load_pose( input_args.native_pdb_filename )
        
        # calculate all base scores for comparison
        self.WT_score = self.sf( self.WT_pose )
        self.WT_elec_score = get_score_by_scoretype( self.sf, "fa_elec", self.WT_pose )
        self.WT_interface_score = get_interface_score( input_args.interface_jump_num, self.sf, self.WT_pose )
        self.WT_hbonds = count_hbonds( self.WT_pose )
        self.WT_interface_contacts, self.WT_interface_contacts_list = count_interface_atomic_contacts( input_args.interface_jump_num, self.WT_pose, verbose = input_args.verbose )
        self.WT_interface_SASA = analyze_interface( input_args.interface_jump_num, self.WT_pose, pack_separated = True )
        
        # instantiate lists to save data for data frame
        mutation_made = []
        pose_score = []
        dG_score = []
        elec_score = []
        dG_elec_score = []
        interface_score = []
        ddG_interface = []
        hbonds = []
        dHbonds = []
        interface_contacts = []
        dinterface_contacts = []
        new_interface_contacts = []
        tot_dinterface_contacts = []
        interface_SASA = []
        dInterface_SASA = []
        rmsd = []
        
        # analyze dat data!
        for pdb in structures:
            # load mutated pdb
            pose = load_pose( pdb )
            if input_args.verbose:
                print "Working on", pdb.split( '/' )[-1]
            
            # collect data
            mut_score = self.sf( pose )
            mut_elec_score = get_score_by_scoretype( self.sf, "fa_elec", pose )
            mut_interface_score = get_interface_score( input_args.interface_jump_num, self.sf, pose )
            mut_hbonds = count_hbonds( pose )
            mut_interface_contacts, mut_interface_contacts_list = count_interface_atomic_contacts( input_args.interface_jump_num, pose, verbose = input_args.verbose )
            mut_interface_SASA = analyze_interface( input_args.interface_jump_num, pose, pack_separated = True )
            mut_rmsd = CA_rmsd( self.WT_pose, pose )
            rmsd.append( mut_rmsd )
            
            # add to lists - unless it was the native, I don't want that there
            if not pose.pdb_info().name() == self.WT_pose.pdb_info().name():
                mutation_made.append( pdb.split( '/' )[-1] )
                pose_score.append( mut_score )
                dG_score.append( mut_score - self.WT_score )
                elec_score.append( mut_elec_score )
                dG_elec_score.append( mut_elec_score - self.WT_elec_score )
                interface_score.append( mut_interface_score )
                ddG_interface.append( mut_interface_score - self.WT_interface_score )
                hbonds.append( mut_hbonds )
                dHbonds.append( mut_hbonds - self.WT_hbonds )
                interface_SASA.append( mut_interface_SASA )
                dInterface_SASA.append( mut_interface_SASA - self.WT_interface_SASA )
                
                native_contacts = 0
                new_contacts = 0
                for contact in mut_interface_contacts_list:
                    if contact in self.WT_interface_contacts_list:
                        native_contacts += 1
                    else:
                        new_contacts += 1
                interface_contacts.append( mut_interface_contacts )
                dinterface_contacts.append( native_contacts - self.WT_interface_contacts )
                new_interface_contacts.append( new_contacts )
                tot_dinterface_contacts.append( ( native_contacts - self.WT_interface_contacts ) + new_contacts )
                
        # return to working dir
        os.chdir( working_dir )
        
        # check out the passed result filename for .csv extension
        filename = input_args.resulting_filename
        if not filename.endswith( ".csv" ):
            filename = filename + ".csv"
        
        # make and dump pandas dataframe as a csv
        self.df = pd.DataFrame( index=mutation_made )
        
        self.df["Pose Score"] = pose_score
        self.df["dG Scores"] = dG_score
        self.df["Elec Score"] = elec_score
        self.df["dG Elec Scores"] = dG_elec_score
        self.df["Interface Score"] = interface_score
        self.df["dG Interface"] = ddG_interface
        self.df["Num Hbonds"] = hbonds
        self.df["delta Hbonds"] = dHbonds
        self.df["interface contacts"] = interface_contacts
        self.df["lost native interface contacts"] = dinterface_contacts
        self.df["new interface contacts"] = new_interface_contacts
        self.df["total change interface contacts"] = tot_dinterface_contacts
        self.df["Interface SASA"] = interface_SASA
        self.df["delta Interface SASA"] = dInterface_SASA
        self.df["rmsd"] = rmsd
        
        print "Dumping", filename
        self.df.to_csv( filename, index=True, index_label = "Mutation Made" )


# runs the program
my_obj = Analyze()
my_obj.go()
