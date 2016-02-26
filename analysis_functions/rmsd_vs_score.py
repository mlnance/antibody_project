#!/usr/bin/python

import argparse

parser = argparse.ArgumentParser(description="Use Rosetta to calculate RMSD between a native pose and a directory of structures")
parser.add_argument("native_pdb_filename", help="the filename of the PDB structure to serve as the native structure")
parser.add_argument("structure_dir", help="where do the structures to which I am comparing the native live?")
input_args = parser.parse_args()

from mutational_analysis import *

# check the structure directory
working_dir = os.getcwd() + '/'
try:
    if os.path.isdir( input_args.structure_dir ):
        structure_dir = input_args.structure_dir
        if not structure_dir[-1] == '/':
            structure_dir = structure_dir + '/'
        files = os.listdir( structure_dir )
        pdbs = []
        for f in files:
            if f.endswith( ".pdb" ):
                f = structure_dir + f
                pdbs.append( os.path.abspath( f ) )
except:
    print "It appears", input_args.structure_dir, "is not a valid directory path. Exiting"
    sys.exit()

# check and load native pose
try:
    if os.path.isfile( input_args.native_pdb_filename ):
        native = load_pose( input_args.native_pdb_filename )
except:
    print "It appears", input_args.native_pdb_filename, "is not a valid pdb file. Exiting"
    sys.exit()
    

# make a scorefunction
sf = get_fa_scorefxn()
sf = apply_sugar_constraints_to_sf( sf, native )

# collect the data
rmsds = []
scores = []
pdb_names = []
for pdb in pdbs:
    mutant = load_pose( pdb )
    
    name = pdb.split( '/' )[-1]
    pdb_names.append( name )
    rmsd = CA_rmsd( native, mutant )
    rmsds.append( rmsd )
    score = sf( mutant )
    scores.append( score )


# dump the data
df = pd.DataFrame( pdb_names )
df["scores"] = scores
df["rmsd"] = rmsds
print df
df.to_csv( "Protocol_on_Lowest_E_Pack_Min_3ay4_again_RMSD_vs_Score.csv", index=False )
