#!/usr/bin/python

def get_E_per_residue_from_file( pdb_filename ):
    from rosetta import get_fa_scorefxn
    from antibody_functions import load_pose
    
    pose = load_pose( pdb_filename )
    sf = get_fa_scorefxn()
    sf( pose )
    
    for residue in pose:
        energy = pose.energies().residue_total_energy( residue.seqpos() )
        if energy > 1.5:
            print residue.name(), '\t', pose.pdb_info().pose2pdb( residue.seqpos() ), '\t', energy



def get_E_per_residue( pose ):
    from rosetta import get_fa_scorefxn
    
    sf = get_fa_scorefxn()
    sf( pose )
    
    for residue in pose:
        energy = pose.energies().residue_total_energy( residue.seqpos() )
        if energy > 4:
            print residue.name(), '\t', pose.pdb_info().pose2pdb( residue.seqpos() ), '\t', energy

if __name__ == "__main__":
    # parse and store args
    import argparse
    parser = argparse.ArgumentParser(description="Use PyRosetta to score a Pose per residue")
    parser.add_argument("pdb_filename", type=str, help="the filename of the PDB structure")
    input_args = parser.parse_args()
    
    get_E_per_residue_from_file( input_args.pdb_filename )

