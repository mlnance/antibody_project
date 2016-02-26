#/usr/bin/python


def make_all_mutations( pdb_filename, mutant_list_file, structure_dir ):
    import os
    import sys
    from antibody_functions import make_my_new_antibody, get_fa_scorefxn, load_pose, Pose
    
    try:
        mutant_list = []
        f = open( mutant_list_file, 'rb' )
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            if line != '':
                if line[0] != '#':
                    mutant_list.append( line )
    except IOError:
        print mutant_list_file, "didn't work. Check your input"
        sys.exit()
    except:
        print "Unexpected error"
        raise
    
        
    if not os.path.isdir( structure_dir ):
        os.mkdir( structure_dir )
    if not structure_dir.endswith( '/' ):
        structure_dir += '/'
    
    sf = get_fa_scorefxn()
    orig_pose = Pose()
    orig_pose.assign( load_pose( pdb_filename ) )
    
    for mut in mutant_list:
        mut_file_name = structure_dir + mut + ".pdb"
        if not os.path.isfile( mut_file_name ):
            print mut
            make_my_new_antibody( mut, sf, orig_pose, structure_dir )
        
        
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Use Rosetta to make point mutations.")
    parser.add_argument("pdb_filename", help="the filename of the PDB structure to serve as the native structure")
    parser.add_argument("mutant_list_file", help="give me the list of mutations you want to have made")
    parser.add_argument("structure_dir", help="where do you want your mutants to be dumped?")
    input_args = parser.parse_args()
    
    make_all_mutations( input_args.pdb_filename, input_args.mutant_list_file, input_args.structure_dir )
