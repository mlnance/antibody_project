#!/usr/bin/python

# This takes a fasc file and reads it and then moves x number of the files into the given directory. Also makes a file containing the path to the lowest energy file as well as the x number of low energy files

def main( fasc_file, low_e_dir, return_this_many_structures ):
    import os
    import sys
    from fasc_reader import read_fasc
    
    # get current directory
    cur_dir = os.getcwd()
    
    # check the validity of the directory and file
    if not os.path.isdir( low_e_dir ):
        print "Something is wrong with your low_e directory. Check your input. Exiting."
        sys.exit()
    else:
        low_e_dir = low_e_dir
    if not os.path.isfile( fasc_file ):
        print "Something is wrong with your fasc file. Exiting."
        sys.exit()
    else:
        fasc_file = fasc_file
        
    # get the absolute paths to everything
    os.chdir( low_e_dir )
    low_e_dir = os.path.abspath( os.getcwd() ) + '/'
    os.chdir( cur_dir )
    
    # use fasc_reader.py to get the relevant file paths
    fasc_files = read_fasc( fasc_file, return_this_many_structures )
    
    # get pdb filename by doing some fancy spliting (ex. /path/to/native_structure_100.pdb --> native_structure)
    split_on_pdb = fasc_files[0].split( ".pdb" )[0]
    split_on_underscore = split_on_pdb.split( '_' )
    pdb_name_with_path = '_'.join( split_on_underscore[ : len( split_on_underscore ) - 1 ] )
    pdb_name = pdb_name_with_path.split( '/' )[-1]
    
    # make a directory for the pdb in the low_e_dir if needed called dump_dir
    dump_dir = low_e_dir + pdb_name + '/'
    if not os.path.isdir( dump_dir ):
        os.mkdir( dump_dir )
        
    # write the top x structures into a file in the dump_dir
    lowest_x_file_name = dump_dir + "lowest_" + str( return_this_many_structures ) + "_total_score"
    with open( lowest_x_file_name, 'wb' ) as fh:
        for f_name in fasc_files:
            line = f_name + "\n"
            fh.write( line )
            
    # get and write the lowest energy structure into a file in the dump_dir
    lowest_e_file = read_fasc( fasc_file, 1 )[0]
    lowest_e_file_name = dump_dir + "lowest_e_total_score"
    with open( lowest_e_file_name, 'wb' ) as fh:
        fh.write( lowest_e_file + "\n" )
        
    # copy the lowest energy structure file into all_lowest_E_structs dir in the low_e_dir
    all_lowest_E_structs = low_e_dir + "all_lowest_E_structs" + '/'
    if not os.path.isdir( all_lowest_E_structs ):
        os.mkdir( all_lowest_E_structs )
    os.chdir( all_lowest_E_structs )
    cmd = "cp %s %s" %( lowest_e_file, all_lowest_E_structs )
    os.popen( cmd )
    os.chdir( cur_dir )
    
    # for each file name given, move it to the specified directory
    os.chdir( dump_dir )
    for f in fasc_files:
        # copy over the files
        cmd = "cp %s %s" %( f, dump_dir )
        os.popen( cmd )
    os.chdir( cur_dir )
    
    # if this is an imported function, return the path to the lowest E pose
    if __name__ != "__main__":
        return lowest_e_file


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Use Python to move x number of low-energy structures into a directory based on its fasc file")
    parser.add_argument("fasc_file", type=str, help="the path to the fasc file for the structure of interest")
    parser.add_argument("low_e_dir", type=str, help="the path to the directory where the structures should be copied to")
    parser.add_argument("return_this_many_structures", type=int, help="how many names do you want this program to return? The top 3 lowest energy? Top 5? Give me a number")
    input_args = parser.parse_args()
    
    main( input_args.fasc_file, input_args.low_e_dir, input_args.return_this_many_structures )
