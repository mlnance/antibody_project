#!/usr/bin/python
__author__ = "morganlnance"



# import needed functions
if __name__ == "__main__":
    print
    print "Importing modules..."
from rosetta import init, pose_from_file
import os, sys



# main function
def main( pose_obj_or_file, Angstroms_around, pickle_dir = None ):
    # load in the pose
    try:
        if __name__ == "__main__":
            init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records" )
            
            # load the pose from the filename
            pose = pose_from_file( pose_obj_or_file )
            
        # else it's an imported function, so load the pose
        else:
            # if the pose object passed is a str, import the pose
            if isinstance( pose_obj_or_file, str ):
                try:
                    pose = pose_from_file( pose_obj_or_file )
                except:
                    init( extra_options="-mute basic -mute core -mute protocols -include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records" )
                    pose = pose_from_file( pose_obj_or_file )

            else:
                pose = pose_obj_or_file
    except:
        print
        print "There was something wrong with the PDB file you gave me:", pose_obj_or_file, "Please check your input. Exiting"
        raise
    
    # instantiate a dictionary to store the data
    res_nums_around_res_given_X_cutoff = {}
    
    for res in pose:
        # get the center of the residue
        center = pose.residue( res.seqpos() ).nbr_atom_xyz()
        
        # holds the res nums that are around this particular residue number
        residues_around_this_res = []
        
        # for every other residue in the pose
        for surrounding_res in pose:
            # if it's not the current residue we're looking at
            if surrounding_res.seqpos() != res.seqpos():
                surrounding_center = pose.residue( surrounding_res.seqpos() ).nbr_atom_xyz()
                
                # check the distance between the center points
                if center.distance( surrounding_center ) <= Angstroms_around:
                    residues_around_this_res.append( surrounding_res.seqpos() )
                
        # add the residues to the dictionary according to residue number
        res_nums_around_res_given_X_cutoff[ res.seqpos() ] = residues_around_this_res
        
        # pull out the name of this PDB to make an appropriate pickle filename
        pdb_name = pose.pdb_info().name().split( '/' )[-1]


    # depending on how this program was run, return the dictionary or dump a pickle file
    if __name__ == "__main__":
        import pickle
        
        # if this was run as a standalone program, dump a pickle file
        # if no pickle_dir was passed, drop it in the current working directory
        if pickle_dir is None:
            pickle_dir = os.getcwd() + '/'
                
        # if a pickle_dir was passed but is not a valid direcotry, drop it in the current working dir
        if pickle_dir is not None and not os.path.isdir( pickle_dir ):
            pickle_dir = os.getcwd() + '/'
                
        # add a trailing '/' to the pickle_dir if needed
        if not pickle_dir.endswith( '/' ):
            pickle_dir += '/'
                
        # make the pickle file name given the pdb_name
        pickle_filename = ''.join( pickle_dir + pdb_name + "_residues_around_" + str( int( Angstroms_around ) ) + "_A.pickle" )
        print "Dumping:", pickle_filename
        
        # drop the dictionary as a pickle into the passed pickle directory
        pickle.dump( res_nums_around_res_given_X_cutoff, open( pickle_filename, "wb" ) )
    # else, this program was imported, so return the dictionary
    else:
        return res_nums_around_res_given_X_cutoff




if __name__ == "__main__":
    import argparse
    
    # parse and store input arguments
    parser = argparse.ArgumentParser(description="Use PyRosetta to glycosylate a pose and find a low E structure")
    parser.add_argument("pose_obj_or_file", type=str, help="/path/to/the PDB file to be used.")
    parser.add_argument("Angstroms_around", type=float, help="how many Angstroms around defines your neighest neighbors?")
    parser.add_argument("--pickle_dir", type=str, default=None, help="/path/to/dir where you would like a pickle file to be dumped.")
    input_args = parser.parse_args()


    main( input_args.pose_obj_or_file, input_args.Angstroms_around, input_args.pickle_dir )
