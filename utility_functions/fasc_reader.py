#!/usr/bin/python

# This can be ran in the command line or in a program. This will take a specified fasc file and find the top [insert number here] files and return or print their path names

def read_fasc( fasc_filename, return_this_many_structures ):
    # try opening and reading the file
    try:
        with open( fasc_filename, 'rb' ) as fh:
            lines = fh.readlines()
    except IOError:
        print "There seems to be something wrong with the fasc filename you gave me.\n"
        raise
    except:
        print "I'm not sure what happened. Returning the error.\n"
        raise
            
            
    data_dict = {}
    data_list = []
    return_these_file_names = []
    
    # read and interpret the lines
    for line in lines:
        # first line in every fasc file - ignore it
        if not line.startswith( "pdb name" ):
            # if there wasn't a mistake and this isn't just a blank line with only a carriage return
            # or a single character with a carriage return
            if len( line ) != 1 and len( line ) != 2:
                # strip off the carriage return
                line = line.rstrip()
                
                try:
                    # get the pdb name
                    pdb_name = line.split( "filename: " )[1].split( ' ' )[0]
                    
                    # pull total score out
                    tot_score = float( line.split( "total_score: " )[1].split( ' ' )[0] )
                    
                    # add the data to a dictionary
                    data_dict[ pdb_name ] = tot_score
                    
                    # and add the score to a list
                    data_list.append( tot_score )
                except:
                    pass
                
                
    # sort the list
    data_list.sort()


    # get as many top scoring structures as the user requested
    top_scoring_values = data_list[ : return_this_many_structures ]
    
    # make the counter the number of top-scoring structures to return
    counter = return_this_many_structures
    
    # print out the names of the files that satisfy the user's request
    for key, value in data_dict.items():
        if counter != 0:
            # if the value in the data dictionary is the same as a value in the top scoring structures
            if value in top_scoring_values:
                
                # print the key -- the filename
                if __name__ == "__main__":
                    print key
                else:
                    return_these_file_names.append( key )
                    
                # decrease counter
                counter -= 1
            
    if __name__ != "__main__":
        return return_these_file_names
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Read a fasc file and return the names of the lowest energy PDBs")
    parser.add_argument("fasc_filename", type=str, help="the name of the fasc file")
    parser.add_argument("return_this_many_structures", type=int, help="how many names do you want this program to return? The top 3 lowest energy? Top 5? Give me a number")
    input_args = parser.parse_args()

    read_fasc( input_args.fasc_filename, input_args.return_this_many_structures )
