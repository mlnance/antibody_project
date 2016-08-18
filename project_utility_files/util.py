#!/usr/bin/python

from gzip import GzipFile
from StringIO import StringIO
import os
import urllib2
import string, random
from io import BytesIO
from os.path import isfile



def download_pdb(pdb_id, dest_dir):
    '''
    downloads pdbs, from yifans rosettacm py script
    :param pdb_id: 4letter accession code
    :param dest_dir:
    :return:
    '''
    try:
        pdb_name= pdb_id.split('_')[0]
    except:
        pass

    dest = '%s.pdb.gz' % ( pdb_name.lower() )
    filename = dest_dir + "/" + dest
    if not isfile(filename[:-3]):
        url = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' % ( pdb_name.upper() )
        #curl_cmd = 'curl -o %s/%s %s' % ( dest_dir, dest, url )
        wget_cmd = 'wget --quiet %s -O %s' % ( url, filename )
        #try:
        #    log_lines = os.popen( curl_cmd ).readlines()
        #except:
        log_lines = os.popen( wget_cmd ).readlines()

        filename = dest_dir + "/" + dest
        print filename
        os.system("gunzip -f %s" %filename)

        if ( os.path.isfile( filename[:-3] )):
            return dest
        else:
            print "Error: cannot download PDB %s!"%pdb_name

    return dest[:-3]



def calc_distance( vec1, vec2 ):
    """
    Calculates the distances between two points in 3D space
    :param vec1: list( x, y, z coordinates ) of point 1
    :param vec2: list( x, y, z coordinates ) of point 2
    :return: float( distance )
    """
    from math import sqrt, pow

    # assert that two lists were passed for coordinates
    if not isinstance( vec1, list ):
        print
        print "vec1 isn't a list - give me a list of xyz coordinates. Exiting"
        sys.exit()
    if not isinstance( vec2, list ):
        print
        print "vec2 isn't a list - give me a list of xyz coordinates. Exiting"
        sys.exit()

    # takes in two lists of xyz coordinates
    x1 = vec1[0]
    y1 = vec1[1]
    z1 = vec1[2]
    x2 = vec2[0]
    y2 = vec2[1]
    z2 = vec2[2]

    # calculate the distance
    dist = sqrt( pow( x2 - x1, 2 ) + pow( y2 - y1, 2 ) + pow( z2 - z1, 2 ) )

    return dist



def renumber_PDB_file( pdb_filename ):
    """
    Renumber each residue starting with residue 1 and chain A to residue X and chain Y
    :param pdb_filename: str( /path/to/pdb filename )
    :return: str( new pdb filename )
    """
    import sys
    from line_definitions import ATOM_line, HETNAM_line, \
        LINK_line, SSBOND_line


    # try to open up the PDB file
    try:
        f = open( pdb_filename, 'r' )
        pdb = f.readlines()
        f.close()
    except IOError:
        print pdb_filename, "doesn't exist in this directory. Did you mean to download it? Exiting."
        sys.exit()

    # create a dictionary of old unique names to new unique names (resname_reschain_resnum)
    old_to_new_dict = {}
    last_chain_seen = ''
    new_res_chain = 'A'
    new_res_num = 1

    for line in pdb:
        if line.startswith( "ATOM" ):
            ATOM = ATOM_line( line )

            # use the first line to create needed information
            if last_chain_seen == '':
                last_chain_seen = ATOM.res_chain

            # update the last_chain_seen
            if last_chain_seen != ATOM.res_chain:
                last_chain_seen = ATOM.res_chain
                new_res_chain = chr( ord( new_res_chain ) + 1 )

            # connect the new residue name with its old one
            if ATOM.uniq_res_name not in old_to_new_dict.keys():
                # create the new residue name
                new_res_name = '_'.join( [ ATOM.res_name, new_res_chain, str( new_res_num ) ] )

                # add the new name to the dictionary
                old_to_new_dict[ ATOM.uniq_res_name ] = new_res_name

                # increase the new_res_num
                new_res_num += 1

    # renumber each ATOM, HETNAM, SSBOND, and LINK line using the old_to_new_dict
    new_pdb = []
    for line in pdb:
        # otherwise collect relevant information depending on what kind of line this is
        if line.startswith( "ATOM" ):
            LINE = ATOM_line( line )
            new_uniq_res_name = old_to_new_dict[ LINE.uniq_res_name ]
            new_res_chain = new_uniq_res_name.split( '_' )[1]
            new_res_num = new_uniq_res_name.split( '_' )[2]
            new_line = LINE.renumber_line( new_res_chain, new_res_num )
            new_pdb.append( new_line + "\n" )

        elif line.startswith( "HETNAM" ):
            LINE = HETNAM_line( line )
            new_uniq_res_name = old_to_new_dict[ LINE.uniq_res_name ]
            new_res_chain = new_uniq_res_name.split( '_' )[1]
            new_res_num = new_uniq_res_name.split( '_' )[2]
            new_line = LINE.renumber_line( new_res_chain, new_res_num )
            new_pdb.append( new_line + "\n" )

        elif line.startswith( "LINK" ):
            LINE = LINK_line( line )

            # have to renumber two parts of this line
            new_uniq_res1_name = old_to_new_dict[ LINE.uniq_res1_name ]
            new_res1_chain = new_uniq_res1_name.split( '_' )[1]
            new_res1_num = new_uniq_res1_name.split( '_' )[2]
            new_uniq_res2_name = old_to_new_dict[ LINE.uniq_res2_name ]
            new_res2_chain = new_uniq_res2_name.split( '_' )[1]
            new_res2_num = new_uniq_res2_name.split( '_' )[2]

            new_line = LINE.renumber_line( new_res1_chain, new_res1_num, new_res2_chain, new_res2_num )
            new_pdb.append( new_line + "\n" )

        elif line.startswith( "SSBOND" ):
            LINE = SSBOND_line( line )

            # have to renumber two parts of this line
            new_uniq_res1_name = old_to_new_dict[ LINE.uniq_res1_name ]
            new_res1_chain = new_uniq_res1_name.split( '_' )[1]
            new_res1_num = new_uniq_res1_name.split( '_' )[2]
            new_uniq_res2_name = old_to_new_dict[ LINE.uniq_res2_name ]
            new_res2_chain = new_uniq_res2_name.split( '_' )[1]
            new_res2_num = new_uniq_res2_name.split( '_' )[2]

            new_line = LINE.renumber_line( new_res1_chain, new_res1_num, new_res2_chain, new_res2_num )
            new_pdb.append( new_line + "\n" )

        elif line.startswith( "TER" ):
            new_pdb.append( line )

        # else I don't care for the line
        else:
            pass

    return new_pdb



def renumber_PDB( pdb_filename, new_filename = None, reverse  = False ):
    # try to open up the PDB file
    try:
        f = open( pdb_filename, 'r' )
        pdb = f.readlines()
        f.close()

    except IOError:
        print pdb_filename, "doesn't exist in this directory. Did you mean to download it? Exiting."
        sys.exit()

    # make a smart whitespace checker  --  needed whitespace is as per PDB file formatting
    full_whitespace_needed = 4

    # start the residue counters
    res_num = 1
    num_of_uniq_residues = 1  # starts at one because indexed at 1

    # get the unique name of the first residue to begin the renumbering process
    for line in pdb:
        # make sure the line starts with ATOM or HETATM
        if line[0:4] == "ATOM" or line[0:6] == "HETATM":
            # get first unique residue identifier (resname + reschain + resnum)
            uniq_res_name = line[17:20].replace( ' ', '' ) + line[21:22].replace( ' ', '' ) + line[22:26].replace( ' ', '' )
            num_of_uniq_residues += 1
            break

    # count how many uniq residues there are if the user wants to number in reverse
    if reverse:
        for line in pdb:
            # make sure the line starts with ATOM or HETATM
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                # determine if this is the same residue as from before, or a new res
                cur_res_name = line[17:20].replace( ' ', '' ) + line[21:22].replace( ' ', '' ) + line[22:26].replace( ' ', '' )
                if cur_res_name != uniq_res_name:
                    uniq_res_name = cur_res_name
                    num_of_uniq_residues += 1

        # change the starting res_num to the number of unique residues
        res_num = num_of_uniq_residues

    # if a new filename was not given, make new_filename = pdb_filename
    if new_filename is None:
        new_filename = pdb_filename
        
    # now open the file to write into it
    with open( new_filename, 'wb' ) as pdb_fh:
        # renumber each residue, ignoring chain, from 1 - to end of protein
        for line in pdb:
            # make sure the line starts with ATOM or HETATM
            if line[0:4] == "ATOM" or line[0:6] == "HETATM":
                # determine if this is the same residue as from before, or a new res
                cur_res_name = line[17:20].replace( ' ', '' ) + line[21:22].replace( ' ', '' ) + line[22:26].replace( ' ', '' )
                if cur_res_name != uniq_res_name:
                    uniq_res_name = cur_res_name
                    if reverse:
                        # going to 1
                        res_num -= 1
                    else:
                        # going up to max
                        res_num += 1

                # get first unique residue identifier (resname + reschain + resnum)
                whitespace_needed = full_whitespace_needed - len( str( res_num ) )

                # update the line's residue number
                line = line[0:22] + ( ( ' ' * whitespace_needed ) + str( res_num ) ) + line[26:]

                # print line to file
                pdb_fh.write( line )

            else:
                pdb_fh.write( line )



def dump_pdb_by_chain( filename, pose, chains, decoy_num, dump_dir = None ):
    """
    Dump a .pdb file of <pose> including only the specified <chains>
    Sample input -- dump_pdb_by_chain( "test_out.pdb", my_pose, [ 'A', 'C', 'E' ]
    Sample input -- dump_pdb_by_chain( "test_out.pdb", my_pose, 'D' )
    :param filename: str( filename for the dumped PDB )
    :param pose: Pose()
    :param chains: list( chains of interest ) or str( chain )
    :param decoy_num: int( the number of the decoy for use when dumping the temporary PDB )
    :param dump_dir: str( /path/to/dump_dir for the temp pdb files made. Files will be deleted ). Default = None = current working directory
    :return: bool( True if dump successful, False if not )
    """
    #################
    #### IMPORTS ####
    #################

    import os, sys
    from line_definitions import PDB_line, HETNAM_line, \
        LINK_line, SSBOND_line

    # check validity of input args
    if not ( isinstance( chains, list ) or isinstance( chains, str ) ):
        print
        print "You didn't pass me a list or a string for the <chains> argument. Exiting."
        sys.exit()

    # make a random, 6-character suffix for the dump files because Jazz gets sassy?
    id = id_generator()

    # dump the given pose to get access to its PDB lines
    # if given a dump_dir
    if dump_dir is not None:
        if dump_dir.endswith( '/' ):
            dump_name = dump_dir + "%s_dump_%s.pdb" %( id, str( decoy_num ) )
        else:
            dump_name = dump_dir + "/%s_dump_%s.pdb" %( id, str( decoy_num ) )
    # else dump in the current working directory
    else:
        dump_name = os.getcwd() + "/%s_dump_%s.pdb" %( id, str( decoy_num ) )
    pose.dump_pdb( dump_name )

    # grab the dumped PDB lines
    with open( dump_name, "rb" ) as fh:
        pdb_lines = fh.readlines()

    # for each line, check the residues chain id and decide whether to keep it
    keep_these_lines = []
    for line in pdb_lines:
        # if this is a residue line
        if line.startswith( "ATOM" ):
            _PDB_line = PDB_line( line )

            # if you were given a list of chain id's to keep
            if isinstance( chains, list ):
                if _PDB_line.res_chain in chains:
                    keep_these_lines.append( line )

            # if you were given one chain to keep
            elif isinstance( chains, str ):
                if _PDB_line.res_chain == chains:
                    keep_these_lines.append( line )

        # otherwise it's some other kind of line, just keep it
        else:
            # if you were given a list of chain id's to keep
            if isinstance( chains, list ):
                # HETNAM lines
                if line.startswith( "HETNAM" ):
                    _HETNAM_line = HETNAM_line( line )
                    if _HETNAM_line.res_chain in chains:
                        keep_these_lines.append( line )

                # SSBOND lines
                if line.startswith( "SSBOND" ):
                    _SSBOND_line = SSBOND_line( line )
                    if _SSBOND_line.res1_chain in chains and _SSBOND_line.res2_chain in chains:
                        keep_these_lines.append( line )

                # LINK lines
                if line.startswith( "LINK" ):
                    _LINK_line = LINK_line( line )
                    if _LINK_line.res1_chain in chains and _LINK_line.res2_chain in chains:
                        keep_these_lines.append( line )

    # remove the dumped pdb
    cmd = "rm %s" %dump_name
    try:
        os.popen( cmd )
    except:
        pass

    # if no residues were selected, return False
    if len( keep_these_lines ) == 0:
        print "No residues matched your chain ID. I won't be dumping a file for you."
        return False

    # otherwise dump the selected residues into a PDB file
    try:
        with open( filename, "wb" ) as fh:
            fh.writelines( keep_these_lines )
        return True
    except:
        raise


    

def fancy_dump_pdb_by_chain( filename, pose, chains ):
    # TODO: figure this one out
    """
    Dump a .pdb file of <pose> including only the specified <chains>
    Sample input -- dump_pdb_by_chain( "test_out.pdb", my_pose, [ 'A', 'C', 'E' ]
    Sample input -- dump_pdb_by_chain( "test_out.pdb", my_pose, 'D' )
    :param filename: str( filename for the dumped PDB )
    :param pose: Pose()
    :param chains: list( chains of interest ) or str( chain )
    :return: bool( True if dump successful, False if not )
    """
    #################
    #### IMPORTS ####
    #################

    from rosetta.utility import vector1_Size, ostream



    # check validity of input args
    if not ( isinstance( chains, list ) or isinstance( chains, str ) ):
        print
        print "You didn't pass me a list or a string for the <chains> argument. Exiting."
        sys.exit()

    # for each residue, check its chain id and keep it in a vector1_Size list if desired
    keep_these_residues = vector1_Size()
    pdb_info = pose.pdb_info()
    for res in pose:
        chain_id = pdb_info.pose2pdb( res.seqpos() ).split( ' ' )[1]

        # if you were given a list of chains to keep
        if isinstance( chains, list ):
            if chain_id in chains:
                keep_these_residues.append( res.seqpos() )
        # if you were given one chain to keep
        elif isinstance( chains, str ):
            if chain_id == chains:
                keep_these_residues.append( res.seqpos() )

    # if no residues were selected, return False
    if len( keep_these_residues ) == 0:
        print "No residues matched your chain ID. I won't be dumping a file for you."
        return False

    # otherwise dump the selected residues into a PDB file
    fh = open( filename, "wb" )
    try:
        pose.dump_pdb( ostream( fh ), keep_these_residues )
        fh.close()        
        return True

    except:
        fh.close()    
        return False



def id_generator( size=6, chars=string.ascii_uppercase + string.digits ):
    '''
    Return an id of length <size> using the characters from <chars>
    :param size: int( desired length of id )
    :param chars: string( characters with which to choose from to create id )
    :return: string( unique id of size <size> )
    '''
    # if SystemRandom is available
    try:
        id = ''.join( random.SystemRandom().choice( chars ) for x in range( size ) )
    # otherwise, use standard random.choice
    except:
        id = ''.join( random.choice( chars ) for x in range( size ) )

    return id



class fasc_dict( dict ):
    """
    Allows you to create your own version of a dictionary so you can add attributes to it
    """
    pass



def read_fasc_file( fasc_filename ):
    """
    Read the <fasc_filename> and return a dictionary of the scoring data
    :param fasc_filename: str( /path/to/.fasc file
    :return: dict( .fasc scoring data )
    """
    #################
    #### IMPORTS ####
    #################

    import re
    
    try:
        import pandas
    except:
        import csv
    
    
    ##########################
    #### FASC FILE READER ####
    ##########################
     
    # create the dictionary that will hold the fasc_data_dict and its corresponding attributes
    fasc_data_dict = fasc_dict()
    fasc_data_dict.nstruct = None
    fasc_data_dict.pdb_name = None
    
    # open up the fasc_file
    with open( fasc_filename, "rb" ) as fh:
        for line in fh:
            # remove newline characters
            line = line.rstrip()
            
            # this should be the very first line of the .fasc file only
            if line.startswith( "pdb" ):
                # need to replace this specific space for clarity
                line = line.replace( "pdb name", "pdb_name" )
                
                # replace X amount of spaces that follow ':'
                arr = re.split( "[:\s]+", line )
            
    # turn the list of line data into an iterator
    iterator = iter( arr )
    
    # cycle through the first iterator to get the pdb_name and the nstruct
    # first iterator should always have this data
    for ii in iterator:
        key = str( ii )
        try:
            value = str( next( iterator ) )
            if key == "pdb_name":
                fasc_data_dict.pdb_name = value
            elif key == "nstruct":
                fasc_data_dict.nstruct = value
        except:
            pass
        
    # will hold valid decoy numbers in the fasc_data_dict object for easy access
    fasc_data_dict.decoy_nums = []
        
    # cycle through the .fasc file again to get the scoring terms and their values
    with open( fasc_filename, 'r' ) as fh:
        for line in fh:
            line = line.rstrip()

            # if the line is NOT the first line of the .fasc file
            if not line.startswith( "pdb" ):
                # split on multiple spaces after the ':' again
                arr = re.split( "[:\s]+", line )
                
                iterator = iter( arr )
                data_holder = {}
                for ii in iterator:
                    key = str( ii )
                    try:
                        # put the decoy data in the temporary holder
                        value = str( next( iterator ) )
                        data_holder[ key ] = value
                        
                        # pull out the filename decoy number
                        if key == "filename":
                            decoy_filename = value.split( '/' )[ -1 ]
                            if "sugar" in decoy_filename:
                                decoy_num = '_'.join( decoy_filename.split( '_' )[ -5: ] )
                            elif "LCM" or "GRM" in decoy_filename:
                                decoy_num = '_'.join( decoy_filename.split( '_' )[ -3: ] )
                            else: # SmallMoves
                                decoy_num = '_'.join( decoy_filename.split( '_' )[ -4: ] )
                            #decoy_num = int( value.replace( ".pdb", '' ).split( '_' )[-1] )
                            fasc_data_dict.decoy_nums.append( decoy_num )
                    except:
                        pass
                    
                # add the decoy data to the fasc_data_dict object given the decoy number
                fasc_data_dict[ decoy_num ] = data_holder
                
    # just because I'm picky - sort the decoy_nums
    fasc_data_dict.decoy_nums.sort()
    
    return fasc_data_dict



def get_score_term_from_fasc_data_dict( fasc_data_dict, score_term ):
    """
    After getting a <fasc_data_dict> from running read_fasc_file, pull out the <score_term> for each decoy in the .fasc file
    :param fasc_data_dict: dict( data dictionary from .fasc file )
    :param score_term: str( the score term you want back for each decoy
    :return: dict( decoy : score_term value )
    """
    # iterate through the dictionary skipping decoys that don't have that term
    score_term_dict = {}
    for decoy_num in fasc_data_dict.decoy_nums:
        try:
            score_term_dict[ decoy_num ] = float( fasc_data_dict[ decoy_num ][ score_term ] )
        # if it can only be a string
        except ValueError:
            score_term_dict[ decoy_num ] = fasc_data_dict[ decoy_num ][ score_term ]
        except:
            pass
        
    return score_term_dict



def get_mean( data ):
    """
    Return the mean of the integers and/or floats in <data>
    :param data: list( numbers )
    :return float( mean of data )
    """
    data = [ float( x ) for x in data ]
    return sum( data ) / float( len( data ) )
