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

    # dump the given pose to get access to its PDB lines
    # make a random, 4-character suffix for the dump files because Jazz gets sassy?

    # if given a dump_dir
    if dump_dir is not None:
        if dump_dir.endswith( '/' ):
            dump_name = dump_dir + "dump_%s.pdb" %str( decoy_num )
        else:
            dump_name = dump_dir + "/dump_%s.pdb" %str( decoy_num )
    # else dump in the current working directory
    else:
        dump_name = os.getcwd() + "/dump_%s.pdb" %str( decoy_num )
    pose.dump_pdb( dump_name )

    # grab the dumped PDB lines
    with open( dump_name, "rb" ) as fh:
        pdb_lines = fh.readlines()

    # remove the dumped pdb
    cmd = "rm %s" %dump_name
    try:
        os.popen( cmd )
    except:
        pass

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



def id_generator( size=4, chars=string.ascii_uppercase + string.digits ):
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
