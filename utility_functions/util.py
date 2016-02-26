from gzip import GzipFile
from StringIO import StringIO
import os
import urllib2
import string
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
