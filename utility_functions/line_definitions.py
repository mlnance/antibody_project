#!/usr/bin/python
__author__ = "morganlnance"


##############################
#### PDB LINE DEFINITIONS ####
##############################

class ATOM_line:
    def __init__( self, line ):
        '''
        # only works for ATOM lines
        line locations taken from the PDB databases's pdb format description page
        ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        page 187 for ATOM
        '''
        self.line = line.rstrip( "\n" )
        self.record_name = str( self.line[0:6].replace( ' ', '' ) )
        self.atom_num = int( self.line[6:11].replace( ' ', '' ) )
        self.atom_name = str( self.line[12:16].replace( ' ', '' ) )
        self.alt_loc = str( self.line[16:17].replace( ' ', '' ) )
        self.res_name = str( self.line[17:20].replace( ' ', '' ) )
        self.res_chain = str( self.line[21:22] )
        self.res_num = int( self.line[22:26].replace( ' ', '' ) )
        self.i_code = str( self.line[26:27] )
        self.x_coord = float( self.line[30:38].replace( ' ', '' ) )
        self.y_coord = float( self.line[38:46].replace( ' ', '' ) )
        self.z_coord = float( self.line[46:54].replace( ' ', '' ) )
        self.occupancy = float( self.line[54:60].replace( ' ', '' ) )
        self.temp_factor = float( self.line[60:66].replace( ' ', '' ) )
        self.element = str( self.line[76:78].replace( ' ', '' ) )
        self.charge = str( self.line[78:80].replace( ' ', '' ) )

        # store the string version of these specific floats as the trailing zero accuracy gets lost
        self.str_x_coord = str( self.line[30:38].replace( ' ', '' ) )
        self.str_y_coord = str( self.line[38:46].replace( ' ', '' ) )
        self.str_z_coord = str( self.line[46:54].replace( ' ', '' ) )
        self.str_occupancy = str( self.line[54:60].replace( ' ', '' ) )
        self.str_temp_factor = str( self.line[60:66].replace( ' ', '' ) )

        # self-defined information here
        # uniq_res_name = resname_reschain_resnum
        self.uniq_res_name = '_'.join( [ self.res_name, self.res_chain, str( self.res_num ) ] )


    def rebuild_ATOM_line( self ):
        rebuild_line = ''

        rebuild_line += ( self.record_name + ' ' * ( 6 - len( self.record_name ) ) )

        # atom numbers are right justified
        rebuild_line += "%s%s" %( ( ' ' * ( 5 - len( str( self.atom_num ) ) ) ), str( self.atom_num ) )

        # add the spacing between atom number and atom name
        rebuild_line += ' '

        # atom name is centered on the main atom letter (C, N, O, S, etc), so this is a case-by-case basis
        if len( self.atom_name ) == 1:
            # if the atom is one letter long, then it is ' X  '
            rebuild_line += " %s  " %self.atom_name
        elif len( self.atom_name ) == 2:
            # if the atom is two letters long, then it depends if there is a number on the left or a letter or number on the right
            # number case (typically hydrogens) '1X  '
            try:
                # check if the first character (the character on the left of the atom name) is a number
                int( self.atom_name[0] )
                rebuild_line += "%s  " %self.atom_name
            except:
                # or the letter case ' XB ' or the other number case ' X6 '
                rebuild_line += " %s " %self.atom_name
        elif len( self.atom_name ) == 3:
        # if the atom is three letters long, then there are two cases depending if there is a number to the left or not
        # number case '1XB '
            try:
                # check if the first character (the character on the left of the atom name) is a number
                int( self.atom_name[0] )
                rebuild_line += "%s " %self.atom_name
            except:
                # or the letters/numbers to the right only case ' XB1'
                rebuild_line += " %s" %self.atom_name
        else:
            # if the atom is four letters long, there is only one orientation '1XB2'
            rebuild_line += self.atom_name

        # if alternate location is blank, add a space instead
        if self.alt_loc == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.alt_loc

        # residue name is right justified
        rebuild_line += ( ' ' * ( 3 - len( self.res_name ) ) ) + self.res_name

        # add a space between residue name and the chain id
        rebuild_line += ' ' + self.res_chain

        # residue number is right justified
        rebuild_line += ( ' ' * ( 4 - len( str( self.res_num ) ) ) ) + str( self.res_num )

        # if insertion code is blank, add a space instead
        if self.i_code == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.i_code

        # add three spaces between insertion code and x coordinate
        rebuild_line += ' ' * 3

        # xyz coordinates are right justified and spaced out according to the length of the coordinates
        # the trailing 0('s) will be lost when doing float --> str, so use the stored string versions
        rebuild_line += ( ' ' * ( 8 - len( self.str_x_coord ) ) ) + self.str_x_coord
        rebuild_line += ( ' ' * ( 8 - len( self.str_y_coord ) ) ) + self.str_y_coord
        rebuild_line += ( ' ' * ( 8 - len( self.str_z_coord ) ) ) + self.str_z_coord

        # occupancy and temperature factors are right justified as well, with no intential spaces in between
        # the trailing 0('s) will be lost when doing float --> str, so use the stored string versions
        rebuild_line += ( ' ' * ( 6 - len( self.str_occupancy ) ) ) + self.str_occupancy
        rebuild_line += ( ' ' * ( 6 - len( self.str_temp_factor ) ) ) + self.str_temp_factor

        # add the spacing between the temperature factor and the element
        rebuild_line += ' ' * 10

        # add the element which is right justified
        rebuild_line += ( ' ' * ( 2 - len( self.element ) ) ) + self.element

        # the charge is two characters long minimum (1- or 2+, for example)
        rebuild_line += self.charge

        # add the rest of the trailing spaces
        rebuild_line += ' ' * ( 81 - len( rebuild_line ) )

        return rebuild_line


    def renumber_line( self, new_res_chain, new_res_num ):
        """
        Renumber a specific PDB line using a <new_res_num> and <new_res_chain> identification
        :param new_res_num: int( new residue number )
        """
        self.res_chain = new_res_chain
        self.res_num = new_res_num
        self.line = self.rebuild_ATOM_line()

        return self.line




class HETNAM_line:
    def __init__( self, line ):
        # only works for Rosetta sugar HETNAM lines
        self.line = line.rstrip( "\n" )
        self.record_name = str( self.line[0:6].replace( ' ', '' ) )
        self.res_name = str( self.line[11:14].replace( ' ', '' ) )
        self.res_chain = str( self.line[15:16] )
        self.res_num = int( self.line[16:20].replace( ' ', '' ) )
        self.i_code = str( self.line[20:21] )
        self.sugar_name = str( self.line[22:81].replace( ' ', '' ) )

        # self-defined information here
        # uniq_res_name = resname_reschain_resnum
        self.uniq_res_name = '_'.join( [ self.res_name, self.res_chain, str( self.res_num ) ] )


    def rebuild_Rosetta_HETNAM_line( self ):
        rebuild_line = ''

        rebuild_line += ( self.record_name + ' ' * ( 6 - len( self.record_name ) ) )

        # add the spaces between the record name and the residue name
        rebuild_line += ' ' * 5

        # residue name is right justified
        rebuild_line += ( ' ' * ( 3 - len( self.res_name ) ) ) + self.res_name

        # add a space between residue name and the chain id
        rebuild_line += ' ' + self.res_chain

        # residue number is right justified
        rebuild_line += ( ' ' * ( 4 - len( str( self.res_num ) ) ) ) + str( self.res_num ) 

        # add the insertion code or, if blank, a space
        if self.i_code == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.i_code

        # add a space between insertion code and the sugar name
        rebuild_line += ' ' + self.sugar_name

        # add the rest of the trailing spaces
        rebuild_line += ' ' * ( 81 - len( rebuild_line ) )

        return rebuild_line


    def renumber_line( self, new_res_chain, new_res_num ):
        """
        Renumber a specific PDB line using a <new_res_num> and <new_res_chain> identification
        :param new_res_num: int( new residue number )
        """
        self.res_chain = new_res_chain
        self.res_num = new_res_num
        self.line = self.rebuild_Rosetta_HETNAM_line()

        return self.line



class PDB_line:
    def __init__( self, line ):
        '''
        # only works for ATOM, HETATM, TER, and ANISOU lines
        line locations taken from the PDB databases's pdb format description page
        ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        page 187 for ATOM
        page 190 for ANISOU
        page 192 for TER
        page 194 for HETATM
        '''
        # found in ATOM, HETATM, TER, and ANISOU records
        self.line = line
        self.atom_num = int( self.line[6:11].replace( ' ', '' ) )
        self.atom_name = str( self.line[12:16].replace( ' ', '' ) )
        self.alt_loc = str( self.line[16:17].replace( ' ', '' ) )
        self.res_name = str( self.line[17:20].replace( ' ', '' ) )
        self.res_chain = str( self.line[21:22] )
        self.res_num = int( self.line[22:26].replace( ' ', '' ) )
        self.i_code = str( self.line[26:27] )

        # found in ATOM, HETATM, and ANISOU records
        if not line.startswith( "TER" ):
            self.element = str( self.line[76:78].replace( ' ', '' ) )
            self.charge = str( self.line[78:80].replace( ' ', '' ) )
            
        # found in ATOM and HETATM records
        if line.startswith( "ATOM" ) or line.startswith( "HETATM" ):
            self.x_coord = float( self.line[30:38].replace( ' ', '' ) )
            self.y_coord = float( self.line[38:46].replace( ' ', '' ) )
            self.z_coord = float( self.line[46:54].replace( ' ', '' ) )
            self.occupancy = float( self.line[54:60].replace( ' ', '' ) )
            self.temp_factor = float( self.line[60:66].replace( ' ', '' ) )
            self.element = str( self.line[76:78].replace( ' ', '' ) )
            self.charge = str( self.line[78:80].replace( ' ', '' ) )

        # found in ANISOU records
        if line.startswith( "ANISOU" ):
            self.u00 = int( self.line[28:35].replace( ' ', '' ) )
            self.u11 = int( self.line[35:42].replace( ' ', '' ) )
            self.u22 = int( self.line[42:49].replace( ' ', '' ) )
            self.u01 = int( self.line[49:56].replace( ' ', '' ) )
            self.u02 = int( self.line[56:63].replace( ' ', '' ) )
            self.u12 = int( self.line[63:70].replace( ' ', '' ) )



class LINK_line:
    def __init__( self, line ):
        # only works for LINK lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 173
        self.line = line.rstrip( "\n" )
        self.record_name = str( self.line[0:6].replace( ' ', '' ) )
        self.atom_name1 = str( self.line[12:16].replace( ' ', '' ) )
        self.alt_loc1 = str( self.line[16:17].replace( ' ', '' ) )
        self.res1_name = str( self.line[17:20].replace( ' ', '' ) )
        self.res1_chain = str( self.line[21:22] )
        self.res1_num = int( self.line[22:26].replace( ' ', '' ) )
        self.res1_i_code = str( self.line[26:27] )
        self.atom_name2 = str( self.line[42:46].replace( ' ', '' ) )
        self.alt_loc2 = str( self.line[46:47].replace( ' ', '' ) )
        self.res2_name = str( self.line[47:50].replace( ' ', '' ) )
        self.res2_chain = str( self.line[51:52] )
        self.res2_num = int( self.line[52:56].replace( ' ', '' ) )
        self.res2_i_code = str( self.line[56:57] )
        self.res1_sym = str( self.line[59:65].replace( ' ', '' ) )
        self.res2_sym = str( self.line[66:72].replace( ' ', '' ) )
        self.link_distance = float( self.line[73:78].replace( ' ', '' ) )

        # store the floats as strings as the trailing zeros get lost with float-->string conversion
        self.str_link_distance = str( self.line[73:78].replace( ' ', '' ) )

        # self-defined information here
        # uniq_res_name = resname1_reschain1_resnum1+resname2_reschain2_resnum2
        self.uniq_res1_name = '_'.join( [ self.res1_name, self.res1_chain, str( self.res1_num ) ] )
        self.uniq_res2_name = '_'.join( [ self.res2_name, self.res2_chain, str( self.res2_num ) ] )
        self.uniq_res_name = '+'.join( [ self.uniq_res1_name, self.uniq_res2_name ] )


    def rebuild_LINK_line( self ):
        rebuild_line = ''

        rebuild_line += ( self.record_name + ' ' * ( 6 - len( self.record_name ) ) )

        # add the spaces between the record name and the atom name
        rebuild_line += ' ' * 6

        # add the atom name
        # atom name is centered on the main atom letter (C, N, O, S, etc), so this is a case-by-case basis
        if len( self.atom_name1 ) == 1:
            # if the atom is one letter long, then it is ' X  '
            rebuild_line += " %s  " %self.atom_name1
        elif len( self.atom_name1 ) == 2:
            # if the atom is two letters long, then it depends if there is a number on the left or a letter or number on the right
            # number case (typically hydrogens) '1X  '
            try:
                # check if the first character (the character on the left of the atom name) is a number
                int( self.atom_name1[0] )
                rebuild_line += "%s  " %self.atom_name1
            except:
                # or the letter case ' XB ' or the other number case ' X6 '
                rebuild_line += " %s " %self.atom_name1
        elif len( self.atom_name1 ) == 3:
        # if the atom is three letters long, then there are two cases depending if there is a number to the left or not
        # number case '1XB '
            try:
                # check if the first character (the character on the left of the atom name) is a number
                int( self.atom_name1[0] )
                rebuild_line += "%s " %self.atom_name1
            except:
                # or the letters/numbers to the right only case ' XB1'
                rebuild_line += " %s" %self.atom_name1
        else:
            # if the atom is four letters long, there is only one orientation '1XB2'
            rebuild_line += self.atom_name1

        # add the alternate location or, if none, a space
        if self.alt_loc1 == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.alt_loc1

        # residue name is right justified
        rebuild_line += ( ' ' * ( 3 - len( self.res1_name ) ) ) + self.res1_name

        # add the chain id with the appropriate space
        rebuild_line += ' ' + self.res1_chain

        # residue number is right justified
        rebuild_line += ( ' ' * ( 4 - len( str( self.res1_num ) ) ) ) + str( self.res1_num ) 

        # add the insertion code or, if blank, a space
        if self.res1_i_code == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.res1_i_code

        # add the spaces in between the insertion code and the second atom name
        rebuild_line += ' ' * 15

        # add the second atom name
        # atom name is centered on the main atom letter (C, N, O, S, etc), so this is a case-by-case basis
        if len( self.atom_name2 ) == 1:
            # if the atom is one letter long, then it is ' X  '
            rebuild_line += " %s  " %self.atom_name2
        elif len( self.atom_name2 ) == 2:
            # if the atom is two letters long, then it depends if there is a number on the left or a letter or number on the right
            # number case (typically hydrogens) '1X  '
            try:
                # check if the first character (the character on the left of the atom name) is a number
                int( self.atom_name2[0] )
                rebuild_line += "%s  " %self.atom_name2
            except:
                # or the letter case ' XB ' or the other number case ' X6 '
                rebuild_line += " %s " %self.atom_name2
        elif len( self.atom_name2 ) == 3:
        # if the atom is three letters long, then there are two cases depending if there is a number to the left or not
        # number case '1XB '
            try:
                # check if the first character (the character on the left of the atom name) is a number
                int( self.atom_name2[0] )
                rebuild_line += "%s " %self.atom_name2
            except:
                # or the letters/numbers to the right only case ' XB1'
                rebuild_line += " %s" %self.atom_name2
        else:
            # if the atom is four letters long, there is only one orientation '1XB2'
            rebuild_line += self.atom_name2

        # add the alternate location or, if none, a space
        if self.alt_loc2 == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.alt_loc2

        # residue name is right justified
        rebuild_line += ( ' ' * ( 3 - len( self.res2_name ) ) ) + self.res2_name

        # add the chain id with the appropriate space
        rebuild_line += ' ' + self.res2_chain

        # residue number is right justified
        rebuild_line += ( ' ' * ( 4 - len( str( self.res2_num ) ) ) ) + str( self.res2_num ) 

        # add the insertion code or, if blank, a space
        if self.res2_i_code == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.res2_i_code

        # add the spaces between the second insertion code and and the symmetry operators
        rebuild_line += ' ' * 2

        # add the symmetry operators, or spaces if there are none
        # symmetry operator 1
        if self.res1_sym == '':
            rebuild_line += ' ' * 6
        else:
            rebuild_line += ( ' ' * ( 6 - len( self.res1_sym ) ) ) + self.res1_sym

        # add the space in between the two symmetry operators
        rebuild_line += ' '

        # symmetry operator 2
        if self.res2_sym == '':
            rebuild_line += ' ' * 6
        else:
            rebuild_line += ( ' ' * ( 6 - len( self.res2_sym ) ) ) + self.res2_sym

        # add the space between the last symmetry operator and the link distance
        rebuild_line += ' '

        # add the link distance
        rebuild_line += ( ' ' * ( 5 - len( self.str_link_distance ) ) ) + self.str_link_distance

        # add the rest of the trailing spaces
        rebuild_line += ' ' * ( 81 - len( rebuild_line ) )

        return rebuild_line


    def renumber_line( self, new_res1_chain, new_res1_num, new_res2_chain, new_res2_num ):
        """
        Renumber a specific PDB line using a <new_res_num> and <new_res_chain> identification
        :param new_res_num: int( new residue number )
        """
        self.res1_chain = new_res1_chain
        self.res1_num = new_res1_num
        self.res2_chain = new_res2_chain
        self.res2_num = new_res2_num
        self.line = self.rebuild_LINK_line()

        return self.line



class MODRES_line:
    def __init__( self, line ):
        # only works for MODRES lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 157
        self.line = line.rstrip( "\n" )
        self.id_code = str( self.line[7:11].replace( ' ', '' ) )
        
        # the residue name used in the PDB file (post-modification)
        self.res_name = str( self.line[12:15].replace( ' ', '' ) )
        
        self.res_chain = str( self.line[16:17] )
        self.res_num = int( self.line[18:22].replace( ' ', '' ) )
        self.i_code = str( self.line[22:23] )

        # what the standard residue name is (pre-modification)
        self.std_res_name = str( self.line[24:27].replace( ' ', '' ) )
        
        self.comment = str( self.line[29:70] )



class SSBOND_line:
    def __init__( self, line ):
        # only works for SSBOND line
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 171
        self.line = line.rstrip( "\n" )
        self.record_name = str( self.line[0:6].replace( ' ', '' ) )
        try:
            self.serial_number = int( self.line[7:10].replace( ' ', '' ) )
        except:
            self.serial_number = None
        self.res1_name = str( self.line[11:14].replace( ' ', '' ) )
        self.res1_chain = str( self.line[15:16] )
        self.res1_num = int( self.line[17:21].replace( ' ', '' ) )
        self.res1_i_code = str( self.line[21:22] )

        self.res2_name = str( self.line[25:28].replace( ' ', '' ) )
        self.res2_chain = str( self.line[29:30] )
        self.res2_num = int( self.line[31:35].replace( ' ', '' ) )
        self.res2_i_code = str( self.line[35:36] )

        self.res1_sym = str( self.line[59:65].replace( ' ', '' ) )
        self.res2_sym = str( self.line[66:72].replace( ' ', '' ) )
        self.bond_distance = float( self.line[73:78].replace( ' ', '' ) )

        # since trailing zeros get lost in float-->string, store the string version of some variables
        self.str_bond_distance = str( self.line[73:78].replace( ' ', '' ) )

        # self-defined information here
        # uniq_res_name = resname1_reschain1_resnum1+resname2_reschain2_resnum2
        self.uniq_res1_name = '_'.join( [ self.res1_name, self.res1_chain, str( self.res1_num ) ] )
        self.uniq_res2_name = '_'.join( [ self.res2_name, self.res2_chain, str( self.res2_num ) ] )
        self.uniq_res_name = '+'.join( [ self.uniq_res1_name, self.uniq_res2_name ] )


    def rebuild_SSBOND_line( self ):
        rebuild_line = ''

        rebuild_line += ( self.record_name + ' ' * ( 6 - len( self.record_name ) ) )

        # add the space between the record name and the serial number
        rebuild_line += ' '

        # add the serial number, if any. Otherwise, add the appropriate number of spaces
        if self.serial_number is not None:
            rebuild_line += self.serial_number
        else:
            rebuild_line += ' ' * 3

        # add the space between the serial number and the residue name
        rebuild_line += ' '

        # residue name is right justified
        rebuild_line += ( ' ' * ( 3 - len( self.res1_name ) ) ) + self.res1_name

        # add the chain id with the appropriate spaces
        rebuild_line += ' ' + self.res1_chain + ' '

        # residue number is right justified
        rebuild_line += ( ' ' * ( 4 - len( str( self.res1_num ) ) ) ) + str( self.res1_num ) 

        # add the insertion code or, if blank, a space
        if self.res1_i_code == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.res1_i_code

        # add the spaces in between the insertion code and the residue 2 name
        rebuild_line += ' ' * 3

        # residue name is right justified
        rebuild_line += ( ' ' * ( 3 - len( self.res2_name ) ) ) + self.res2_name

        # add the chain id with the appropriate spaces
        rebuild_line += ' ' + self.res2_chain + ' '

        # residue number is right justified
        rebuild_line += ( ' ' * ( 4 - len( str( self.res2_num ) ) ) ) + str( self.res2_num ) 

        # add the insertion code or, if blank, a space
        if self.res2_i_code == '':
            rebuild_line += ' '
        else:
            rebuild_line += self.res2_i_code

        # add the spaces between the second insertion code and the symmetry operator
        rebuild_line += ' ' * 24

        # add the symmetry operators, or spaces if there are none
        # symmetry operator 1
        if self.res1_sym == '':
            rebuild_line += ' ' * 6
        else:
            rebuild_line += ( ' ' * ( 6 - len( self.res1_sym ) ) ) + self.res1_sym

        # add the space in between the two symmetry operators
        rebuild_line += ' '

        # symmetry operator 2
        if self.res2_sym == '':
            rebuild_line += ' ' * 6
        else:
            rebuild_line += ( ' ' * ( 6 - len( self.res2_sym ) ) ) + self.res2_sym

        # add the bond distance
        rebuild_line += ( ' ' * ( 5 - len( self.str_bond_distance ) ) ) + self.str_bond_distance

        # add the rest of the trailing spaces
        rebuild_line += ' ' * ( 81 - len( rebuild_line ) )

        return rebuild_line


    def renumber_line( self, new_res1_chain, new_res1_num, new_res2_chain, new_res2_num ):
        """
        Renumber a specific PDB line using a <new_res_num> and <new_res_chain> identification
        :param new_res_num: int( new residue number )
        """
        self.res1_chain = new_res1_chain
        self.res1_num = new_res1_num
        self.res2_chain = new_res2_chain
        self.res2_num = new_res2_num
        self.line = self.rebuild_SSBOND_line()

        return self.line




##############################
#### CIF LINE DEFINITIONS ####
##############################

class cif_struct_conn_lines:
    def __init__( self, cif_filename ):
        """
        Give this class an object (list) of the lines in the cif file
        """
        # only works for the struct_conn line in an mmcif file (ie a xxxx.cif file)
        self.cif_filename = cif_filename
    
    def check_if_has_mmcif_dict( self ):
        try:
            self.mmcif_dict = MMCIF2Dict.MMCIF2Dict( self.cif_filename )
            self.mmcif_dict[ "_struct_conn.conn_type_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr1_auth_asym_id" ] 
            self.mmcif_dict[ "_struct_conn.ptnr1_auth_comp_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr1_auth_seq_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr2_auth_asym_id" ] 
            self.mmcif_dict[ "_struct_conn.ptnr2_auth_comp_id" ]
            self.mmcif_dict[ "_struct_conn.ptnr2_auth_seq_id" ]
            return True
        except:
            # print in color (if possible) letting the user know the cif file didn't work
            try:
                text = "~~Something wrong with the _struct_conn keys in " + self.cif_filename.split( '/' )[-1] + " using LINK records instead"
                print(Fore.YELLOW + text + Style.RESET_ALL)
            except:
                print "~~Something wrong with the _struct_conn keys in", self.cif_filename.split( '/' )[-1], "using LINK records instead"
            return False
        
    def connection_types( self ):
        return self.mmcif_dict[ "_struct_conn.conn_type_id" ]
    
    def partner1_chain_ids( self ):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_asym_id" ] 
   
    def partner1_res_names( self ):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_comp_id" ]
    
    def partner1_res_nums( self ):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_seq_id" ]
    
    def partner2_chain_ids( self ):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_asym_id" ]
    
    def partner2_res_names( self ):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_comp_id" ]
    
    def partner2_res_nums( self ):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_seq_id" ]
    
    def get_uniq_connection_names( self ):
        connection_types = self.connection_types()
        res1_names = self.partner1_res_names()
        res1_chains = self.partner1_chain_ids()
        res1_nums = self.partner1_res_nums()
        res2_names = self.partner2_res_names()
        res2_chains = self.partner2_chain_ids()
        res2_nums = self.partner2_res_nums()
        
        unique_connection_names = []
        for ii in range( len( connection_types ) ):
            if connection_types[ii] == "covale" or connection_types[ii] == "metalc":
                res1_unique_name = res1_names[ii] + '_' + res1_chains[ii] + '_' + res1_nums[ii]
                res2_unique_name = res2_names[ii] + '_' + res2_chains[ii] + '_' + res2_nums[ii]
                uniq_connection_name = res1_unique_name + '+' + res2_unique_name
                if uniq_connection_name not in unique_connection_names:
                    unique_connection_names.append( uniq_connection_name )
        return unique_connection_names        
    


class cif_atom_site_lines:
    def __init__( self, cif_filename):
        """
        Give this class an object (list) of the lines in the cif file
        """
        # only works for the atom_site line in an mmcif file (ie a xxx.cif file)
        try:
            self.cif_filename = cif_filename
            self.mmcif_dict = MMCIF2Dict.MMCIF2Dict( self.cif_filename )
            self.mmcif_dict[ "_atom_site.group_PDB" ]
            self.mmcif_dict[ "_atom_site.auth_asym_id" ]
            self.mmcif_dict[ "_atom_site.auth_comp_id" ]
            self.mmcif_dict[ "_atom_site.auth_seq_id" ]
        except:
            pass
        
    def atom_types( self ):
        """
        Returns a list of "ATOM" and "HETATM" strings
        """
        return self.mmcif_dict[ "_atom_site.group_PDB" ]
    
    def chain_ids( self ):
        return self.mmcif_dict[ "_atom_site.auth_asym_id" ]
    
    def res_names( self ):
        return self.mmcif_dict[ "_atom_site.auth_comp_id" ]
    
    def res_nums( self ):
        return self.mmcif_dict[ "_atom_site.auth_seq_id" ]
    
    def cif_uniq_pro_names( self ):
        atom_types = self.atom_types()
        res_names = self.res_names()
        res_chains = self.chain_ids()
        res_nums = self.res_nums()
        
        unique_protein_names = []
        for ii in range( len( atom_types ) ):
            if atom_types[ii] == "ATOM":
                unique_name = res_names[ii] + '_' + res_chains[ii] + '_' + res_nums[ii]
                if unique_name not in unique_protein_names:
                    unique_protein_names.append( unique_name )
        return unique_protein_names

    def cif_uniq_het_names( self ):
        atom_types = self.atom_types()
        res_names = self.res_names()
        res_chains = self.chain_ids()
        res_nums = self.res_nums()
        
        unique_hetatm_names = []
        for ii in range( len( atom_types ) ):
            if atom_types[ii] == "HETATM":
                if res_names[ii] != "HOH" and res_names[ii] != "DOD":
                    unique_name = res_names[ii] + '_' + res_chains[ii] + '_' + res_nums[ii]
                    if unique_name not in unique_hetatm_names:
                        unique_hetatm_names.append( unique_name )
        return unique_hetatm_names
