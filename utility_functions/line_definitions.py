#!/usr/bin/python
__author__ = "morganlnance"


##############################
#### PDB LINE DEFINITIONS ####
##############################

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
    def __init__(self, line):
        # only works for LINK lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 173
        self.line = line.rstrip( '/n' )
        self.atom_name1 = line = str( self.line[ 12:16 ] )
        self.alt_loc1 = str( self.line[ 16:17] )
        self.res1_name = str( self.line[17:20] )
        self.res1_chain = str( self.line[21:22] )
        self.res1_seq = int( self.line[22:26] )
        self.i_code1 = str( self.line[26:27] )
        self.atom_name2 = str( self.line[42:46] )
        self.alt_loc2 = str( self.line[46:47] )
        self.res2_name = str( self.line[47:50] )
        self.res2_chain = str( self.line[51:52] )
        self.res2_seq = int( self.line[52:56] )
        self.i_code2 = str( self.line[56:57] )
        self.link_dist = float( self.line[73:78] )



class MODRES_line:
    def __init__(self, line):
        # only works for MODRES lines
        # line locations taken from the PDB databases's pdb format description page
        # ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        # page 157
        self.line = line.rstrip( '/n' )
        self.id_code = str( self.line[7:11].replace( ' ', '' ) )
        
        # the residue name used in the PDB file (post-modification)
        self.res_name = str( self.line[12:15].replace( ' ', '' ) )
        
        self.res_chain = str( self.line[16:17] )
        self.res_num = int( self.line[18:22].replace( ' ', '' ) )
        self.i_code = str( self.line[22:23] )

        # what the standard residue name is (pre-modification)
        self.std_res_name = str( self.line[24:27].replace( ' ', '' ) )
        
        self.comment = str( self.line[29:70] )



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
        
    def connection_types(self):
        return self.mmcif_dict[ "_struct_conn.conn_type_id" ]
    
    def partner1_chain_ids(self):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_asym_id" ] 
   
    def partner1_res_names(self):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_comp_id" ]
    
    def partner1_res_nums(self):
        return self.mmcif_dict[ "_struct_conn.ptnr1_auth_seq_id" ]
    
    def partner2_chain_ids(self):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_asym_id" ]
    
    def partner2_res_names(self):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_comp_id" ]
    
    def partner2_res_nums(self):
        return self.mmcif_dict[ "_struct_conn.ptnr2_auth_seq_id" ]
    
    def get_uniq_connection_names(self):
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
    def __init__(self, cif_filename):
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
        
    def atom_types(self):
        """
        Returns a list of "ATOM" and "HETATM" strings
        """
        return self.mmcif_dict[ "_atom_site.group_PDB" ]
    
    def chain_ids(self):
        return self.mmcif_dict[ "_atom_site.auth_asym_id" ]
    
    def res_names(self):
        return self.mmcif_dict[ "_atom_site.auth_comp_id" ]
    
    def res_nums(self):
        return self.mmcif_dict[ "_atom_site.auth_seq_id" ]
    
    def cif_uniq_pro_names(self):
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

    def cif_uniq_het_names(self):
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
