import random

AA_one_letter_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
AA_three_letter_list = ["Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Try"]
nucleotide_list = ['A', 'C', 'G', 'U']

# list of codons for each AA 
Ala = ["GCU", "GCC", "GCA", "GCG"]
Cys = ["UGU", "UGC"]
Asp = ["GAU", "GAC"]
Glu = ["GAA", "GAG"]
Phe = ["UUU", "UUC"]
Gly = ["GGU", "GGC", "GGA", "GGG"]
His = ["CAU", "CAC"]
Ile = ["AUU", "AUC", "AUA"]
Lys = ["AAA", "AAG"]
Leu = ["CUU", "CUC", "CUA", "CUG", "UUA", "UUG"]
Met = ["AUG"]
Asn = ["AAU", "AAC"]
Pro = ["CCU", "CCC", "CCA", "CCG"]
Gln = ["CAA", "CAG"]
Arg = ["CGU", "CGC", "CGA", "CGG"]
Ser = ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"]
Thr = ["ACU", "ACC", "ACA", "ACG"]
Val = ["GUU", "GUC", "GUA", "GUG"]
Trp = ["UGG"]
Tyr = ["UAU", "UAC"]

codon_dict = {"UUU" : 'F', "UUC" : 'F', "UUA" : 'L', "UUG" : 'L', "CUU" : 'L', "CUC" : 'L', "CUA" : 'L', "CUG" : 'L', "AUU" : 'I', "AUC" : 'I', "AUA" : 'I', "AUG" : 'M', "GUU" : 'V', "GUC" : 'V', "GUA" : 'V', "GUG" : 'V', "UCU" : 'S', "UCC" : 'S', "UCA" : 'S', "UCG" : 'S', "CCU" : 'P', "CCC" : 'P', "CCA" : 'P', "CCG" : 'P', "ACU" : 'T', "ACC" : 'T', "ACA" : 'T', "ACG" : 'T', "GCU" : 'A', "GCC" : 'A', "GCA" : 'A', "GCG" : 'A', "UAU" : 'Y', "UAC" : 'Y', "UAA" : '*', "UAG" : '*', "CAU" : 'H', "CAC" : 'H', "CAA" : 'Q', "CAG" : 'Q', "AAU" : 'N', "AAC" : 'N', "AAA" : 'K', "AAG" : 'K', "GAU" : 'D', "GAC" : 'D', "GAA" : 'E', "GAG" : 'E', "UGU" : 'C', "UGC" : 'C', "UGA" : '*', "UGG" : 'W', "CGU" : 'R', "CGC" : 'R', "CGA" : 'R', "CGG" : 'R', "AGU" : 'S', "AGC" : 'S', "AGA" : 'R', "AGG" : 'R', "GGU" : 'G', "GGC" : 'G', "GGA" : 'G', "GGG" : 'G'}

# links the codons to the single-letter amino acid codes
AA_dict = {'A' : Ala, 'C' : Cys, 'D' : Asp, 'E' : Glu, 'F' : Phe, 'G' : Gly, 'H' : His, 'I' : Ile, 'K' : Lys, 'L' : Leu, 'M' : Met, 'N' : Asn, 'P' : Pro, 'Q' : Gln, 'R' : Arg, 'S' : Ser, 'T' : Thr, 'V' : Val, 'W' : Trp, 'Y' : Tyr}


def initialize():
    AA_seq = "TCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK"
    codon_seq_len = len(AA_seq) * 3
    all_dna_seq = []
    magic_num = 1000
    
    for i in range(magic_num):
        rando_dna = ""
        for ii in AA_seq:
            rando_dna = rando_dna + AA_dict[ii][random.randrange( 0, len(AA_dict[ii]) )]
        all_dna_seq.append(rando_dna)
        
    rand_choice = random.randrange(0, magic_num)
    orig_dna_seq = all_dna_seq[rand_choice]
        
    #print "Your amino acid sequence is:", AA_seq
    #print "Your DNA sequence is:", orig_dna_seq
    
    
def pick_random_AA( AA_seq, verbose = False ):
    # pick a random amino acid to mutate
    random_AA_choice = random.randrange( 0, len( AA_seq ) )
    
    # tell user, if desired
    if verbose:
        print "Mutating", AA_seq[ random_AA_choice ],
        print "at position", random_AA_choice
    
    return random_AA_choice



def get_codon_from_AA_position( dna_seq, AA_position ):
    # get the codon from the dna sequence given AA position
    codon_pos_start = AA_pos * 3
    codon_pos_end = codon_pos_start + 3
    orig_codon = dna[ codon_pos_start : codon_pos_end ]
    
    return orig_codon



def get_random_codon_given_AA( AA ):
    # make a random codon given an AA type
    rando_codon = random.choice( AA_dict[ AA ] )
    
    return rando_codon



def pick_random_codon_position( num_point_mutations ):
    # num_point_mutations should be 1, 2, or 3
    mutation_pos_choices = []
    ii = 0
    
    while ii < num_point_mutations:
        # pick a random position to mutate in the codon, only if it hasn't already been chosen
        random_pos_choice = random.randrange( 0, 3 )
        if random_pos_choice not in mutation_pos_choices:
            mutation_pos_choices.append( random_pos_choice )
            ii += 1
            
    return mutation_pos_choices



def mutate_codon( codon, mutation_pos_choices ):
    # for each position in the codon, mutate to a random nucleotide
    for mutation_position in mutation_pos_choices:
        current_single_nucleotide = codon[ mutation_position ]
        
        # get the three index values left of possible nucleotide_list to mutate to after removing current nucleotide as an option
        r = range( 4 )
        r.remove( nucleotide_list.index( current_single_nucleotide ) )
        
        # choose a random mutation and make it
        mutation = nucleotide_list[ random.choice( r ) ]
        # this is a round-about way of replacing a single nucleotide in the codon string
        codon = list( codon )
        codon[ mutation_position ] = mutation
        codon = "".join( codon )
        
    return codon



def random_codon_mutation( orig_codon, num_point_mutations, verbose = False ):
    # initialize a counter and some copies
    ii = 0
    jj = 0
    orig_AA = codon_dict[ orig_codon ]
    
    # it's a while True loop because some random codons could be non-coding, thus need to try again
    while True:
        # pick [a] random position[s] to mutate in the codon
        mutation_pos_choices = pick_random_codon_position( num_point_mutations )
        
        # mutate the codon as many times as specified
        codon = mutate_codon( orig_codon, mutation_pos_choices )
        
        # if mutation is non-coding, pick a random amino acid instead
        if codon_dict[ codon ] == '*':
            if verbose:
                print "    * Mutation led to a non-coding codon, picking a random mutation instead"
                
            while True:
                random_AA_choice = random.choice( AA_one_letter_list )
                if random_AA_choice != orig_AA:
                    codon = random.choice( AA_dict[ random_AA_choice ] )
                    break
            
            if verbose:
                print "    Original codon:", orig_codon,
                print "is now:", codon
            return codon
                
        # if mutation resulted in the original amino acid, pick a random amino acid instead
        if codon_dict[ codon ] == codon_dict[ orig_codon ]:
            if verbose:
                print "    * Mutation resulted in the original amino acid, picking a random mutation instead"
            while True:
                random_AA_choice = random.choice(AA_one_letter_list)
                if random_AA_choice != orig_AA:
                    codon = AA_dict[random_AA_choice][random.randrange( 0, len(AA_dict[random_AA_choice]) )]
                    break
            
            if verbose:
                print "    Original codon:", orig_codon,
                print "is now:", codon
            return codon
            
        # otherwise the codon is okay, so return it
        if verbose:
            print "    Original codon:", orig_codon,
            print "is now:", codon
            return codon
        
        
        
        
def mutation_time( dna ):
    # pick a random AA to mutate and return its codon
    rando_AA_choice_pos = pick_random_AA()
    codon_choice = get_codon_from_AA_sequence(dna, rando_AA_choice_pos)
    
    # make 1, 2, or 3 mutations in the codon, based on a random probability of each happening
    single_point = [1, 2, 3, 4, 5, 6] # 60% chaince
    double_point = [7, 8, 9]          # 30% chance
    triple_point = [10]               # 10% chance
    random_number = random.randint(1, 10)

    # make single, double, or triple point mutation based on the random number
    if random_number in single_point:
        #print "  Making a single point mutation"
        new_codon = random_codon_mutation(codon_choice, 1)
        
    if random_number in double_point:
        #print "  Making a double point mutation"
        new_codon = random_codon_mutation(codon_choice, 2)
        
    if random_number in triple_point:
        #print "  Making a triple point mutation"
        new_codon = random_codon_mutation(codon_choice, 3)
                
    # replace old codon with new one
    dna = list(dna)
    dna[codon_pos_start:codon_pos_end] = new_codon
    dna = "".join(dna)
    
    # return the amino acid to mutate to
    return codon_dict[ new_codon ]



def decide_my_mutation( AA ):
    # func(AA letter)
    codon_choice = get_random_codon_given_AA( AA )
    
    # make 1, 2, or 3 mutations in the codon, based on a random probability of each happening
    single_point = [1, 2, 3, 4, 5, 6] # 60% chaince
    double_point = [7, 8, 9]          # 30% chance
    triple_point = [10]               # 10% chance
    random_number = random.randint(1, 10)
    
    # make single, double, or triple point mutation based on the random number
    if random_number in single_point:
        #print "  Making a single point mutation"
        new_codon = random_codon_mutation(codon_choice, 1)
        
    if random_number in double_point:
        #print "  Making a double point mutation"
        new_codon = random_codon_mutation(codon_choice, 2)
            
    if random_number in triple_point:
        #print "  Making a triple point mutation"
        new_codon = random_codon_mutation(codon_choice, 3)
        
        
    # return the amino acid to mutate to
    return codon_dict[new_codon]
