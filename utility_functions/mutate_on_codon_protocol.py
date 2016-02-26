from mutate_on_codon_functions import *
import argparse
import os
import sys
import pandas as pd

parser = argparse.ArgumentParser( description="Use Python to make codon mutations based off of a DNA sequence" )
parser.add_argument( "dna_seq_file", type=str, help="the filename of the DNA sequence to be mutated" )
#parser.add_argument( "single_point_mutation_probability", type=float, help="how likely do you want a single point codon mutation to be? Give me a float from 0 to 100")
#parser.add_argument( "double_point_mutation_probability", type=float, help="how likely do you want a double point codon mutation to be? Give me a float from 0 to 100")
#parser.add_argument( "triple_point_mutation_probability", type=float, help="how likely do you want a triple point codon mutation to be? Give me a float from 0 to 100")
input_args = parser.parse_args()


# check the validity of the dna_seq_file
if not os.path.isfile( input_args.dna_seq_file ):
    print "You gave me an invalid path to your DNA sequence file. Exiting."
    sys.exit()
else:
    dna_seq_file = input_args.dna_seq_file

# make sure the probability numbers passed are from 0 to 100
'''
if not ( 0 <= input_args.single_point_mutation_probability <= 100 ):
    print "You gave me an inappropriate probability for your single point mutation likelihood. Exiting."
    sys.exit()
else:
    single_point_prob = float( input_args.single_point_mutation_probability / 100 )
if not ( 0 <= input_args.double_point_mutation_probability <= 100 ):
    print "You gave me an inappropriate probability for your double point mutation likelihood. Exiting."
    sys.exit()
else:
    double_point_prob = float( input_args.double_point_mutation_probability / 100 )
if not ( 0 <= input_args.triple_point_mutation_probability <= 100 ):
    print "You gave me an inappropriate probability for your triple point mutation likelihood. Exiting."
    sys.exit()
else:
    triple_point_prob = float( input_args.triple_point_mutation_probability / 100 )

# check to make sure that the single, double, and triple point mutation likelihood adds up to 100
if ( input_args.single_point_mutation_probability + input_args.double_point_mutation_probability + input_args.triple_point_mutation_probability ) != 100:
    print "Your probability of single, double, and triple point mutations didn't add up to 100. Exiting."
    sys.exit()
'''

# read the file to get the DNA sequence
dna_seq = []
with open( dna_seq_file, 'rb' ) as fh:
    lines = fh.readlines()
    for line in lines:
        line = line.strip()
        for char in line:
            dna_seq.append( char )
# turn the list into a string
dna_seq = ''.join( dna_seq )

# check to make sure the length of the DNA is divisible by 3
if len( dna_seq ) % 3 != 0:
    print "Your DNA sequence is not a multiple of three. Exiting."
    sys.exit()


# prepare a Pandas DataFrame to hold the data
# columns = three-letter code for each amino acid possibility
columns = AA_three_letter_list
# holds the data to be returned
df = pd.DataFrame( columns = columns)

# for each codon, determine the likelihood of that amino acid mutating into the other 19
codon_start_positions = range( 0, len( dna_seq ), 3 )
for codon_start in codon_start_positions:
    codon_stop = codon_start + 3
    codon = dna[ codon_start : codon_stop ]
    
    
