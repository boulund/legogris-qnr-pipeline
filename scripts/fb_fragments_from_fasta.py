#!/usr/bin/env python
# Fredrik Boulund 20110630

from sys import argv, exit
from os import path, system
from random import randint, random, choice, sample
from FB_functions import read_fasta
from optparse import OptionParser
import re


desc="""Script cuts sequences into a multitude of pieces.
\nFredrik Boulund 20110630"""

usage = "script.py [options] sequences.fasta"

parser = OptionParser(usage=usage,description=desc)

parser.add_option("-p", "--probability", dest="probability", type="float",
                  default=0, help="set probability of amino acid mutation, 0<=p<1 [default=%default]")
parser.add_option("-r", "--replicates", dest="replicates",
                  help="number of replicates [default=%default]",type="int", default=1000)
parser.add_option("-l", "--fraglength", dest="fraglength", type="int",
                  help="maximum fragment length, will be reduced to longest sequence in source file if too long [default=%default]",
                  default=210)
parser.add_option("-f", "--fixedfraglength", dest="fixed", type="int",
                  help="Set a fixed fraglength to create fragments of [default=%default]",
                  default=False)
parser.add_option("-o", "--outputfile", dest="outputfile", type="string", default="fragments.pfa",
                  help="Filename for output file")
parser.add_option("-t", "--sequence-type", dest="sequence_type", type="string", default="protein",
                  help="Sequence input type [dna/protein]")

(options, args) = parser.parse_args()

if len(args)<1:
    print parser.print_help()
    print "ERROR: Need fasta file to work on!"
    exit()



# Read QNR-sequences from fasta and store sequences in list of tuples
seqlist = read_fasta(args[0])

##---------------------------------------------------------------------------##
##                       DATA FOR FRAGMENT CREATION                          ##
##---------------------------------------------------------------------------##

# The number of different fragment lengths to be created:
MAXIMUM_FRAGMENT_LENGTH = options.fraglength  # max value is ~210 (213) for Qnr

# Number of fragments per fragment length:
NUMBER_OF_FRAGMENTS_PER_LENGTH = options.replicates

# Determine the minimum sequence length
# (sets upper limit for fragment starting point)
SEQUENCE_MINLENGTH = 5000 # sufficiently large number
for seqid,sequence in seqlist:
    if len(sequence) < SEQUENCE_MINLENGTH:
        SEQUENCE_MINLENGTH = len(sequence)

# SEQUENCE_MINLENGTH is now not longer than minimum sequence length in src file
if MAXIMUM_FRAGMENT_LENGTH > SEQUENCE_MINLENGTH:
    MAXIMUM_FRAGMENT_LENGTH = SEQUENCE_MINLENGTH-1

# Read probability of per amino acid error from command line option:
if options.probability > 1.0:
    print "Probability was incorrectly entered, needs 0<p<1"
    exit()
else:
    P_MUTATION = options.probability

# Set of valid amino acid residues:
if options.sequence_type == 'dna':
    aa_list = set(['A', 'T', 'C', 'G'])
elif options.sequence_type == 'protein':
    aa_list = set(['A','C','D','E','F','G','H','I','K','L','M',
                'N','P','Q','R','S','T','V','W','Y'])
else:
    print "Illegal sequence type. Possible values are dna and protein"

if options.fixed:
    fragment_lengths = [options.fixed]
else:
    # Create a list of fragment sizes, starting from 1 amino acid
    # up to MAXIMUM_FRAGMENT_LENGTH using stepsize 1 amino acid.
    fragment_lengths = range(10, MAXIMUM_FRAGMENT_LENGTH, 1)


##---------------------------------------------------------------------------##
##                              FRAGMENT CREATION                            ##
##---------------------------------------------------------------------------##

# Output file
outfile = open(options.outputfile,"w")

# Loop over each fragment length
for fragment_length in fragment_lengths:

    # Loop over current fragment length
    for fragment_number in xrange(1,NUMBER_OF_FRAGMENTS_PER_LENGTH+1):

        # Randomly select one source sequence from the seqlist
        random_src_seq = randint(0,len(seqlist)-1)


        # Unzip the current (randomly selected) sequence
        cur_seqid, cur_seq = seqlist[random_src_seq]
        cur_seqlength = len(cur_seq)
        # Create random (valid) starting point within sequence length
        random_starting_point = randint(0, cur_seqlength - fragment_length)
        # Compute the endpoint according to current fragment length
        endpoint = random_starting_point + fragment_length
        # Cut fragment from random_starting_point -> endpoint
        cur_fragment = cur_seq[random_starting_point:endpoint]

        # If mismatches are wanted; for each amino acid in fragment
        # throw random number and if less than P_MUTATION we
        # randomly select one amino acid residue in its place
        # If P_MUTATION == 0 there can never be mutations
        if P_MUTATION > 0.0:
            cur_fragment_mut = ''
            for aa in cur_fragment:
                utfall = random()
                if utfall < P_MUTATION:
                    # Sample a random aa that is not the current
                    cur_fragment_mut = ''.join([cur_fragment_mut,
                                            ''.join((sample(aa_list-set(aa),1)))])
                else:
                    # No mutation occurs
                    cur_fragment_mut = ''.join([cur_fragment_mut, aa])
            cur_fragment = cur_fragment_mut

        # Write fragment to disk with new identifier:
        # >fragment_length__random_starting_point__sequenceID
        outfile.write(''.join([">",str(fragment_length),"__",str(random_starting_point),"__",cur_seqid[1:],"\n",
                       cur_fragment,"\n"]))

# Close file
outfile.close()

