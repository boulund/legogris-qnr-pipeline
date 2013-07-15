# Lots of useful functions, Fredrik Boulund 2011
# read_fasta
# format_fasta

__doc__ = """   Fredrik Boulund function module 2011
    Read individual function descriptions for more details"""


##---------------------------------------------------------------------------##
##                      READ FASTA FILE INTO PYTHON                          ##
##---------------------------------------------------------------------------##
def read_fasta(filename):
    """
    Takes a filename and reads FASTA sequence information.

    Stores information in a list with (sequenceID,sequence)-tuples.

    Input:
        filename    name or path of the FASTA file to be read.
   
    Output:
        seqlist     list of tuples with sequenceID and sequence
                    read from the input FASTA file.
                    Note that sequence_identifier does not end in 
                    a newline character! This needs to be added 
                    to print sequences in correct format.
    
    Errors:
        (none)

    Author: Fredrik Boulund, 2011
    """

    fasta_file = open(filename) # Open file at filename
    line = fasta_file.readline() # Read first line
    first = True # Boolean determining if it is the first sequence in file
    sequence = "" # Stores sequence information
    seqlist = [] # List to store tuples of seqid,sequence information
    while fasta_file:
        if first and line.startswith(">"):
            # If first occurrence of ">", treat it a bit different
            sequence_identifier = line.rstrip()
            line = fasta_file.readline()
            first = False
        
        elif line.startswith(">") and not first:
            # Store complete (previous) sequence and id in tuple into list
            sequence = "".join(sequence.split('\n')) # Remove embedded newlines
            seqlist.append((sequence_identifier,sequence)) 
            
            sequence = ""  # Reset sequence string
            sequence_identifier = line.rstrip() # Set new seqid
            line = fasta_file.readline()
        
        elif line == "": #EOF
            # Store complete (last) sequence and id in tuple into list
            sequence = "".join(sequence.split('\n')) # Remove embedded newlines
            seqlist.append((sequence_identifier,sequence))
            break
        
        else: # This is sequence territory
            sequence = ''.join([sequence,line])
            line = fasta_file.readline()

    return seqlist
######################### END read_fasta




##-----------------------------------------------##
##               FORMAT FASTA FILES              ##
##-----------------------------------------------##
def format_fasta(sequences, linelength=80 ):
    '''
    Takes a list with sequences (output from read_fasta) and restructures
    the output, limiting the line length to whatever the user desires.
    
    Input:
        sequences   list of sequences, each in one
                    complete string with \n markers 
                    between identifier line and sequence.
        linelength  an integer with the maximum number of 
                    columns per row.
    Output:
        outsequences   list of sequences with hopefully
                       better formatting.
    Errors:
         (none)     probably a lot... :)
    '''

    from math import ceil

    linelength = int(linelength)

    # Holds the output sequences
    outsequences = []
    for identifier, sequence in sequences:
         
        # Determine how to split the sequence string
        number_of_rows = int(ceil(len(sequence) / float(linelength)))
        
        seq = []
       
        if number_of_rows > 1:
            seq.append(identifier)
            for row in xrange(0,number_of_rows):
                seq.append(sequence[:linelength])
                sequence = sequence[linelength:]
            seq.append(sequence)
            seq = '\n'.join(seq)
        else:
            seq = '\n'.join([identifier,sequence])
        outsequences.append(seq)

    return outsequences
############## END format_fasta


if __name__ == "__main__": 
    print "Unit testing"+"-"*25
    print "Nothing to test"



##OOOOOOOOOOOOOOOOOOOOOOLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLDDDDDDDDDDDDDDDDDDDDDD

###---------------------------------------------------------------------------##
###                      CLASSIFY SEQUENCES AS QNR                            ##
###---------------------------------------------------------------------------##
#def classify_qnr(sequence_length, domain_score, func="", longseqcutoff=55, longseqthresh=85):
#    """
#    Classifies a sequence as Qnr or not.

#    Using the domain_score and a hardcoded (or user defined function)
#    to classify a given sequence as Qnr or not.

#    Input:
#        sequence        a string with sequence information.
#        domain_score    a float with the domain score for this sequence.
#        func            an optional function to determine classification.
#        longseqthresh   the threshold score cutoff for long sequences
#    Output:
#        classification  a boolean determining whether it should be classified
#                        as Qnr or not.
#    Errors:
#        (none)
#    """
#    
#    # Define a hardcoded function if none given
#    if func == "":
#        func = lambda L: 0.41*L + 18
#    
#    # Pretty self-explanatory. Has a range in which the classification
#    # function is used, determined by the first if-statement 
#    if (int(sequence_length) >= longseqthresh) and (float(domain_score) > longseqcutoff):
#        return True
#    elif int(sequence_length) < longseqthresh:
#        if float(domain_score) > func(float(sequence_length)):
#            return True
#        else:
#            return False
#    else:
#        return False

########################## END classify_qnr




