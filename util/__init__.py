def fixfasta(sequence):

    from math import ceil

    splitstring = sequence.split("\n",1)
    number_of_rows = int(ceil(len(splitstring[1]) / 80.0))
    seq = []
    if number_of_rows > 1:
        for row in xrange(1,number_of_rows):
            seq.append(splitstring[1][:80] + "\n")
            splitstring[1] = splitstring[1][80:]
        seq.append(splitstring[1]) # + "\n")
        seq.insert(0,splitstring[0] + "\n")
        return ''.join(seq)
    else:
        return sequence

##-----------------------------------------------##
##            FIX FASTA FORMATTING               ##
##-----------------------------------------------##
def fixfastas(sequences):
    '''
    Takes a list of sequences and tries to correct their
    formatting.

    Input::

        sequences   list of sequences, each in one
                    complete string with \n markers
                    between identifier line and sequence.

    Returns::

        outsequences   list of sequences with hopefully
                    better formatting.

    Errors::

        (none)
    '''

    return [fixfasta(seq) for seq in sequences]
############## END fixfasta

def sequence_to_fasta(id, sequence):
    return fixfasta(''.join(['>', id, '\n', sequence, '\n']))

##-----------------------------------------------##
##                CUSTOM EXCEPTIONS              ##
##-----------------------------------------------##

# PathError       General error with path
class PathError(Exception):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class ParseError(Exception):
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)