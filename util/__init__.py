##-----------------------------------------------##
##            FIX FASTA FORMATTING               ##
##-----------------------------------------------##
def fixfastas(sequences):
    '''
    Takes a list of sequences and tries to correct their
    formatting.

    Designed to be used only in the
    the function retrieve_sequences_from_hmmsearch.

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

def fragment_to_fasta(fragment, id, sequence='protein'):
    seq = (fragment['dna'] if sequence == 'dna' else fragment['protein'])
    return fixfasta(''.join(['>', id, '\n', seq, '\n']))

def parse_fasta(filename):
    sequences = []
    tempseqid = ''
    tempseq = []
    try:
        seqfile = open(filename,'r')
        for line in seqfile:
            if line.startswith('>'):
                seq = { 'id': tempseqid, 'dna': ''.join(tempseq) }
                sequences.append(seq)
                tempseqid = line[1::]
                tempseq = []
            elif not line.startswith(">"):
                tempseq.append(line)
        seq = { 'id': tempseqid, 'dna': ''.join(tempseq) }
        sequences.append(seq)
        return sequences
    except OSError:
        raise PathError(''.join(['ERROR: cannot open', refseqpath]))
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
