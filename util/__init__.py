from itertools import izip_longest

def fixfasta(sequence):
    """
    Take a FASTA-formatted sequence and tries to correct its
    formatting by splitting it on lines of maximum length 80.

    Args::

        sequence (str)  FASTA sequence, in one
                        complete string with \\n markers
                        between identifier line and sequence.

    Returns::

        sequence (str)  Sequence with hopefully
                        better formatting.
    """

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

def sequence_to_fasta(id, sequence):
    """
    Returns the supplied `id` and `sequence` (both strings) as a sequence string in FASTA format.
    """
    return fixfasta(''.join(['>', id, '\n', sequence, '\n']))

def grouper(iterable, n, fillvalue=None):
    """
    Partition data into fixed-length tuples. Recipe from itertools docs.

    >>> grouper('ABCDEFG', 3, 'x')
    [('A', 'B', 'C'),
    ('D', 'E', 'F'),
    ('G', 'x', 'x')]
    """
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

class PathError(Exception):
    """Raised when an invalid or otherwise not usable path has been supplied."""
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)

class ParseError(Exception):
    """Raised when an error has occured while parsing an input string or file."""
    def __init__(self,message):
        self.message = message
    def __str__(self):
        return repr(self.message)
