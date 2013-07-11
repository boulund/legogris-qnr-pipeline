from util import PathError
from itertools import grouper

class FASTAParser:
    def __init__(self, logfile):
        self.logfile = logfile

    #for iterating
    def parse_fasta(self, filename):
        id = ''
        desc = ''
        tempseq = []
        try:
            seqfile = open(filename,'r')
            for lines in grouper(seqfile):
                idline = lines[0]
                dna = lines[1]
                #optional: quality = lines[3]
                if ' ' in lines[0]:
                    (id, desc) = lines[0][1::].split(' ', 1)
                    yield { 'id': id.strip(), 'desc': desc.strip(), 'dna': dna.strip() }
                else:
                    id = lines[0][1::]
                    yield { 'id': id.strip(), 'dna': dna.strip() }
        except OSError:
            raise PathError(''.join(['ERROR: cannot open', refseqpath]))

    #Returns list of sequences
    def fasta_to_list(self, filename):
        return [f for f in self.parse_fasta(filename)]

    #Collect data into fixed-length chunks or blocks. Recipe from itertools docs.
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    def grouper(iterable, n, fillvalue=None):
        args = [iter(iterable)] * n
        return izip_longest(fillvalue=fillvalue, *args)
