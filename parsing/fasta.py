from util import PathError

class FASTAParser:
    def __init__(self, logfile):
        self.logfile = logfile


    #for iterating
    def parse_fasta(self, filename):
        tempseqid = ''
        tempseq = []
        try:
            seqfile = open(filename,'r')
            for line in seqfile:
                if line.startswith('>'):
                    if tempseq or tempseqid:
                        yield { 'id': tempseqid, 'dna': ''.join(tempseq) }
                    tempseqid = line[1::].strip()
                    tempseq = []
                elif not line.startswith(">"):
                    tempseq.append(line.strip())
            if tempseq or tempseqid:
                yield { 'id': tempseqid, 'dna': ''.join(tempseq) }
        except OSError:
            raise PathError(''.join(['ERROR: cannot open', refseqpath]))

    #Returns list of sequences
    def fasta_to_list(self, filename):
        return [f for f in self.parse_fasta(filename)]
