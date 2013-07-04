from util import PathError

class FASTAParser:
    def __init__(self, logfile):
        self.logfile = logfile

    #Reinserts sequences with hmm score and dscore in indb
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
