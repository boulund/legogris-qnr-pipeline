from util import PathError

class FASTAParser:
    def __init__(self, logfile):
        self.logfile = logfile

    #TODO: gzip like FASTQ
    def parse(self, filename):
        return self.parse_fasta(filename)

    #for iterating
    def parse_fasta(self, filename):
        id = ''
        desc = ''
        tempseq = []
        try:
            seqfile = open(filename,'r')
            for line in seqfile:
                if line.startswith('>'):
                    if not id is '':
                        yield { 'id': id.strip(), 'desc': desc.strip(), 'dna': ''.join(tempseq) }
                    if ' ' in line:
                        (id, desc) = line[1::].split(' ', 1)
                    else:
                        id = line[1::].strip()
                        desc = ''
                    tempseq = []
                elif not line.startswith('>'):
                    tempseq.append(line.rstrip())
            if not id is '':
                yield { 'id': id.strip(), 'desc': desc.strip(), 'dna': ''.join(tempseq) }
        except OSError:
            raise PathError(''.join(['ERROR: cannot open', refseqpath]))
