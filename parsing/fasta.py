from util import PathError

class FASTAParser:
    def __init__(self, logfile):
        self.logfile = logfile

    #TODO: gzip like FASTQ
    def parse(self, filename):
        """
        Calls parse_fasta.
        """
        return self.parse_fasta(filename)

    def parse_fasta(self, filename):
        """
        Translate input in file at path defined by `filename` from FASTA format.

        Returns:
            Generator function of dictionaries with the following keys:
                * id (str): Identifier for the sequence.
                * desc (str, optional): If the original identifier contains at least one space, desc will be the string after the last one.
                * dna (str): The sequence DNA string.

        Raises:
            PathError
        """
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
