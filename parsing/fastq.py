from util import PathError
from itertools import izip_longest
import zlib

class FASTQParser:
    def __init__(self, logfile, gzip=False):
        self.logfile = logfile
        self.gzip = gzip

    def parse(self, filename):
        if self.gzip:
            return self.parse_gzip(filename)
        return self.parse_fastq(filename)

    def parse_gzip(self, filename):
        gzfile = open(filename, 'rb')
        d = zlib.decompressobj(16+zlib.MAX_WBITS)
        eof = False
        data = ''
        while True:
            while True:
                cdata = gzfile.read(2000000)
                if not cdata:
                    eof = True
                    break
                data = ''.join([data, d.decompress(cdata)])
                lines = data.split('\n')
                num = len(lines) - 1
                if num >= 4:
                    slines = lines[0:num - (num%4)]
                    data = '\n'.join(lines[num - (num%4)::])
                    break
            if eof:
                break
            for ls in grouper(slines, 4):
                yield self.parse_lines(ls)

    #for iterating
    def parse_fastq(self, filename):
        id = ''
        desc = ''
        tempseq = []
        try:
            seqfile = open(filename,'r')
            for lines in grouper(seqfile, 4):
                yield self.parse_lines(lines)
        except OSError:
            raise PathError(''.join(['ERROR: cannot open', refseqpath]))

    def parse_lines(self, lines):
        idline = lines[0]
        dna = lines[1]
        #optional: quality = lines[3]
        if ' ' in lines[0]:
            (id, desc) = lines[0][1::].split(' ', 1)
            return { 'id': id.strip(), 'desc': desc.strip(), 'dna': dna.strip() }
        else:
            id = lines[0][1::]
            return { 'id': id.strip(), 'dna': dna.strip() }

#Collect data into fixed-length chunks or blocks. Recipe from itertools docs.
# grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)
