from util import PathError, grouper
import zlib

class FASTQParser:
    def __init__(self, logfile, gzip=False):
        """
        Kwargs:
            gzip (bool): If set to True, all calls to parse will first attempt to perform gzip decompression on the file.
        """
        self.logfile = logfile
        self.gzip = gzip

    def parse(self, filename):
        """Perform `parse_gzip` or `parse_fastq` depending on if gzip decompression is enabled or not."""
        if self.gzip:
            return self.parse_gzip(filename)
        return self.parse_fastq(filename)

    def parse_gzip(self, filename):
        """
        Translate input in file at path defined by `filename` from FASTQ format after gzip decompression.

        Returns:
            List of dictionaries with the following keys:
                * id (str): Identifier for the sequence.
                * desc (str, optional): If the original identifier contains at least one space, desc will be the string after the last one.
                * dna (str): The sequence DNA string.
        """
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
                yield self._parse_lines(ls)

    #for iterating
    def parse_fastq(self, filename):
        """
        Translate input in file at path defined by `filename` from FASTQ format.

        Returns:
            List of dictionaries with the following keys:
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
            for lines in grouper(seqfile, 4):
                yield self._parse_lines(lines)
        except OSError:
            raise PathError(''.join(['ERROR: cannot open', refseqpath]))

    def _parse_lines(self, lines):
        idline = lines[0]
        dna = lines[1]
        #optional: quality = lines[3]
        if ' ' in lines[0]:
            (id, desc) = lines[0][1::].split(' ', 1)
            return { 'id': id.strip(), 'desc': desc.strip(), 'dna': dna.strip() }
        else:
            id = lines[0][1::]
            return { 'id': id.strip(), 'dna': dna.strip() }

