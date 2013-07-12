#!/usr/bin/python
from parsing.fastq import FASTQParser

parser = FASTQParser(None)
for seq in parser.parse_gzip('tutorial/database/india2.fastq.gz'):
    pass
    print seq
