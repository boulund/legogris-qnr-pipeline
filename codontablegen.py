#!/bin/python
import sys
#TODO: Loop over codon table and generate dict with strings, should be miuch faster than regex
f = open('codons.txt', 'r')
for l in f:
    s = l.split('"')
    (rex, acid) = s[1].split(' ')
    (first, last) = rex.split('[')
    last = last[:-1]
    for c in last:
        sys.stdout.write("'%s%s': '%s'," % (first, c, acid))
    sys.stdout.write("\n")


