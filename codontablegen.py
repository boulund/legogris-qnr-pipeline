#!/bin/python
#TODO: Loop over codon table and generate dict with strings, should be miuch faster than regex
f = open('codons.txt', 'r')
for l in f:
    s = l.split('"')
    (rex, acid) = s[1].split(' ')
    print rex

