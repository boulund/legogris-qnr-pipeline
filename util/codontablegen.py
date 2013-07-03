#!/bin/python
import sys
def oldtablegen():
    f = open('codons.txt', 'r')
    for l in f:
        s = l.split('"')
        (rex, acid) = s[1].split(' ')
        (first, last) = rex.split('[')
        last = last[:-1]
        for c in last:
            sys.stdout.write("'%s%s': '%s'," % (first, c, acid))
        sys.stdout.write("\n")

def tablegen():
    for a in 'AGCTRYSWKMBVDHNX':
        print("GT['%s'] = <char **>calloc(89, sizeof(char *))" % a)
        for b in 'AGCTRYSWKMBVDHNX':
            print("GT['%s']['%s'] = <char *>calloc(89, sizeof(char))" % (a, b))
tablegen()
