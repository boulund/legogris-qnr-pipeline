import uuid
import json
import random
from collections import defaultdict
from cpython.ref cimport PyObject
#---------------------------------------------------------------------------
cdef extern from "Python.h":
    PyObject* PyString_FromString(char *v)
cdef extern from "stdlib.h":
  ctypedef unsigned int size_t
  size_t strlen(char *s)
  void *malloc(size_t size)
  void *calloc(size_t n, size_t size)
  void free(void *ptr)
  int strcmp(char *a, char *b)
  char * strcpy(char *a, char *b)
cdef char[89] _COMPLEMENTS
_COMPLEMENTS['A'] = 'T'
_COMPLEMENTS['C'] = 'G'
_COMPLEMENTS['G'] = 'C'
_COMPLEMENTS['T'] = 'A'
_COMPLEMENTS['Y'] = 'R'
_COMPLEMENTS['R'] = 'Y'
_COMPLEMENTS['S'] = 'S'
_COMPLEMENTS['W'] = 'W'
_COMPLEMENTS['K'] = 'M'
_COMPLEMENTS['M'] = 'K'
_COMPLEMENTS['B'] = 'V'
_COMPLEMENTS['V'] = 'B'
_COMPLEMENTS['D'] = 'H'
_COMPLEMENTS['H'] = 'D'
_COMPLEMENTS['N'] = 'N'
_COMPLEMENTS['X'] = 'X'
cdef char*** GT = <char ***>calloc(90, sizeof(char**))
GT['A'] = <char **>calloc(90, sizeof(char *))
GT['A']['A'] = <char *>calloc(90, sizeof(char))
GT['A']['G'] = <char *>calloc(90, sizeof(char))
GT['A']['C'] = <char *>calloc(90, sizeof(char))
GT['A']['T'] = <char *>calloc(90, sizeof(char))
GT['A']['R'] = <char *>calloc(90, sizeof(char))
GT['A']['Y'] = <char *>calloc(90, sizeof(char))
GT['A']['S'] = <char *>calloc(90, sizeof(char))
GT['A']['W'] = <char *>calloc(90, sizeof(char))
GT['A']['K'] = <char *>calloc(90, sizeof(char))
GT['A']['M'] = <char *>calloc(90, sizeof(char))
GT['A']['B'] = <char *>calloc(90, sizeof(char))
GT['A']['V'] = <char *>calloc(90, sizeof(char))
GT['A']['D'] = <char *>calloc(90, sizeof(char))
GT['A']['H'] = <char *>calloc(90, sizeof(char))
GT['A']['N'] = <char *>calloc(90, sizeof(char))
GT['A']['X'] = <char *>calloc(90, sizeof(char))
GT['G'] = <char **>calloc(90, sizeof(char *))
GT['G']['A'] = <char *>calloc(90, sizeof(char))
GT['G']['G'] = <char *>calloc(90, sizeof(char))
GT['G']['C'] = <char *>calloc(90, sizeof(char))
GT['G']['T'] = <char *>calloc(90, sizeof(char))
GT['G']['R'] = <char *>calloc(90, sizeof(char))
GT['G']['Y'] = <char *>calloc(90, sizeof(char))
GT['G']['S'] = <char *>calloc(90, sizeof(char))
GT['G']['W'] = <char *>calloc(90, sizeof(char))
GT['G']['K'] = <char *>calloc(90, sizeof(char))
GT['G']['M'] = <char *>calloc(90, sizeof(char))
GT['G']['B'] = <char *>calloc(90, sizeof(char))
GT['G']['V'] = <char *>calloc(90, sizeof(char))
GT['G']['D'] = <char *>calloc(90, sizeof(char))
GT['G']['H'] = <char *>calloc(90, sizeof(char))
GT['G']['N'] = <char *>calloc(90, sizeof(char))
GT['G']['X'] = <char *>calloc(90, sizeof(char))
GT['C'] = <char **>calloc(90, sizeof(char *))
GT['C']['A'] = <char *>calloc(90, sizeof(char))
GT['C']['G'] = <char *>calloc(90, sizeof(char))
GT['C']['C'] = <char *>calloc(90, sizeof(char))
GT['C']['T'] = <char *>calloc(90, sizeof(char))
GT['C']['R'] = <char *>calloc(90, sizeof(char))
GT['C']['Y'] = <char *>calloc(90, sizeof(char))
GT['C']['S'] = <char *>calloc(90, sizeof(char))
GT['C']['W'] = <char *>calloc(90, sizeof(char))
GT['C']['K'] = <char *>calloc(90, sizeof(char))
GT['C']['M'] = <char *>calloc(90, sizeof(char))
GT['C']['B'] = <char *>calloc(90, sizeof(char))
GT['C']['V'] = <char *>calloc(90, sizeof(char))
GT['C']['D'] = <char *>calloc(90, sizeof(char))
GT['C']['H'] = <char *>calloc(90, sizeof(char))
GT['C']['N'] = <char *>calloc(90, sizeof(char))
GT['C']['X'] = <char *>calloc(90, sizeof(char))
GT['T'] = <char **>calloc(90, sizeof(char *))
GT['T']['A'] = <char *>calloc(90, sizeof(char))
GT['T']['G'] = <char *>calloc(90, sizeof(char))
GT['T']['C'] = <char *>calloc(90, sizeof(char))
GT['T']['T'] = <char *>calloc(90, sizeof(char))
GT['T']['R'] = <char *>calloc(90, sizeof(char))
GT['T']['Y'] = <char *>calloc(90, sizeof(char))
GT['T']['S'] = <char *>calloc(90, sizeof(char))
GT['T']['W'] = <char *>calloc(90, sizeof(char))
GT['T']['K'] = <char *>calloc(90, sizeof(char))
GT['T']['M'] = <char *>calloc(90, sizeof(char))
GT['T']['B'] = <char *>calloc(90, sizeof(char))
GT['T']['V'] = <char *>calloc(90, sizeof(char))
GT['T']['D'] = <char *>calloc(90, sizeof(char))
GT['T']['H'] = <char *>calloc(90, sizeof(char))
GT['T']['N'] = <char *>calloc(90, sizeof(char))
GT['T']['X'] = <char *>calloc(90, sizeof(char))
GT['R'] = <char **>calloc(90, sizeof(char *))
GT['R']['A'] = <char *>calloc(90, sizeof(char))
GT['R']['G'] = <char *>calloc(90, sizeof(char))
GT['R']['C'] = <char *>calloc(90, sizeof(char))
GT['R']['T'] = <char *>calloc(90, sizeof(char))
GT['R']['R'] = <char *>calloc(90, sizeof(char))
GT['R']['Y'] = <char *>calloc(90, sizeof(char))
GT['R']['S'] = <char *>calloc(90, sizeof(char))
GT['R']['W'] = <char *>calloc(90, sizeof(char))
GT['R']['K'] = <char *>calloc(90, sizeof(char))
GT['R']['M'] = <char *>calloc(90, sizeof(char))
GT['R']['B'] = <char *>calloc(90, sizeof(char))
GT['R']['V'] = <char *>calloc(90, sizeof(char))
GT['R']['D'] = <char *>calloc(90, sizeof(char))
GT['R']['H'] = <char *>calloc(90, sizeof(char))
GT['R']['N'] = <char *>calloc(90, sizeof(char))
GT['R']['X'] = <char *>calloc(90, sizeof(char))
GT['Y'] = <char **>calloc(90, sizeof(char *))
GT['Y']['A'] = <char *>calloc(90, sizeof(char))
GT['Y']['G'] = <char *>calloc(90, sizeof(char))
GT['Y']['C'] = <char *>calloc(90, sizeof(char))
GT['Y']['T'] = <char *>calloc(90, sizeof(char))
GT['Y']['R'] = <char *>calloc(90, sizeof(char))
GT['Y']['Y'] = <char *>calloc(90, sizeof(char))
GT['Y']['S'] = <char *>calloc(90, sizeof(char))
GT['Y']['W'] = <char *>calloc(90, sizeof(char))
GT['Y']['K'] = <char *>calloc(90, sizeof(char))
GT['Y']['M'] = <char *>calloc(90, sizeof(char))
GT['Y']['B'] = <char *>calloc(90, sizeof(char))
GT['Y']['V'] = <char *>calloc(90, sizeof(char))
GT['Y']['D'] = <char *>calloc(90, sizeof(char))
GT['Y']['H'] = <char *>calloc(90, sizeof(char))
GT['Y']['N'] = <char *>calloc(90, sizeof(char))
GT['Y']['X'] = <char *>calloc(90, sizeof(char))
GT['S'] = <char **>calloc(90, sizeof(char *))
GT['S']['A'] = <char *>calloc(90, sizeof(char))
GT['S']['G'] = <char *>calloc(90, sizeof(char))
GT['S']['C'] = <char *>calloc(90, sizeof(char))
GT['S']['T'] = <char *>calloc(90, sizeof(char))
GT['S']['R'] = <char *>calloc(90, sizeof(char))
GT['S']['Y'] = <char *>calloc(90, sizeof(char))
GT['S']['S'] = <char *>calloc(90, sizeof(char))
GT['S']['W'] = <char *>calloc(90, sizeof(char))
GT['S']['K'] = <char *>calloc(90, sizeof(char))
GT['S']['M'] = <char *>calloc(90, sizeof(char))
GT['S']['B'] = <char *>calloc(90, sizeof(char))
GT['S']['V'] = <char *>calloc(90, sizeof(char))
GT['S']['D'] = <char *>calloc(90, sizeof(char))
GT['S']['H'] = <char *>calloc(90, sizeof(char))
GT['S']['N'] = <char *>calloc(90, sizeof(char))
GT['S']['X'] = <char *>calloc(90, sizeof(char))
GT['W'] = <char **>calloc(90, sizeof(char *))
GT['W']['A'] = <char *>calloc(90, sizeof(char))
GT['W']['G'] = <char *>calloc(90, sizeof(char))
GT['W']['C'] = <char *>calloc(90, sizeof(char))
GT['W']['T'] = <char *>calloc(90, sizeof(char))
GT['W']['R'] = <char *>calloc(90, sizeof(char))
GT['W']['Y'] = <char *>calloc(90, sizeof(char))
GT['W']['S'] = <char *>calloc(90, sizeof(char))
GT['W']['W'] = <char *>calloc(90, sizeof(char))
GT['W']['K'] = <char *>calloc(90, sizeof(char))
GT['W']['M'] = <char *>calloc(90, sizeof(char))
GT['W']['B'] = <char *>calloc(90, sizeof(char))
GT['W']['V'] = <char *>calloc(90, sizeof(char))
GT['W']['D'] = <char *>calloc(90, sizeof(char))
GT['W']['H'] = <char *>calloc(90, sizeof(char))
GT['W']['N'] = <char *>calloc(90, sizeof(char))
GT['W']['X'] = <char *>calloc(90, sizeof(char))
GT['K'] = <char **>calloc(90, sizeof(char *))
GT['K']['A'] = <char *>calloc(90, sizeof(char))
GT['K']['G'] = <char *>calloc(90, sizeof(char))
GT['K']['C'] = <char *>calloc(90, sizeof(char))
GT['K']['T'] = <char *>calloc(90, sizeof(char))
GT['K']['R'] = <char *>calloc(90, sizeof(char))
GT['K']['Y'] = <char *>calloc(90, sizeof(char))
GT['K']['S'] = <char *>calloc(90, sizeof(char))
GT['K']['W'] = <char *>calloc(90, sizeof(char))
GT['K']['K'] = <char *>calloc(90, sizeof(char))
GT['K']['M'] = <char *>calloc(90, sizeof(char))
GT['K']['B'] = <char *>calloc(90, sizeof(char))
GT['K']['V'] = <char *>calloc(90, sizeof(char))
GT['K']['D'] = <char *>calloc(90, sizeof(char))
GT['K']['H'] = <char *>calloc(90, sizeof(char))
GT['K']['N'] = <char *>calloc(90, sizeof(char))
GT['K']['X'] = <char *>calloc(90, sizeof(char))
GT['M'] = <char **>calloc(90, sizeof(char *))
GT['M']['A'] = <char *>calloc(90, sizeof(char))
GT['M']['G'] = <char *>calloc(90, sizeof(char))
GT['M']['C'] = <char *>calloc(90, sizeof(char))
GT['M']['T'] = <char *>calloc(90, sizeof(char))
GT['M']['R'] = <char *>calloc(90, sizeof(char))
GT['M']['Y'] = <char *>calloc(90, sizeof(char))
GT['M']['S'] = <char *>calloc(90, sizeof(char))
GT['M']['W'] = <char *>calloc(90, sizeof(char))
GT['M']['K'] = <char *>calloc(90, sizeof(char))
GT['M']['M'] = <char *>calloc(90, sizeof(char))
GT['M']['B'] = <char *>calloc(90, sizeof(char))
GT['M']['V'] = <char *>calloc(90, sizeof(char))
GT['M']['D'] = <char *>calloc(90, sizeof(char))
GT['M']['H'] = <char *>calloc(90, sizeof(char))
GT['M']['N'] = <char *>calloc(90, sizeof(char))
GT['M']['X'] = <char *>calloc(90, sizeof(char))
GT['B'] = <char **>calloc(90, sizeof(char *))
GT['B']['A'] = <char *>calloc(90, sizeof(char))
GT['B']['G'] = <char *>calloc(90, sizeof(char))
GT['B']['C'] = <char *>calloc(90, sizeof(char))
GT['B']['T'] = <char *>calloc(90, sizeof(char))
GT['B']['R'] = <char *>calloc(90, sizeof(char))
GT['B']['Y'] = <char *>calloc(90, sizeof(char))
GT['B']['S'] = <char *>calloc(90, sizeof(char))
GT['B']['W'] = <char *>calloc(90, sizeof(char))
GT['B']['K'] = <char *>calloc(90, sizeof(char))
GT['B']['M'] = <char *>calloc(90, sizeof(char))
GT['B']['B'] = <char *>calloc(90, sizeof(char))
GT['B']['V'] = <char *>calloc(90, sizeof(char))
GT['B']['D'] = <char *>calloc(90, sizeof(char))
GT['B']['H'] = <char *>calloc(90, sizeof(char))
GT['B']['N'] = <char *>calloc(90, sizeof(char))
GT['B']['X'] = <char *>calloc(90, sizeof(char))
GT['V'] = <char **>calloc(90, sizeof(char *))
GT['V']['A'] = <char *>calloc(90, sizeof(char))
GT['V']['G'] = <char *>calloc(90, sizeof(char))
GT['V']['C'] = <char *>calloc(90, sizeof(char))
GT['V']['T'] = <char *>calloc(90, sizeof(char))
GT['V']['R'] = <char *>calloc(90, sizeof(char))
GT['V']['Y'] = <char *>calloc(90, sizeof(char))
GT['V']['S'] = <char *>calloc(90, sizeof(char))
GT['V']['W'] = <char *>calloc(90, sizeof(char))
GT['V']['K'] = <char *>calloc(90, sizeof(char))
GT['V']['M'] = <char *>calloc(90, sizeof(char))
GT['V']['B'] = <char *>calloc(90, sizeof(char))
GT['V']['V'] = <char *>calloc(90, sizeof(char))
GT['V']['D'] = <char *>calloc(90, sizeof(char))
GT['V']['H'] = <char *>calloc(90, sizeof(char))
GT['V']['N'] = <char *>calloc(90, sizeof(char))
GT['V']['X'] = <char *>calloc(90, sizeof(char))
GT['D'] = <char **>calloc(90, sizeof(char *))
GT['D']['A'] = <char *>calloc(90, sizeof(char))
GT['D']['G'] = <char *>calloc(90, sizeof(char))
GT['D']['C'] = <char *>calloc(90, sizeof(char))
GT['D']['T'] = <char *>calloc(90, sizeof(char))
GT['D']['R'] = <char *>calloc(90, sizeof(char))
GT['D']['Y'] = <char *>calloc(90, sizeof(char))
GT['D']['S'] = <char *>calloc(90, sizeof(char))
GT['D']['W'] = <char *>calloc(90, sizeof(char))
GT['D']['K'] = <char *>calloc(90, sizeof(char))
GT['D']['M'] = <char *>calloc(90, sizeof(char))
GT['D']['B'] = <char *>calloc(90, sizeof(char))
GT['D']['V'] = <char *>calloc(90, sizeof(char))
GT['D']['D'] = <char *>calloc(90, sizeof(char))
GT['D']['H'] = <char *>calloc(90, sizeof(char))
GT['D']['N'] = <char *>calloc(90, sizeof(char))
GT['D']['X'] = <char *>calloc(90, sizeof(char))
GT['H'] = <char **>calloc(90, sizeof(char *))
GT['H']['A'] = <char *>calloc(90, sizeof(char))
GT['H']['G'] = <char *>calloc(90, sizeof(char))
GT['H']['C'] = <char *>calloc(90, sizeof(char))
GT['H']['T'] = <char *>calloc(90, sizeof(char))
GT['H']['R'] = <char *>calloc(90, sizeof(char))
GT['H']['Y'] = <char *>calloc(90, sizeof(char))
GT['H']['S'] = <char *>calloc(90, sizeof(char))
GT['H']['W'] = <char *>calloc(90, sizeof(char))
GT['H']['K'] = <char *>calloc(90, sizeof(char))
GT['H']['M'] = <char *>calloc(90, sizeof(char))
GT['H']['B'] = <char *>calloc(90, sizeof(char))
GT['H']['V'] = <char *>calloc(90, sizeof(char))
GT['H']['D'] = <char *>calloc(90, sizeof(char))
GT['H']['H'] = <char *>calloc(90, sizeof(char))
GT['H']['N'] = <char *>calloc(90, sizeof(char))
GT['H']['X'] = <char *>calloc(90, sizeof(char))
GT['N'] = <char **>calloc(90, sizeof(char *))
GT['N']['A'] = <char *>calloc(90, sizeof(char))
GT['N']['G'] = <char *>calloc(90, sizeof(char))
GT['N']['C'] = <char *>calloc(90, sizeof(char))
GT['N']['T'] = <char *>calloc(90, sizeof(char))
GT['N']['R'] = <char *>calloc(90, sizeof(char))
GT['N']['Y'] = <char *>calloc(90, sizeof(char))
GT['N']['S'] = <char *>calloc(90, sizeof(char))
GT['N']['W'] = <char *>calloc(90, sizeof(char))
GT['N']['K'] = <char *>calloc(90, sizeof(char))
GT['N']['M'] = <char *>calloc(90, sizeof(char))
GT['N']['B'] = <char *>calloc(90, sizeof(char))
GT['N']['V'] = <char *>calloc(90, sizeof(char))
GT['N']['D'] = <char *>calloc(90, sizeof(char))
GT['N']['H'] = <char *>calloc(90, sizeof(char))
GT['N']['N'] = <char *>calloc(90, sizeof(char))
GT['N']['X'] = <char *>calloc(90, sizeof(char))
GT['X'] = <char **>calloc(90, sizeof(char *))
GT['X']['A'] = <char *>calloc(90, sizeof(char))
GT['X']['G'] = <char *>calloc(90, sizeof(char))
GT['X']['C'] = <char *>calloc(90, sizeof(char))
GT['X']['T'] = <char *>calloc(90, sizeof(char))
GT['X']['R'] = <char *>calloc(90, sizeof(char))
GT['X']['Y'] = <char *>calloc(90, sizeof(char))
GT['X']['S'] = <char *>calloc(90, sizeof(char))
GT['X']['W'] = <char *>calloc(90, sizeof(char))
GT['X']['K'] = <char *>calloc(90, sizeof(char))
GT['X']['M'] = <char *>calloc(90, sizeof(char))
GT['X']['B'] = <char *>calloc(90, sizeof(char))
GT['X']['V'] = <char *>calloc(90, sizeof(char))
GT['X']['D'] = <char *>calloc(90, sizeof(char))
GT['X']['H'] = <char *>calloc(90, sizeof(char))
GT['X']['N'] = <char *>calloc(90, sizeof(char))
GT['X']['X'] = <char *>calloc(90, sizeof(char))

for a in 'AGCTRYSWKMBVDHNX':
    for b in 'AGCTRYSWKMBVDHNX':
        for c in 'AGCTRYSWKMBVDHNX':
            GT[<char>ord(a)][<char>ord(b)][<char>ord(c)] = ord('X')

GT['T']['T']['T'] = 'F'
GT['T']['T']['C'] = 'F'
GT['T']['T']['Y'] = 'F'

GT['T']['T']['A'] = 'L'
GT['T']['T']['G'] = 'L'
GT['T']['T']['R'] = 'L'

GT['C']['T']['T'] = 'L'
GT['C']['T']['C'] = 'L'
GT['C']['T']['A'] = 'L'
GT['C']['T']['G'] = 'L'
GT['C']['T']['M'] = 'L'
GT['C']['T']['R'] = 'L'
GT['C']['T']['W'] = 'L'
GT['C']['T']['S'] = 'L'
GT['C']['T']['Y'] = 'L'
GT['C']['T']['K'] = 'L'
GT['C']['T']['V'] = 'L'
GT['C']['T']['H'] = 'L'
GT['C']['T']['D'] = 'L'
GT['C']['T']['B'] = 'L'
GT['C']['T']['X'] = 'L'
GT['C']['T']['N'] = 'L'

GT['A']['T']['T'] = 'I'
GT['A']['T']['C'] = 'I'
GT['A']['T']['A'] = 'I'
GT['A']['T']['Y'] = 'I'
GT['A']['T']['W'] = 'I'
GT['A']['T']['M'] = 'I'
GT['A']['T']['H'] = 'I'
GT['A']['T']['G'] = 'M'

GT['G']['T']['T'] = 'V'
GT['G']['T']['C'] = 'V'
GT['G']['T']['A'] = 'V'
GT['G']['T']['G'] = 'V'
GT['G']['T']['M'] = 'V'
GT['G']['T']['R'] = 'V'
GT['G']['T']['W'] = 'V'
GT['G']['T']['S'] = 'V'
GT['G']['T']['Y'] = 'V'
GT['G']['T']['K'] = 'V'
GT['G']['T']['V'] = 'V'
GT['G']['T']['H'] = 'V'
GT['G']['T']['D'] = 'V'
GT['G']['T']['B'] = 'V'
GT['G']['T']['X'] = 'V'
GT['G']['T']['N'] = 'V'

GT['T']['C']['T'] = 'S'
GT['T']['C']['C'] = 'S'
GT['T']['C']['A'] = 'S'
GT['T']['C']['G'] = 'S'
GT['T']['C']['M'] = 'S'
GT['T']['C']['R'] = 'S'
GT['T']['C']['W'] = 'S'
GT['T']['C']['S'] = 'S'
GT['T']['C']['Y'] = 'S'
GT['T']['C']['K'] = 'S'
GT['T']['C']['V'] = 'S'
GT['T']['C']['H'] = 'S'
GT['T']['C']['D'] = 'S'
GT['T']['C']['B'] = 'S'
GT['T']['C']['X'] = 'S'
GT['T']['C']['N'] = 'S'

GT['C']['C']['T'] = 'P'
GT['C']['C']['C'] = 'P'
GT['C']['C']['A'] = 'P'
GT['C']['C']['G'] = 'P'
GT['C']['C']['M'] = 'P'
GT['C']['C']['R'] = 'P'
GT['C']['C']['W'] = 'P'
GT['C']['C']['S'] = 'P'
GT['C']['C']['Y'] = 'P'
GT['C']['C']['K'] = 'P'
GT['C']['C']['V'] = 'P'
GT['C']['C']['H'] = 'P'
GT['C']['C']['D'] = 'P'
GT['C']['C']['B'] = 'P'
GT['C']['C']['X'] = 'P'
GT['C']['C']['N'] = 'P'

GT['A']['C']['T'] = 'T'
GT['A']['C']['C'] = 'T'
GT['A']['C']['A'] = 'T'
GT['A']['C']['G'] = 'T'
GT['A']['C']['M'] = 'T'
GT['A']['C']['R'] = 'T'
GT['A']['C']['W'] = 'T'
GT['A']['C']['S'] = 'T'
GT['A']['C']['Y'] = 'T'
GT['A']['C']['K'] = 'T'
GT['A']['C']['V'] = 'T'
GT['A']['C']['H'] = 'T'
GT['A']['C']['D'] = 'T'
GT['A']['C']['B'] = 'T'
GT['A']['C']['X'] = 'T'
GT['A']['C']['N'] = 'T'

GT['G']['C']['T'] = 'A'
GT['G']['C']['C'] = 'A'
GT['G']['C']['A'] = 'A'
GT['G']['C']['G'] = 'A'
GT['G']['C']['M'] = 'A'
GT['G']['C']['R'] = 'A'
GT['G']['C']['W'] = 'A'
GT['G']['C']['S'] = 'A'
GT['G']['C']['Y'] = 'A'
GT['G']['C']['K'] = 'A'
GT['G']['C']['V'] = 'A'
GT['G']['C']['H'] = 'A'
GT['G']['C']['D'] = 'A'
GT['G']['C']['B'] = 'A'
GT['G']['C']['X'] = 'A'
GT['G']['C']['N'] = 'A'

GT['T']['A']['T'] = 'Y'
GT['T']['A']['C'] = 'Y'
GT['T']['A']['Y'] = 'Y'

GT['T']['A']['A'] = '*'
GT['T']['A']['G'] = '*'
GT['T']['A']['R'] = '*'

GT['C']['A']['T'] = 'H'
GT['C']['A']['C'] = 'H'
GT['C']['A']['Y'] = 'H'

GT['C']['A']['A'] = 'Q'
GT['C']['A']['G'] = 'Q'
GT['C']['A']['R'] = 'Q'

GT['A']['A']['T'] = 'N'
GT['A']['A']['C'] = 'N'
GT['A']['A']['Y'] = 'N'

GT['A']['A']['A'] = 'K'
GT['A']['A']['G'] = 'K'
GT['A']['A']['R'] = 'K'

GT['G']['A']['T'] = 'D'
GT['G']['A']['C'] = 'D'
GT['G']['A']['Y'] = 'D'
GT['G']['A']['A'] = 'E'
GT['G']['A']['G'] = 'E'
GT['G']['A']['R'] = 'E'

GT['T']['G']['T'] = 'C'
GT['T']['G']['C'] = 'C'
GT['T']['G']['Y'] = 'C'
GT['T']['G']['A'] = '*'
GT['T']['G']['G'] = 'W'

GT['C']['G']['T'] = 'R'
GT['C']['G']['C'] = 'R'
GT['C']['G']['A'] = 'R'
GT['C']['G']['G'] = 'R'
GT['C']['G']['M'] = 'R'
GT['C']['G']['R'] = 'R'
GT['C']['G']['W'] = 'R'
GT['C']['G']['S'] = 'R'
GT['C']['G']['Y'] = 'R'
GT['C']['G']['K'] = 'R'
GT['C']['G']['V'] = 'R'
GT['C']['G']['H'] = 'R'
GT['C']['G']['D'] = 'R'
GT['C']['G']['B'] = 'R'
GT['C']['G']['X'] = 'R'
GT['C']['G']['N'] = 'R'

GT['A']['G']['T'] = 'S'
GT['A']['G']['C'] = 'S'
GT['A']['G']['Y'] = 'S'

GT['A']['G']['A'] = 'R'
GT['A']['G']['G'] = 'R'
GT['A']['G']['R'] = 'R'

GT['G']['G']['T'] = 'G'
GT['G']['G']['C'] = 'G'
GT['G']['G']['A'] = 'G'
GT['G']['G']['G'] = 'G'
GT['G']['G']['M'] = 'G'
GT['G']['G']['R'] = 'G'
GT['G']['G']['W'] = 'G'
GT['G']['G']['S'] = 'G'
GT['G']['G']['Y'] = 'G'
GT['G']['G']['K'] = 'G'
GT['G']['G']['V'] = 'G'
GT['G']['G']['H'] = 'G'
GT['G']['G']['D'] = 'G'
GT['G']['G']['B'] = 'G'
GT['G']['G']['X'] = 'G'
GT['G']['G']['N'] = 'G'

cpdef bytes frame_sequence(char *sequence, int frame):
    cdef int i, j, k, l
    l = len(sequence)
    cdef char c
    cdef char* dseq = <char *>calloc(l + 1, sizeof(char))
    try:
        if frame > 2:
            #Reverse and frame adjust
            j = 0
            for i in range(l+2-frame, -1, -1):
                c = sequence[i]
                if c != 10: #skip newline
                #Complement
                    dseq[j] = _COMPLEMENTS[c]
                    j += 1
            dseq[j] = 0
        else:
            #Performance gain no biggie here
            j = 0
            for i in range(frame, l):
                c = sequence[i]
                if c != 10:
                    dseq[j] = c
                    j += 1
            dseq[j] = 0
        return dseq
    finally:
        free(dseq)


#Translates the supplied DNA string in all 6 reading frames and stores the result in a FASTA format text file as well as in serialized JSON in a supplied key/value store.
def translate_sequence(char *id, char *name, char *desc, char *sequence):
    global GT
    cdef char x = ord('X')
    cdef int i, j, frame, k, l
    cdef char *dna
    cdef char *protein
    cdef char c
    cdef char* s = NULL
    l = len(sequence)
    cdef char* dseq = <char *>calloc(l + 1, sizeof(char))
    cdef char* pseq = <char *>calloc(l/3 + 2, sizeof(char))
    if not dseq or not pseq:
        raise MemoryError()
    #Local variables = less overhead
    try:
        result = []
        for frame in range(0,6):
            #First 3 frames are normal, following 3 are reverse complements
            if frame > 2:
                #Reverse and frame adjust
                j = 0
                for i in range(l+2-frame, -1, -1):
                    c = sequence[i]
                    if c != 10: #skip newline
                    #Complement
                        dseq[j] = _COMPLEMENTS[c]
                        j += 1
                dseq[j] = 0
            else:
                #Performance gain no biggie here
                j = 0
                for i in range(frame, l):
                    c = sequence[i]
                    if c != 10:
                        dseq[j] = c
                        j += 1
                dseq[j] = 0
            dna  = dseq
            #Translate dna codons to protein
            j = 0
            k = len(dna)
            for i in range(0, k, 3):
                c = GT[<char>dna[i]][<char>dna[i+1] if i < k-1 else x ][<char>dna[i+2] if i < k-2 else x ]
                pseq[j] = c
                j += 1
            pseq[j] = 0
            protein = pseq
            seq = {
                'protein': protein,
                'name': name,
                'description': desc,
                'frame': frame
            }
            out = ''.join(['>', id, '_', str(frame), '\n', protein, '\n'])
            try:
                result.append((frame, json.dumps(seq), out))
            except:
                print("* dna: ", dna)
                print(seq)
                exit(1)
    finally:
        free(dseq)
        free(pseq)
    return result
