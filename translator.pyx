import uuid
import json
import random
from collections import defaultdict
#---------------------------------------------------------------------------
cdef extern from "stdlib.h":
  ctypedef unsigned int size_t
  size_t strlen(char *s)
  void *malloc(size_t size)
  void *calloc(size_t n, size_t size)
  void free(void *ptr)
  int strcmp(char *a, char *b)
  char * strcpy(char *a, char *b)
cdef int[89] _COMPLEMENTS
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
_GENCODE = defaultdict(lambda: 'X', {
    'TTT': 'F','TTC': 'F','TTY': 'F',
    'TTA': 'L','TTG': 'L','TTR': 'L',
    'CTT': 'L','CTC': 'L','CTA': 'L','CTG': 'L','CTM': 'L','CTR': 'L','CTW': 'L','CTS': 'L','CTY': 'L','CTK': 'L','CTV': 'L','CTH': 'L','CTD': 'L','CTB': 'L','CTX': 'L','CTN': 'L',
    'ATT': 'I','ATC': 'I','ATA': 'I','ATY': 'I','ATW': 'I','ATM': 'I','ATH': 'I', 'ATG': 'M',
    'GTT': 'V','GTC': 'V','GTA': 'V','GTG': 'V','GTM': 'V','GTR': 'V','GTW': 'V','GTS': 'V','GTY': 'V','GTK': 'V','GTV': 'V','GTH': 'V','GTD': 'V','GTB': 'V','GTX': 'V','GTN': 'V',
    'TCT': 'S','TCC': 'S','TCA': 'S','TCG': 'S','TCM': 'S','TCR': 'S','TCW': 'S','TCS': 'S','TCY': 'S','TCK': 'S','TCV': 'S','TCH': 'S','TCD': 'S','TCB': 'S','TCX': 'S','TCN': 'S',
    'CCT': 'P','CCC': 'P','CCA': 'P','CCG': 'P','CCM': 'P','CCR': 'P','CCW': 'P','CCS': 'P','CCY': 'P','CCK': 'P','CCV': 'P','CCH': 'P','CCD': 'P','CCB': 'P','CCX': 'P','CCN': 'P',
    'ACT': 'T','ACC': 'T','ACA': 'T','ACG': 'T','ACM': 'T','ACR': 'T','ACW': 'T','ACS': 'T','ACY': 'T','ACK': 'T','ACV': 'T','ACH': 'T','ACD': 'T','ACB': 'T','ACX': 'T','ACN': 'T',
    'GCT': 'A','GCC': 'A','GCA': 'A','GCG': 'A','GCM': 'A','GCR': 'A','GCW': 'A','GCS': 'A','GCY': 'A','GCK': 'A','GCV': 'A','GCH': 'A','GCD': 'A','GCB': 'A','GCX': 'A','GCN': 'A',
    'TAT': 'Y','TAC': 'Y','TAY': 'Y',
    'TAA': '*','TAG': '*','TAR': '*',
    'CAT': 'H','CAC': 'H','CAY': 'H',
    'CAA': 'Q','CAG': 'Q','CAR': 'Q',
    'AAT': 'N','AAC': 'N','AAY': 'N',
    'AAA': 'K','AAG': 'K','AAR': 'K',
    'GAT': 'D','GAC': 'D','GAY': 'D',
    'GAA': 'E','GAG': 'E','GAR': 'E',
    'TGT': 'C','TGC': 'C','TGY': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R','CGC': 'R','CGA': 'R','CGG': 'R','CGM': 'R','CGR': 'R','CGW': 'R','CGS': 'R','CGY': 'R','CGK': 'R','CGV': 'R','CGH': 'R','CGD': 'R','CGB': 'R','CGX': 'R','CGN': 'R',
    'AGT': 'S','AGC': 'S','AGY': 'S',
    'AGA': 'R','AGG': 'R','AGR': 'R',
    'GGT': 'G','GGC': 'G','GGA': 'G','GGG': 'G','GGM': 'G','GGR': 'G','GGW': 'G','GGS': 'G','GGY': 'G','GGK': 'G','GGV': 'G','GGH': 'G','GGD': 'G','GGB': 'G','GGX': 'G','GGN': 'G',
    'CT': 'L', 'GT': 'V', 'TC': 'S', 'CC': 'P', 'AC': 'T', 'GC': 'A', 'CG': 'R', 'GG': 'G'
})
#Translates the supplied DNA string in all 6 reading frames and stores the result in a FASTA format text file as well as in serialized JSON in a supplied key/value store.
def translate_sequence(char *name, char *desc, char *sequence):
    cdef int i, j, frame, l
    cdef char *dna
    cdef char *protein
    cdef char c
    cdef char* s = NULL
    l = len(sequence)
    cdef char* dseq = <char *>calloc(l + 1, sizeof(char))
    #Local variables = less overhead
    result = []
    gencode = _GENCODE
    for frame in range(0,6):
        #First 3 frames are normal, following 3 are reverse complements
        if frame > 2:
            #Reverse and frame adjust
            j = 0
            for i in range(l-1, frame-4, -1):
                c = sequence[i]
                if c != 10: #skip newline
                #Complement
                    dseq[i] = _COMPLEMENTS[c]
                    j += 1
            #dseq[j] = 'Z'
            #jd = ''.join(dseq)
            dna  = dseq
            #print('REV')
        else:
            #print('NOREV')
            dd  = sequence[frame::]
            dna = dd
        #Translate dna codons to protein
        pseq = []
        for i in range(0, len(dna), 3):
            pseq.append(gencode[dna[i:i+3]])
        d = ''.join(pseq)
        protein = d
        #protein = ''.join([gencode[dna[i:i+3]] for i in xrange(0, len(dna), 3)])
        #Faster but less secure (wrt collissions) than stock uuid4
        id = uuid.UUID(int=random.getrandbits(128), version=4).hex
        seq = {
            'id': id,
            'dna': dna,
            'protein': protein,
            'name': name,
            'description': desc,
            'frame': frame+1
        }
        out = ''.join(['>', id, '\n', protein, '\n'])
        try:
            #result.append((id, json.dumps(seq), out))
            result.append((id, json.dumps(seq), out))
        except:
            print("* dd: ",dd)
            print("* dna: ", dna)
            print(seq)
            exit(1)
    return result
