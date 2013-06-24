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
_GENCODE = {
    'TTT': ord('F'),'TTC': ord('F'),'TTY': ord('F'),
    'TTA': ord('L'),'TTG': ord('L'),'TTR': ord('L'),
    'CTT': ord('L'),'CTC': ord('L'),'CTA': ord('L'),'CTG': ord('L'),'CTM': ord('L'),'CTR': ord('L'),'CTW': ord('L'),'CTS': ord('L'),'CTY': ord('L'),'CTK': ord('L'),'CTV': ord('L'),'CTH': ord('L'),'CTD': ord('L'),'CTB': ord('L'),'CTX': ord('L'),'CTN': ord('L'),
    'ATT': ord('I'),'ATC': ord('I'),'ATA': ord('I'),'ATY': ord('I'),'ATW': ord('I'),'ATM': ord('I'),'ATH': ord('I'), 'ATG': ord('M'),
    'GTT': ord('V'),'GTC': ord('V'),'GTA': ord('V'),'GTG': ord('V'),'GTM': ord('V'),'GTR': ord('V'),'GTW': ord('V'),'GTS': ord('V'),'GTY': ord('V'),'GTK': ord('V'),'GTV': ord('V'),'GTH': ord('V'),'GTD': ord('V'),'GTB': ord('V'),'GTX': ord('V'),'GTN': ord('V'),
    'TCT': ord('S'),'TCC': ord('S'),'TCA': ord('S'),'TCG': ord('S'),'TCM': ord('S'),'TCR': ord('S'),'TCW': ord('S'),'TCS': ord('S'),'TCY': ord('S'),'TCK': ord('S'),'TCV': ord('S'),'TCH': ord('S'),'TCD': ord('S'),'TCB': ord('S'),'TCX': ord('S'),'TCN': ord('S'),
    'CCT': ord('P'),'CCC': ord('P'),'CCA': ord('P'),'CCG': ord('P'),'CCM': ord('P'),'CCR': ord('P'),'CCW': ord('P'),'CCS': ord('P'),'CCY': ord('P'),'CCK': ord('P'),'CCV': ord('P'),'CCH': ord('P'),'CCD': ord('P'),'CCB': ord('P'),'CCX': ord('P'),'CCN': ord('P'),
    'ACT': ord('T'),'ACC': ord('T'),'ACA': ord('T'),'ACG': ord('T'),'ACM': ord('T'),'ACR': ord('T'),'ACW': ord('T'),'ACS': ord('T'),'ACY': ord('T'),'ACK': ord('T'),'ACV': ord('T'),'ACH': ord('T'),'ACD': ord('T'),'ACB': ord('T'),'ACX': ord('T'),'ACN': ord('T'),
    'GCT': ord('A'),'GCC': ord('A'),'GCA': ord('A'),'GCG': ord('A'),'GCM': ord('A'),'GCR': ord('A'),'GCW': ord('A'),'GCS': ord('A'),'GCY': ord('A'),'GCK': ord('A'),'GCV': ord('A'),'GCH': ord('A'),'GCD': ord('A'),'GCB': ord('A'),'GCX': ord('A'),'GCN': ord('A'),
    'TAT': ord('Y'),'TAC': ord('Y'),'TAY': ord('Y'),
    'TAA': ord('*'),'TAG': ord('*'),'TAR': ord('*'),
    'CAT': ord('H'),'CAC': ord('H'),'CAY': ord('H'),
    'CAA': ord('Q'),'CAG': ord('Q'),'CAR': ord('Q'),
    'AAT': ord('N'),'AAC': ord('N'),'AAY': ord('N'),
    'AAA': ord('K'),'AAG': ord('K'),'AAR': ord('K'),
    'GAT': ord('D'),'GAC': ord('D'),'GAY': ord('D'),
    'GAA': ord('E'),'GAG': ord('E'),'GAR': ord('E'),
    'TGT': ord('C'),'TGC': ord('C'),'TGY': ord('C'), 'TGA': ord('*'), 'TGG': ord('W'),
    'CGT': ord('R'),'CGC': ord('R'),'CGA': ord('R'),'CGG': ord('R'),'CGM': ord('R'),'CGR': ord('R'),'CGW': ord('R'),'CGS': ord('R'),'CGY': ord('R'),'CGK': ord('R'),'CGV': ord('R'),'CGH': ord('R'),'CGD': ord('R'),'CGB': ord('R'),'CGX': ord('R'),'CGN': ord('R'),
    'AGT': ord('S'),'AGC': ord('S'),'AGY': ord('S'),
    'AGA': ord('R'),'AGG': ord('R'),'AGR': ord('R'),
    'GGT': ord('G'),'GGC': ord('G'),'GGA': ord('G'),'GGG': ord('G'),'GGM': ord('G'),'GGR': ord('G'),'GGW': ord('G'),'GGS': ord('G'),'GGY': ord('G'),'GGK': ord('G'),'GGV': ord('G'),'GGH': ord('G'),'GGD': ord('G'),'GGB': ord('G'),'GGX': ord('G'),'GGN': ord('G'),
    'CT': ord('L'), 'GT': ord('V'), 'TC': ord('S'), 'CC': ord('P'), 'AC': ord('T'), 'GC': ord('A'), 'CG': ord('R'), 'GG': ord('G')

}
#Translates the supplied DNA string in all 6 reading frames and stores the result in a FASTA format text file as well as in serialized JSON in a supplied key/value store.
def translate_sequence(char *name, char *desc, char *sequence):
    cdef int i, j, frame, l
    cdef char *dna
    cdef char *protein
    cdef char c
    cdef char* s = NULL
    l = len(sequence)
    cdef char* dseq = <char *>calloc(l + 1, sizeof(char))
    cdef char* pseq = <char *>calloc(l/3 + 2, sizeof(char))
    #Local variables = less overhead
    result = []
    gencode = _GENCODE
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
            #jd = ''.join(dseq)
            #print('REV')
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
        #pseq = []
        #TODO: Optimize this too
        #for i in range(0, len(dna), 3):
        #    pseq.append(gencode[dna[i:i+3]])
        #d = ''.join(pseq)
        #protein = d
        j = 0
        l = len(dna)
        for i in range(0, l-1, 3):
            if i == l-2:
                nuc = dna[i:i+2]
            else:
                nuc = dna[i:i+3]
            try:
                c = gencode[nuc]
                pseq[j] = c
            except KeyError:
                pseq[j] = 'X'
            j += 1
        protein = pseq
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
            print("* dna: ", dna)
            print(seq)
            exit(1)
    return result
