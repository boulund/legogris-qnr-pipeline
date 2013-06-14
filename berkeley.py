#!/usr/bin/python
import json
import uuid
import bsddb


def create_db(name):
    return bsddb.hashopen(name, 'n')

def open_fragments_passed():
    return bsddb.rnopen('fragments_passed.db', 'c')

def open_fragments():
    return open('fragments.db')

def open(name):
    return bsddb.hashopen(name, 'c')

def test_insert():
    d = open()
    fragment = {'id': uuid.uuid4().hex,
                'name': 'test',
                'frame': 1,
                'dna': 'TCGAGCT',
                'aa': 'PCX',
                'passed_hmm': False,
                'hmm_score': 0.9
                }
    d[fragment['id']] = json.dumps(fragment)
    d.close()
    return fragment['id']

def test_get(id):
    d = open()
    fragment = json.loads(d[id])
    print fragment
    d.close()

def test():
    d = create_db()
    d.close()
    id = test_insert()
    test_get(id)


