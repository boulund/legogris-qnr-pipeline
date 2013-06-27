#!/usr/bin/python
import json
import uuid
import bsddb3 as bsddb


def create_db(name):
    db = bsddb.hashopen(name, 'n')
    db.clear()
    return db

def open_fragments_passed(flag='c'):
    return bsddb.rnopen('fragments_passed.db', flag)

def open_dna_input(flag='c'):
    return open('dna-input.db', flag)

def open_fragments(flag='c'):
    return open('fragments.db', flag)

def open_clusters():
    db = bsddb.db.DB()
    db.set_flags(bsddb.db.DB_DUPSORT)
    db.open('clusters.db', bsddb.db.DB_HASH, bsddb.db.DB_CREATE)
    return db

def open(name, flag='c'):
    return bsddb.hashopen(name, flag)

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


