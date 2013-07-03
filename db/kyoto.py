from kyotocabinet import DB
from util import PathError

def open(path, truncate=False):
    return (KDB(path+'.dna', truncate), KDB(path+'.prot', truncate))

class KDB:
    def __init__(self, path, truncate=False):
        self.db = DB()
        self.batchsize = 1000
        self.batch = {}
        dbparams = '.kct#apow=0#bnum=10000000#msiz='+str(2<<30)
        if truncate:
            result = self.db.open(path+dbparams, DB.OWRITER | DB.OCREATE | DB.OTRUNCATE)
        else:
            result = self.db.open(path+dbparams, DB.OWRITER)
        if not result:
            raise PathError('DNA outdb open error: %s ' % self.db.error())
            exit(1)

    def put(self, key, val):
        self.batch[key] = val
        if len(self.batch) >= self.batchsize:
            self.flush()

    def flush(self):
        self.db.set_bulk(self.batch)
        self.batch = {}

    def close(self):
        del self.db
