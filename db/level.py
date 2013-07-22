import leveldb

def open(path, truncate=False, sync=False):
    return (LDB(path+'.dna', truncate, sync), LDB(path+'.prot', truncate, sync))

class LDB:
    def __init__(self, path, truncate=False, sync=True):
        self.batchsize = 1000
        self.batchcount = 0
        self.sync = sync
        self.batch = leveldb.WriteBatch()
        self.db = leveldb.LevelDB(path, write_buffer_size=1024*(2 << 19))

    def get(self, key):
        return self.db.Get(key)

    def put(self, key, val):
        self.batch.Put(key, val)
        self.batchcount += 1
        if self.batchcount >= self.batchsize:
            self.flush()

    def __iter__(self):
        return self.db.RangeIter()

    def flush(self):
        self.db.Write(self.batch, sync=self.sync)
        self.batch = leveldb.WriteBatch()

    def close(self):
        self.flush()
        del self.db

    def truncate(self):
        for x in self.db.RangeIter(include_value=False):
            self.db.Delete(x)

