import leveldb

def open(path, truncate=False):
    return (LDB(path+'.dna'), LDB(path+'.prot'))

class LDB:
    def __init__(self, path, truncate=False, sync=False):
        self.batchsize = 1000
        self.batchcount = 0
        self.sync = sync
        self.batch = leveldb.WriteBatch()
        self.db = leveldb.LevelDB(path, write_buffer_size=1024*(2 << 19))

    def put(self, key, val):
        self.batch.Put(key, val)
        self.batchcount += 1
        if self.batchcount >= self.batchsize:
            self.flush()

    def flush(self):
        self.db.Write(self.batch, sync=self.sync)
        self.batch = leveldb.WriteBatch()

    def close(self):
        del self.db
