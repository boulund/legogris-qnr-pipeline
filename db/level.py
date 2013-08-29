import plyvel

def open(path, truncate=False, sync=False):
    """
    Open a new connection to a pair of DNA and protein databases using the Google leveldb engine.

    Args:
        * path (str): Base path and prefix for database directories.

    Kwargs:
        * truncate (bool, False): If databases should be truncated before opening.
        * sync (bool, False): If writes to the databases should be performed synchronously.

    Returns:
        dbs (tuple): Two LDB instances.
    """
    return (LDB(path+'.dna', truncate, sync), LDB(path+'.prot', truncate, sync))

class LDB:
    def __init__(self, path, truncate=False, sync=True):
        """
        Open a new connection to a database using the Google leveldb engine.

        Args:
            * path (str): Path to database.

        Kwargs:
            * truncate (bool, False): If database should be truncated before opening.
            * sync (bool, False): If writes to the database should be performed synchronously.
        """

        self.batchsize = 1000
        self.batchcount = 0
        self.sync = sync
        self.db = plyvel.DB(path, create_if_missing=True, write_buffer_size=1024*(2 << 19))
        self.batch = self.db.write_batch(sync=sync)

    def get(self, key):
        """
        Retrieve the item with the given `key`.
        """
        return self.db.get(key)

    def put(self, key, val):
        """
        Put `val` at `key`.
        Note that disk writing is done in batches, so be sure to call `close` or `flush` to make sure that values are put into the store.
        """
        self.batch.put(key, val)
        self.batchcount += 1
        if self.batchcount >= self.batchsize:
            self.flush()

    def __iter__(self):
        return self.db.iterator()

    def flush(self):
        """
        Write `put` calls to database.
        """
        self.batch.write()
        self.batch = self.db.write_batch()

    def close(self):
        """
        Flush the database and delete the connection to it.
        """
        self.flush()
        del self.db

    def truncate(self):
        """
        Delete all values in the database. This is done using iteration, which is obviously not very effective, so should be done sparingly.
        """
        for key in self.db.iterator(include_value=False):
            self.db.delete(key)

