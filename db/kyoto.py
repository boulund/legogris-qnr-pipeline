from kyotocabinet import DB
from util import PathError

def open(path, truncate=False):
    """
    Open a new connection to a pair of DNA and protein databases using the Kyoto Cabinet engine.

    Args:
        * path (str): Base path and prefix for database directories.

    Kwargs:
        * truncate (bool, False): If databases should be truncated before opening.

    Returns:
        dbs (tuple): Two LDB instances.
    """
    return (KDB(path+'.dna', truncate), KDB(path+'.prot', truncate))

#TODO: __iter__ and trunc
class KDB:
    def __init__(self, path, truncate=False):
        """
        Open a new connection to a database using the Kyoto Cabinet engine.

        Args:
            * path (str): Path to database.

        Kwargs:
            * truncate (bool, False): If database should be truncated before opening.
        """
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
        self.batch[key] = val
        if len(self.batch) >= self.batchsize:
            self.flush()

    def flush(self):
        """
        Write `put` calls to database.
        """
        self.db.set_bulk(self.batch)
        self.batch = {}

    def close(self):
        """
        Flush the database and delete the connection to it.
        """
        self.flush()
        del self.db
