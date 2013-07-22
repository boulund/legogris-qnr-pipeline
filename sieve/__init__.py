#!/bin/python
from __future__ import print_function

import time

class Sieve(object):
    """
    Abstract base-class for gene-filtering sieve implementations.
    """
    def __init__(self, params, logfile, name, param_names):
        """
        Instantiate a new Sieve. Should not be called directly, but extended from subclass.
        Generally, name and params_names are specified in the extending constructors while params and logfile are passed on from initialization.

        Args:
            params (dict): Supplied parameter to this instance of the sieve. Parameters are defined in param_names.

            logfile (`Logger`): Written to for logging purposes.

            name (str): Name of the Sieve.

            param_names (list): List of parameters. Items in list can be either:
                * String, meaning the parameter is mandatory.
                * A tuple of a string and default value, would the value not be supplied in the params argument.

        Raises: KeyError
        """

        self.param_names = param_names
        self.name = name
        self.logfile = logfile
        for param in self.param_names:
            if isinstance(param, str):
                if not param in params:
                    raise KeyError('Missing mandatory parameter %s for sieve %s' % (param, self.name))
                param_name = param
            else:
                param_name = param[0]
                default_value = param[1]
            if param_name in params:
                setattr(self, param_name, params[param_name])
            else:
                setattr(self, param_name, default_value)

    def run(self, indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath):
        """
        Abstract.

        Args:
            indnadb, inprotdb (`DB`): Key-value stores for looking up DNA strings and protein sequences, respectively.

            infilepath (str): Path to file with input sequences to handle. Implementing subclasses should handle the following formats, determined from file ending:
                * .fa(sta): FASTA
                * .f(ast)q: FASTQ
                * .fa(sta).gz: GZip-compressed FASTA
                * .f(ast)q.gz: GZip-compressed FASTQ

            outdnadb, outprotdb (`DB`): Key-value-stores for storing output DNA strings and protein sequences after processing.

            outfilepath (str): Path to file to store output sequences. Format should be determined in a similar manner as for infilepath. Additionally, the existence of the following substrings in the filename determine output:
                * .n : Nucleotides (DNA)
                * .p : Amino acids (Protein)
        """

        raise NotImplementedError('Please implement this method.')

def run_sieve(sieve, paths, logfile, dbengine):
    """
    Run the supplied Sieve instance with the files and databases at the supplied paths.

    Args:
        sieve (`Sieve`): Sieve subclass instance to run.

        paths (tuple): Tuple of strings specifying the files used for (in order): (indb, infile, outdb, outfile).

        logfile (`Logger`): Log file object

        dbengine (module): Database engine to use for indb and outdb.

    """
    (indbpath, infilepath, outdbpath, outfilepath) = paths
    indnadb = inprotdb =  outdnadb =  outprotdb = None
    try:
        if indbpath:
            (indnadb, inprotdb) = dbengine.open(indbpath, truncate=False)
        if outdbpath:
            (outdnadb, outprotdb) = dbengine.open(outdbpath, truncate=True)
        logfile.writeline('Start: %s at %s' % (sieve.name, time.asctime(time.localtime())))
        return sieve.run(indnadb, inprotdb, infilepath, outdnadb, outprotdb, outfilepath)
    finally:
        if not indnadb is None:
            indnadb.close()
        if not inprotdb is None:
            inprotdb.close()
        if not outdnadb is None:
            outdnadb.close()
        if not outprotdb is None:
            outprotdb.close()
        logfile.writeline('Finish: %s at %s' % (sieve.name, time.asctime(time.localtime())))
        logfile.flush()

def _run_sieves(sieves, dbs, files, logfile, dbengine, startindex=0, endindex=-1):
    if endindex == -1:
        endindex = len(sieves)
    for i in xrange(startindex, endindex):
        if isinstance(sieves[i], tuple):
            (s, params) = sieves[i]
            sieve = s.sieve(params, logfile)
            run_sieve(sieve, (dbs[i], files[i], dbs[i+1], files[i+1]), logfile, dbengine)
        elif isinstance(sieves[i], list):
            for j in xrange(len(sieves[i])):
                s = sieves[i][j]
                sfiles = [files[i]] + [f+str(j) for f in files[i+1::]]
                sdbs = [dbs[i]] + [db+str(j) for db in dbs[i+1::]]
                ssieves = [s]+sieves[i+1::]
                _run_sieves(ssieves, sdbs, sfiles, logfile, dbengine)
            return

def run_sieves(sieves, dbs, files, logfile, dbengine, startindex=0, endindex=-1):
    """
    Run a pipeline of Sieves, using the output of one as the input for the next.

    Args:
        sieves (list): List of `n` tuples of sieve modules and parameters: [(sieve_module1, {'param_one': 1}), ...]

        dbs (list): List of `n`\+1 strings to use as paths for databases.

        files (list): List of `n`\+1 strings to use as paths for input/output files.

        dbengine (module): Database engine to use for databases.

    Kwargs:
        startindex (int): To run a subset of the defined sieves, specify startindex > 0 to use as starting index.

        endindex (int): To run a subset of the defined sieves, specify n >= endindex > startindex to use as exclusive end index.

    """
    _run_sieves(sieves, dbs, files, logfile, dbengine, startindex, endindex)
