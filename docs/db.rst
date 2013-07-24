db
==
In SofT, databases are generally instantiated in a pair using the modules' `open` method, where the first database is used for raw DNA sequences and the second one stores JSON serialiazed dictionaries containing the translated protein sequences. Numerical, incremental IDs should be used for the DNA database. The IDs of entries in the protein database can be derived from the corresponding DNA sequence, see :py:class:`sieve.dnareader`.

.. automodule:: db
    :members:
    :undoc-members:

    level
    =======
    .. automodule:: db.level
        :members:
        :undoc-members:

    kyoto
    =======
    .. automodule:: db.kyoto
        :members:
        :undoc-members:

