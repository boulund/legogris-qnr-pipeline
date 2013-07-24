Installation
============
SofT does not require installation - all that is required is the installation of its dependencies and that a Python pipeline is constructed. Sample pipelines for identification of QNR and beta-lactamase genes are included in the samples directory.

Getting the source code
-----------------------
The source code is hosted at Bitbucket__, where the most recent version can be downloaded.

__ https://bitbucket.org/Legogris/qnr-pipeline/

Dependencies
------------
* Python 2.6
* HMMER__ >= v3.0
* `String Graph Assembler`__ >= v0.10.6
* Database engine. Either one of:
    * leveldb__ >= 1.12 (Recommended)
    * `Kyoto Cabinet`__ with `Python bindings`__

For development, Cython__ >= 0.14.1 is used to build the pyx file(s).

__ http://hmmer.janelia.org/
__ https://github.com/jts/sga
__ https://code.google.com/p/leveldb/
__ http://fallabs.com/kyotocabinet/
__ http://fallabs.com/kyotocabinet/pythonlegacydoc/
__ http://cython.org/#download

