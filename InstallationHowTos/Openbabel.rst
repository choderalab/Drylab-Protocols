How to install Openbabel
************************

Installation route 1
--------------------

Download from Github_ and compile from source.

.. _Github: https://github.com/openbabel/openbabel

Example method
~~~~~~~~~~~~~~
Compling from source when using Anaconda/Miniconda python. Compiling to ensure that the correct python bindings are installed is a little tricky. Below is an example where the user has python located in ~/miniconda2.

1. Ensure ``cmake``, ``swig``, and ``lxml`` are installed. For instance::

        conda install cmake lxml swig

2. Inside the directory containing uncompiled openbabel::

    mkdir build
    cd build
    cmake ../ -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON -DCMAKE_INSTALL_PREFIX=~/miniconda2 -DPYTHON_INCLUDE_DIR=~/miniconda2/include/python2. DCMAKE_LIBRARY_PATH=~/miniconda2/lib -DSWIG_DIR=~/miniconda/share/swig/3.0.2/ -DSWIG_EXECUTABLE=~/miniconda2/bin/swig -DPYTHON_LIBRARY=~/miniconda2/lib/libpython2.7.so
    make install

You have to replace ``~/miniconda2`` with whatever directory Anaconda or miniconda is located.


Installation route 2
--------------------

Install via conda using the omnia channel::

    conda install -c omnia openbabel

Note
~~~~

- 7/25/2016: Python bindings don't work for the above installation, but command line tools work fine.