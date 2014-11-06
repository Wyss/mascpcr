'''
======================================================================
 mascpcr: multiplex allele-specific colony (MASC) PCR design pipeline
======================================================================

``mascpcr`` is a native Python pipeline for designing multiplex allele-
specific colony (MASC) PCR primers for genome engineering applications. 

Setup / installation is fairly simple (the package may be used in place or
may be install in your Python site-packages directory by running this script).

Python dependencies:

    biopython       https://pypi.python.org/pypi/biopython
    bitarray        https://pypi.python.org/pypi/bitarray/
    libnano         https://github.com/Wyss/libnano
    mauve-py        https://github.com/Wyss/mauve-py
    numpy           https://pypi.python.org/pypi/numpy
    primer3-py      https://github.com/benpruitt/primer3-py
    six             https://pypi.python.org/pypi/six


See README.md for more information.

'''

from . import indexing, ioutil, offtarget, pipeline, primercandidate, \
			  primerpartition, genbankfeatures

from .pipeline import generateLUTs, findMascPrimers


__all__ = ['indexing', 'ioutil', 'offtarget', 'pipeline', 'primercandidate', 
		   'primerpartition', 'generateLUTs', 'findMascPrimers',
		   'genbankfeatures']