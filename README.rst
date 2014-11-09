=======
mascpcr
=======
Ben Pruitt; Nick Conway
2014-11-07

``mascpcr`` is a toolkit and configurable pipeline for generating multiplex 
allele-specific colony (MASC) PCR primer sets.

Features:
    * "Minimal input" pipeline -- only requires 2 genbank files and target
      region indices to produce standard MASC PCR primer sets 
    * Module-based and command line interface 
    * Detailed reporting and order-ready .xlsx output 
    * Low level toolkit for more sophisticated MASC PCR-related design tasks

------

Installation
------------

1. Make sure you have ``setuptools`` installed (``pip install setuptools``).
2. Install the necessary dependencies:

pip-installable:
    - `biopython       <https://pypi.python.org/pypi/biopython>`_
    - `bitarray        <https://pypi.python.org/pypi/bitarray/>`_
    - `numpy           <https://pypi.python.org/pypi/numpy>`_
    - `openpyxl        <https://pypi.python.org/pypi/openpyxl>`_
    - `primer3-py      <https://github.com/benpruitt/primer3-py>`_
    - `six             <https://pypi.python.org/pypi/six>`_
non-pip installable (as of the writing of this document):
    - `libnano         <https://github.com/Wyss/libnano>`_
    - `mauve-py        <https://github.com/Wyss/mauve-py>`_

3. Run `setup.py install` to install the package on your PYTHONPATH.

*Use of a virtualenv is encouraged. For more configurable installation options, 
see the docs.*


Usage
-----

The simplest way to use ``mascpcr`` is via the command line interface::

    $ mascpcrcli recoded_genome.gb reference_genome.gb 1202000 1252000
         (1)            (2)                (3)           (4)     (5)
       
    1. CLI script installed in your PATH 
    2. Recoded / modified genome genbank file
    3. Reference genome genbank file
    4. MASC primer design start index
    5. MASC primer design end index (inclusive)


Of course, ``mascpcr`` can also be imported and used as a Python module. See 
``examples/`` for some example use cases, and ``docs/`` for more information.


Future plans
------------
``mascpcr`` functions well for our current needs, so we have no immediate plans
for further development. That being said, we are happy to help you expand / 
adapt the code base for your specific project needs.


License
-------
MIT
