=======
mascpcr
=======

``mascpcr`` is a toolkit and configurable pipeline for generating multiplex 
allele-specific colony (MASC) PCR primer sets.

Features:
    * "Minimal input" pipeline -- requires 2 genbank files and target
      region indices to produce standard MASC PCR primer sets 
    * Module-based and command line interfaces 
    * Detailed reporting and order-ready .xlsx output 
    * Low level toolkit for more sophisticated MASC PCR-related design tasks

------

Installation
------------

Python library dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Make sure you have ``setuptools`` installed (``pip install setuptools``).
2. Install the necessary dependencies in the following order::

    cython
    biopython bitarray numpy openpyxl primer3-py six
    https://github.com/Wyss/libnano
    https://github.com/Wyss/mauve-py

3. (optional) Run `setup.py install` to install the package on your PYTHONPATH.

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

Copyright (C) 2014. Copyright (c) 2014 Wyss Institute at Harvard University
Ben Pruitt & Nick Conway
See LICENSE for full GPLv2 license.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
