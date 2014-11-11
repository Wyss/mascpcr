# Copyright (C) 2014. Ben Pruitt & Nick Conway
# See LICENSE for full GPLv2 license.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

"""
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

"""

from . import indexing, ioutil, offtarget, pipeline, primercandidate, \
			  primerpartition, genbankfeatures

from .pipeline import generateLUTs, findMascPrimers


__all__ = ['indexing', 'ioutil', 'offtarget', 'pipeline', 'primercandidate', 
		   'primerpartition', 'generateLUTs', 'findMascPrimers',
		   'genbankfeatures']