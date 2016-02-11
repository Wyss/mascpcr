#!/usr/bin/env python

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
    libnano         https://github.com/Wyss/libnano
    mauve-py        https://github.com/Wyss/mauve-py
    numpy           https://pypi.python.org/pypi/numpy
    openpyxl        https://pypi.python.org/pypi/openpyxl
    primer3-py      https://github.com/benpruitt/primer3-py
    six             https://pypi.python.org/pypi/six


See README.md for more information.

"""

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='mascpcr',
    version='0.0.1',
    license='GPLv2',
    author='Ben Pruitt, Nick Conway',
    author_email='benjamin.pruitt@wyss.harvard.edu',
    url='https://github.com/wyss/mascpcr',
    description='MASC PCR design pipeline',
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
    ],
    packages=['mascpcr'],
    install_requires=['numpy', 'biopython', 'openpyxl', 'six',
                      'primer3-py', 'mauve-py', 'libnano'],
    test_suite='tests',
    scripts=['scripts/mascpcrcli', 'scripts/mascpcrfeatidx']
)