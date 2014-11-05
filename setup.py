
'''
======================================================================
 mascpcr: multiplex allele-specific colony (MASC) PCR design pipeline
======================================================================

``mascpcr`` is a native Python pipeline for designing multiplex allele-
specific colony (MASC) PCR primers for genome engineering pipelines. 

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

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(
    name='mascpcr',
    version='0.0.1',
    license='GPLv2',
    author='Ben Pruitt, Nick Conway',
    author_email='bpruittvt@gmail.com',
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
    test_suite= "tests"
)