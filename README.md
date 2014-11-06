MASC PCR design pipeline
========================

Mid-development stage pipeline for designing MASC PCR primers for recoded DNA sequences. 

### Installation
You can import and run the pipeline from within the package directory, assuming you have installed all of the necessary dependencies:

Python dependencies:

    biopython       https://pypi.python.org/pypi/biopython
    bitarray        https://pypi.python.org/pypi/bitarray/
    libnano         https://github.com/Wyss/libnano
    mauve-py        https://github.com/Wyss/mauve-py
    numpy           https://pypi.python.org/pypi/numpy
    openpyxl        https://pypi.python.org/pypi/openpyxl
    primer3-py      https://github.com/benpruitt/primer3-py
    six             https://pypi.python.org/pypi/six

If you want to be able to use the pipeline outside of the package directory or take advantage of the included scripts you can install the module by running the include setup.py script:

	$ python setup.py install

### Testing
The package includes a fairly comprehensive test suite, along with the respective input and output files (see `tests/`). The best way to run the tests is using the [`nosetests`](https://nose.readthedocs.org/en/latest/) script from the package root. Alternatively, if you have `setuptools` installed you can simply call:

	$ setup.py test 

### Documentation
We may provide more extensive documentation at some point in the future, but for now we will refer you to the `examples/`, `scripts/`, and `tests/` directories for usage examples and some handy command line scripts (these scripts will be installed in your PATH by the setup.py script so you can use them from within any directory).