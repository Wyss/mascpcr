from distutils.core import setup
from Cython.Build import cythonize

from numpy.distutils.misc_util import get_numpy_include_dirs

setup(
    name =                  'mascpcr',
    ext_modules =           cythonize('indexutils.pyx'),
    include_dirs =          get_numpy_include_dirs(),
)
