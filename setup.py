from distutils.core import setup

from numpy.distutils.misc_util import get_numpy_include_dirs

setup(
    name =                  'mascpcr',
    include_dirs =          get_numpy_include_dirs(),
)
