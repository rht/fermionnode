from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
        name = "pyqm",
        ext_modules = cythonize('pyqm.py'), # accepts a glob pattern
        )

