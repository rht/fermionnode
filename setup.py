from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
        name = "libfermi node",
        ext_modules = cythonize('libferminode.py'), # accepts a glob pattern
        )

