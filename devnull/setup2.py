from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

setup(
        name = "libfermi node",
        ext_modules = cythonize('ferminode.py'), # accepts a glob pattern
        )

