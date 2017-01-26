#==================================================================================================
#                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

# run with:
# python setup.py build_ext --inplace

# before running Python program you will need to export the C++ shared library path:
# export LD_LIBRARY_PATH=/path/to/URT/lib:$LD_LIBRARY_PATH

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
from numpy import get_include

# linking to C++ libURT.so library
ext = Extension('CyURT',
               sources = ['CyURT.pyx'],
               include_dirs = [get_include()],
               libraries = ['URT'],
               extra_compile_args = ['-std=c++11','-Wall','-march=native','-DUSE_BLAZE','-DBLAZE_BLAS_INCLUDE_FILE <mkl_cblas.h>'],
               extra_link_args = ['-L../lib'],
               language='c++')

# NB: in extra_compile_args replace the header <mkl_cblas.h> by the one associated with the BLAS/LAPACK replacement library, for example OpenBLAS will use the header <cblas.h>

ext.cython_directives = {'boundscheck': False,'wraparound': False}
# turn off bounds-checking for entire function
# turn off negative index wrapping for entire function

setup(cmdclass = {'build_ext' : build_ext}, ext_modules = [ext])
