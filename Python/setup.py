#==================================================================================================
#                    Copyright (C) 2016 Olivier Mallet - All Rights Reserved                      
#==================================================================================================

# run with:
# python setup.py build_ext --inplace

# before running python program you need to export the library path:
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
               extra_compile_args = ['-std=c++14','-Wall','-march=native','-DUSE_BLAZE','-DBLAZE_BLAS_INCLUDE_FILE <cblas.h>'],
               extra_link_args = ['-L../lib'],
               language='c++')



ext.cython_directives = {'boundscheck': False,'wraparound': False}
# turn off bounds-checking for entire function
# turn off negative index wrapping for entire function

setup(cmdclass = {'build_ext' : build_ext}, ext_modules = [ext])
