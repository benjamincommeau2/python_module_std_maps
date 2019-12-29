# -*- coding: utf-8 -*-
"""
Created on Fri May 22 00:04:38 2015

@author: linuxfreebird
"""
import sys
nerrors=1
with open('python_inputs.txt') as f:
  lines=f.read()
  try:
    c=lines.split(' ')[0]
    nerrors=int(c)
    if(nerrors<1):
      raise ValueError
  except ValueError:
    sys.exit("Error in "+sys.argv[0]+": First command line argument passed to the bash script calling this python script received invalid parameter nerrors="+c)
if(nerrors<1):
  nerrors=1
from distutils.core import setup, Extension
#import numpy.distutils.misc_util
 
c_ext = Extension("FGT", sources=["FGT.cpp"],
  language='c++',extra_compile_args=["-std=c++11","-fmax-errors="+str(nerrors)])
 
print("finished")
setup(
    ext_modules=[c_ext],
#    include_dirs=numpy.distutils.misc_util.get_numpy_include_dirs(),
)
