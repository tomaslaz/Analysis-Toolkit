# -*- coding: utf-8 -*-

"""
Wrapper to utilities.c

"""
import os

from ctypes import CDLL, c_double, c_int, POINTER

from .numpy_utils import CPtrToDouble, CPtrToInt



################################################################################

# load lib
_lib = CDLL(os.path.join(os.path.dirname(__file__), "_utilities.so"))

################################################################################

# atomicSeparation2 prototype
_lib.atomicSeparation2.restype = c_double
_lib.atomicSeparation2.argtypes = [c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_double, c_int, c_int, c_int]

# atomicSeparation2
def atomicSeparation2(ax, ay, az, bx, by, bz, xdim, ydim, zdim, pbcx, pbcy, pbcz):
    """
    Calculate atomic separation squared.
    
    """
    return _lib.atomicSeparation2(ax, ay, az, bx, by, bz, xdim, ydim, zdim, pbcx, pbcy, pbcz)

################################################################################
