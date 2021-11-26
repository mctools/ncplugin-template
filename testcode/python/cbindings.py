__all__ = ['hooks']

import ctypes
import pathlib
import numpy as np
import random
import os

# Bindings which calls into the compiled C++ library. Specifically, the
# C-mangled functions in ../src/NCForPython.cc are called via ctypes. For
# convenience, the code here wraps a bit further, exposing the functionality via
# normal python functions with numpy support rather than ctypes semantics.

def _findlib(testcodelib = True):
    from . import _cmakeinfo as _info
    libname = _info.ncplugintestcode_libname if testcodelib else _info.ncplugin_libname
    libpath = pathlib.Path(__file__).parent / '../../../../lib/' / libname
    if not libpath.exists():
        raise RuntimeError('Could not find library %s'%libpath)
    return libpath.absolute().resolve()

def _setuptestlibhooks():
    from . import _cmakeinfo as _info
    hooks = {}
    #first load pluginlib by hand (in case environment is not set correctly):
    lib_ncplugin = ctypes.CDLL(_findlib(False))
    fct_getname = lib_ncplugin.ncplugin_getname
    fct_getname.restype = ctypes.c_char_p
    assert fct_getname().decode('ascii') == _info.ncplugin_name, "inconsistent library loaded"

    #Then the testcode library:
    lib = ctypes.CDLL(_findlib())

    #Dig out C-mangled (i.e. 'extern "C"') functions from the library, and
    #provide info about the function signature (i.e. return type and types and
    #numbers of arguments):

    #Use numpy ctypes integration, to transparently accept numpy arrays as double* arguments.
    from numpy.ctypeslib import ndpointer
    npdoubleptr = ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")

    f = lib.nctest_samplemanyscatmu
    f.restype = None # void
    #f.argtypes = [ ctypes.c_int, ctypes.c_double,  ctypes.c_double,ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_uint, npdoubleptr ]
    f.argtypes = [ ctypes.c_char_p, ctypes.c_double, ctypes.c_uint, npdoubleptr ]
    #f.argtypes = [ ctypes.c_int, ctypes.c_double,ctypes.c_double, ctypes.c_uint, npdoubleptr ]
    #f.argtypes = [ ctypes.c_int, ctypes.c_char_p,ctypes.c_double, ctypes.c_uint, npdoubleptr ]
    hooks['samplemanyscatmu'] = f

    f = lib.nctest_getmanyxsvalues
    f.restype = None # void
    #f.argtypes = [ ctypes.c_int,ctypes.c_double, ctypes.c_double,  ctypes.c_double,ctypes.c_double, ctypes.c_double, ctypes.c_uint, npdoubleptr, npdoubleptr ]
    f.argtypes = [ ctypes.c_char_p, ctypes.c_uint, npdoubleptr, npdoubleptr ]
    #f.argtypes = [ ctypes.c_int, ctypes.c_double,ctypes.c_uint, npdoubleptr, npdoubleptr ]
    #f.argtypes = [ ctypes.c_int, ctypes.c_char_p,ctypes.c_uint, npdoubleptr, npdoubleptr ]
    hooks['getmanyxsvalues'] = f

    hooks['ncplugin_register'] = lib.ncplugin_register

    return hooks

hooks = _setuptestlibhooks()

