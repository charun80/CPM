# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 22:32:29 2017

@author: Matthias Höffken
"""

from __future__ import print_function
import os


__all__ = ["writeMatches", "computeCPMFlow", "readMatches" ]


__CPMLibName = 'libctypesCPM.so'
__CPMLibPath = os.path.dirname( os.path.abspath(__file__) )

_cCPMCall = None
_cCPMFreeCall = None


import numpy as np
import ctypes as ct

####################################################################################

class __CPMResult(ct.Structure):
    """ structure for 1-channel image """
    _fields_=[("m_outMatching_pf", np.ctypeslib.ndpointer( dtype=np.float32, ndim=2, flags=('C_CONTIGUOUS','ALIGNED') ) ), \
              ("m_nMatches_i",ct.c_int) ]
    
    
    def __del__(self):
        _cCPMFreeCall( ct.byref(self) )
        for base in self.__class__.__bases__:
            # Avoid problems with diamond inheritance.
            basekey = 'del_' + str(base)
            if not hasattr(self, basekey):
                setattr(self, basekey, 1)
            else:
                continue
            # Call this base class' destructor if it has one.
            if hasattr(base, "__del__"):
                base.__del__(self)
    
    
    def getArray( self ):
        if 0 < self.m_nMatches_i:
            # workaround for bug in numpy.ctypeslib.as_array
            l_arrayPointer = ct.cast( self.m_outMatching_pf, ct.POINTER( ct.c_float ) )
            return np.ctypeslib.as_array( l_arrayPointer, shape=( self.m_nMatches_i,4) ).copy()
        else:
            return np.zeros( (0,4), dtype=np.float32 )




def __loadCPMLibrary():
    cpmlib = np.ctypeslib.load_library( __CPMLibName, __CPMLibPath )
    
    # Caller Return type
    cpmlib.computeCPMFlow.restype = __CPMResult
    
    # Caller parameters
    cpmlib.computeCPMFlow.argtypes = [ \
        np.ctypeslib.ndpointer( dtype=np.float32, flags=('C_CONTIGUOUS','ALIGNED')), \
        np.ctypeslib.ndpointer( dtype=np.float32, flags=('C_CONTIGUOUS','ALIGNED')), \
        ct.c_int, ct.c_int, ct.c_int, \
        ct.c_int  ]
    
    # Caller freeer
    cpmlib.clearFlowResult.argtypes = [ ct.POINTER( __CPMResult ) ]
    
    global _cCPMCall
    global _cCPMFreeCall
    
    _cCPMCall = cpmlib.computeCPMFlow
    _cCPMFreeCall = cpmlib.clearFlowResult
    

####################################################################################

__loadCPMLibrary()

####################################################################################


def computeCPMFlow( img1, img2, n_steps=3 ):
    img1 = np.ascontiguousarray( img1, dtype=np.float32 )
    img2 = np.ascontiguousarray( img2, dtype=np.float32 )
    
    assert( np.all( img1.shape == img2.shape ) )
    
    if 2 != img1.ndim:
        if (3 != img1.ndim) or (3 != img1.shape[2]):
            raise ValueError("Wrong image shapes: %s" % str(img1.shape) )
        nChannels = img1.shape[2]
    else:
        nChannels = 1
    
    mval = max( img1.max(), img2.max() )
    if 1 < mval:
        # normalization
        #nval = 1. / float(mval)
        nval = 1. / 255.
        img1 *= nval
        img2 *= nval
    
    n_steps = max(1, int(n_steps))
    
    # Process images
    cflow_res = _cCPMCall( img1, img2, img1.shape[0], img1.shape[1], nChannels, n_steps )
    
    # parse output
    return cflow_res.getArray()


####################################################################################


def writeMatches( fFname_s, fMatches ):
    # fMatches: [ N x 4 ]
    N = fMatches.shape[0]
    
    lFOut = file( fFname_s, "w" )
    
    for i in xrange(N): 
        lFOut.write( "%.0f %.0f %.0f %.0f\n" % tuple(fMatches[i,:]) )
    
    lFOut.flush()
    lFOut.close()


def readMatches( fFname_s ):
    # fMatches: [ N x 4 ]
    
    lMatchFlines = file( fFname_s, "r" ).readlines()
    N = len( lMatchFlines )
    lMatches = np.empty( (N,4), np.float32 )
    lMatches[:] = float('nan')
    
    for i in xrange(N):
        lMatches[i,:] = map( float, lMatchFlines[i].split(" ") )
    
    return lMatches


####################################################################################

if __name__ == "__main__":
    pass
    

