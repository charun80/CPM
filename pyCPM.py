# -*- coding: utf-8 -*-
"""
Created on Sun Sep 10 14:05:08 2017

@author: Matthias Höffken
"""

from pyCPMFlow import computeCPMFlow, writeMatches, readMatches

import numpy as np
import sys
import cv2


img1Fname  = sys.argv[1]
img2Fname  = sys.argv[2]
matchFname = sys.argv[3]
nSteps     = int( sys.argv[4] )


#img1 = cv2.imread( img1Fname, flags=cv2.IMREAD_GRAYSCALE )
img1 = cv2.imread( img1Fname )
#img2 = cv2.imread( img2Fname, flags=cv2.IMREAD_GRAYSCALE )
img2 = cv2.imread( img2Fname )
print img2.shape

matches = computeCPMFlow( img1, img2, nSteps )
writeMatches( matchFname, matches )
rmatches = readMatches( matchFname )
assert( np.all( (matches == rmatches).ravel() ) )

