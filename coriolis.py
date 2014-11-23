#!/bin/python

import numpy as np
from StencilBuilder import *

uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])

u  = Scalar("u" , uloc)
v  = Scalar("v" , vloc)
w  = Scalar("w" , wloc)

ut = interpx( interpy( v ) )
vt = interpx( interpy( u ) )

print("ut[i] = {0};\n".format(ut.getString(0,0,0,8)))

print("vt[i] = {0};\n".format(vt.getString(0,0,0,8)))
