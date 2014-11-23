#!/bin/python

import numpy as np
from StencilBuilder import *

uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])

u  = Field("u" , uloc)
v  = Field("v" , vloc)
w  = Field("w" , wloc)

fc = Scalar("fc")

ut = fc * interpx( interpy( v ) )
vt = fc * interpx( interpy( u ) )

print("ut[i,j,k] = {0};\n".format(ut.getString(0,0,0,12)))
print("vt[i,j,k] = {0};\n".format(vt.getString(0,0,0,12)))
