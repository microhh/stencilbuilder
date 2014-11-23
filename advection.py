#!/bin/python

import numpy as np
from StencilBuilder import *

uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
zloc  = 0
zhloc = 1

u  = Field("u", uloc)
v  = Field("v", vloc)
w  = Field("w", wloc)

dzi4 = Vector("dzi4", zloc)

ut = gradx( interpx(u) * interpx(u) ) \
   + grady( interpx(v) * interpy(u) ) \
   + gradz( interpx(w) * interpz(u) ) * dzi4

print("ut[i] = {0};\n".format(ut.getString(0,0,0,8)))
