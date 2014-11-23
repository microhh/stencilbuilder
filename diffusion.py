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

dxidxi = Scalar("dxidxi")
dyidyi = Scalar("dxidyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

ut = gradx( gradx(u) ) * dxidxi \
   + grady( grady(u) ) * dyidyi \
   + gradz( gradz(u) * dzhi4 ) * dzi4

print("ut[i,j,k] = {0};\n".format(ut.getString(0,0,0,12)))
