#!/bin/python

import numpy as np
from StencilBuilder import *

uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
sloc = np.array([0,0,0])
zloc  = 0
zhloc = 1

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)
evisc = Field("evisc", sloc)

dxi = Scalar("dxi")
dyi = Scalar("dyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

eviscy = interpy( interpx(evisc) )
eviscz = interpz( interpx(evisc) )

ut = gradx( evisc  * ( gradx(u)*dxi   + gradx(u)*dxi ) ) \
   + grady( eviscy * ( grady(u)*dyi   + gradx(v)*dxi ) ) \
   + gradz( eviscz * ( gradz(u)*dzhi4 + gradx(w)*dxi ) )

print("ut[i] = {0};\n".format(ut.getString(0,0,0,8)))
