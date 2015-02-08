#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

dxidxi = Scalar("dxidxi")
dyidyi = Scalar("dxidyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

ut = gradx( gradx(u) ) * dxidxi \
   + grady( grady(u) ) * dyidyi \
   + gradz( gradz(u) * dzhi4 ) * dzi4

print("ut[i,j,k] = {0};\n".format(ut.getString(0,0,0,12)))
