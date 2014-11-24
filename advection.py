#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

dxi = Scalar("dxi")
dyi = Scalar("dyi")
dzi4 = Vector("dzi4", zloc)

ut = gradx( interpx(u) * interpx(u) ) * dxi \
   + grady( interpx(v) * interpy(u) ) * dyi \
   + gradz( interpx(w) * interpz(u) ) * dzi4

print("ut[i,j,k] = {0};\n".format(ut.getString(0,0,0,12)))
