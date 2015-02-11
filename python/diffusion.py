#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)

dxidxi = Scalar("dxidxi")
dyidyi = Scalar("dxidyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

utrhs = gradx( gradx(u) ) * dxidxi \
      + grady( grady(u) ) * dyidyi \
      + gradz( gradz(u) * dzhi4 ) * dzi4

printStencil(ut, utrhs, "=")
