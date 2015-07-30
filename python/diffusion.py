#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)
wt = Field("wt", wloc)

dxidxi = Scalar("dxidxi")
dyidyi = Scalar("dxidyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

utrhs = gradx( gradx(u) ) * dxidxi \
      + grady( grady(u) ) * dyidyi \
      + gradz( gradz(u) * dzhi4 ) * dzi4

wtrhs = gradx( gradx(w) ) * dxidxi \
      + grady( grady(w) ) * dyidyi \
      + gradz( gradz(w) * dzi4 ) * dzhi4

printStencil(ut, utrhs, "=", "bot")
printStencil(ut, utrhs, "=", "int")
printStencil(ut, utrhs, "=", "top")

printStencil(wt, wtrhs, "=", "bot")
printStencil(wt, wtrhs, "=", "int")
printStencil(wt, wtrhs, "=", "top")
