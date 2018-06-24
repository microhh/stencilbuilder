#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)
wt = Field("wt", wloc)

dxidxi = Scalar("dxidxi")
dyidyi = Scalar("dxidyi")

dzi  = Vector("dzi" , zloc )
dzhi = Vector("dzhi", zhloc)

utrhs = grad2x( grad2x(u) ) * dxidxi \
      + grad2y( grad2y(u) ) * dyidyi \
      + grad2z( grad2z(u) * dzhi ) * dzi

wtrhs = grad2x( grad2x(w) ) * dxidxi \
      + grad2y( grad2y(w) ) * dyidyi \
      + grad2z( grad2z(w) * dzi ) * dzhi

printStencil(ut, utrhs, "=", "bot")
printEmptyLine(3)
printStencil(ut, utrhs, "=", "int")
printEmptyLine(3)
printStencil(ut, utrhs, "=", "top")

printEmptyLine(6)

printStencil(wt, wtrhs, "=", "bot")
printEmptyLine(3)
printStencil(wt, wtrhs, "=", "bot+1")
printEmptyLine(3)
printStencil(wt, wtrhs, "=", "int")
printEmptyLine(3)
printStencil(wt, wtrhs, "=", "top-1")
printEmptyLine(3)
printStencil(wt, wtrhs, "=", "top")
