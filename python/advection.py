#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)
wt = Field("wt", wloc)

dxi = Scalar("dxi")
dyi = Scalar("dyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

utrhs = gradx( interpx(u) * interpx(u) ) * dxi \
      + grady( interpx(v) * interpy(u) ) * dyi \
      + gradz( interpx(w) * interpz(u) ) * dzi4

wtrhs = gradx( interpz(u) * interpx(w) ) * dxi \
      + grady( interpz(v) * interpy(w) ) * dyi \
      + gradz( interpz(w) * interpz(w) ) * dzhi4

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

#printLoop(ut, utrhs, "=", kstart="grid->kstart+1", kend="grid->kend-1"  )
