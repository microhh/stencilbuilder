#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)

dxi = Scalar("dxi")
dyi = Scalar("dyi")
dzi4 = Vector("dzi4", zloc)

utrhs = gradx( interpx(u) * interpx(u) ) * dxi \
      + grady( interpx(v) * interpy(u) ) * dyi \
      + gradz( interpx(w) * interpz(u) ) * dzi4

printStencil(ut, utrhs, "=")
