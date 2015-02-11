#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)
vt = Field("vt", vloc)

fc = Scalar("fc")

utrhs = fc * interpx( interpy( v ) )
vtrhs = fc * interpx( interpy( u ) )

printStencil(ut, utrhs, "=")
printStencil(vt, vtrhs, "=")
