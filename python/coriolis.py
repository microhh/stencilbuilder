#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)
vt = Field("vt", vloc)

ug = Vector("ug", zloc)
vg = Vector("vg", zloc)

fc = Scalar("fc")

ugrid = Scalar("ugrid")
vgrid = Scalar("vgrid")

utrhs = fc * ( interpx( interpy( v ) ) + vgrid - vg )
vtrhs = fc * ( interpx( interpy( u ) ) + ugrid - ug )

printStencil(ut, utrhs, "+=", "int")
printStencil(vt, vtrhs, "+=", "int")

#printEmptyLine()
#printLoop(vt, vtrhs, "-=", "int")
