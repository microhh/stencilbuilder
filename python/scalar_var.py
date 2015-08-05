#!/usr/bin/python

from StencilBuilder import *

b = Field("b", sloc)
w = Field("w", wloc)

bmean = Vector("bmean", zloc)

b2_shear = Field("b2_shear", sloc)
b2_turb  = Field("b2_turb" , sloc)

visc = Scalar("visc")
two  = Scalar("2.")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

rhs_shear = (b-bmean) * interpz(w) * gradz( interpz(bmean) ) * dzi4

rhs_turb = gradz ( interpz(b-bmean) * w ) * dzi4


printStencil(b2_shear, rhs_shear, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(b2_shear, rhs_shear, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(b2_shear, rhs_shear, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(b2_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(b2_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(b2_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)
