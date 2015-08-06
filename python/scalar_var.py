#!/usr/bin/python

from StencilBuilder import *

b = Field("b", sloc)
w = Field("w", wloc)

bmean = Vector("bmean", zloc)

b2_shear = Field("b2_shear", sloc)
b2_turb  = Field("b2_turb" , sloc)
b2_visc  = Field("b2_visc" , sloc)
b2_diss  = Field("b2_diss" , sloc)

visc = Scalar("visc")

dxi = Scalar("cgi*dxi")
dyi = Scalar("cgi*dyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

rhs_shear = 2.*(b-bmean) * interpz(w) * gradz( interpz(bmean) ) * dzi4

rhs_turb = gradz ( interpz(b-bmean)**2 * w ) * dzi4

rhs_visc = visc * gradz( gradz( b-bmean ) * dzhi4 ) * dzi4

rhs_diss = 2.*visc * ( ( gradx( interpx( b-bmean ) ) * dxi  )**2 \
                     + ( grady( interpy( b-bmean ) ) * dyi  )**2 \
                     + ( gradz( interpz( b-bmean ) ) * dzi4 )**2 )

#printStencil(b2_shear, rhs_shear, "-=", "bot", "[k]")
#printEmptyLine(3)
printStencil(b2_shear, rhs_shear, "-=", "int", "[k]")
#printEmptyLine(3)
#printStencil(b2_shear, rhs_shear, "-=", "top", "[k]")

printEmptyLine(6)

#printStencil(b2_turb, rhs_turb, "-=", "bot", "[k]")
#printEmptyLine(3)
printStencil(b2_turb, rhs_turb, "-=", "int", "[k]")
#printEmptyLine(3)
#printStencil(b2_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

#printStencil(b2_visc, rhs_visc, "+=", "bot", "[k]")
#printEmptyLine(3)
printStencil(b2_visc, rhs_visc, "+=", "int", "[k]")
#printEmptyLine(3)
#printStencil(b2_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

#printStencil(b2_diss, rhs_diss, "-=", "bot", "[k]")
#printEmptyLine(3)
printStencil(b2_diss, rhs_diss, "-=", "int", "[k]")
#printEmptyLine(3)
#printStencil(b2_diss, rhs_diss, "-=", "top", "[k]")
