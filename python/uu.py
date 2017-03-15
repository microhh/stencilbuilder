#!/usr/bin/python

from StencilBuilder import *
import numpy as np

b        = Field("b", sloc )
u        = Field("u", uloc )
w        = Field("w", wloc )

umean    = Vector("umean", zloc)
bmean    = Vector("bmean", zloc)

visc     = Scalar("visc")
sinalpha = Scalar("sin(alpha)")
dxi      = Scalar("cgi*dxi")
dyi      = Scalar("cgi*dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

u2_shear = Field("u2_shear", uloc)
u2_turb  = Field("u2_turb" , uloc)
u2_buoy  = Field("u2_buoy" , uloc)
u2_visc  = Field("u2_visc" , uloc)
u2_diss  = Field("u2_diss" , uloc)

rhs_shear = 2.*interpxz(w)*(u-umean) * gradz(interpz(umean))*dzi4

rhs_turb = gradz( interpx(w) * interpz(u-umean)**2) * dzi4

rhs_buoy = 2.*(u-umean)*interpx( b-bmean ) * sinalpha 

rhs_visc = visc * gradz(gradz((u-umean)**2 ) * dzhi4) * dzi4

rhs_diss = 2.*visc * ( gradx( interpx(u-umean) ) * dxi  * gradx( interpx(u-umean) )*dxi \
                     + grady( interpy(u-umean) ) * dyi  * grady( interpy(u-umean) )*dyi \
                     + gradz( interpz(u-umean) ) * dzi4 * gradz( interpz(u-umean) )*dzi4)

printStencil(u2_shear, rhs_shear, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(u2_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(u2_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(u2_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(u2_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(u2_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(u2_buoy, rhs_buoy, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(u2_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(u2_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(u2_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(u2_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(u2_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(u2_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(u2_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(u2_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(u2_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(u2_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(u2_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(u2_diss, rhs_diss, "-=", "top", "[k]")