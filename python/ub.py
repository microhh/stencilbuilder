#!/usr/bin/python

from StencilBuilder import *
import numpy as np

b        = Field("b", sloc )
u        = Field("u", uloc )
v        = Field("v", vloc )
w        = Field("w", wloc )
p        = Field("p", sloc)

umean    = Vector("umean", zloc)
vmean    = Vector("vmean", zloc)
bmean    = Vector("bmean", zloc)
pmean    = Vector("pmean", zloc)

visc     = Scalar("visc")
sinalpha = Scalar("sin(alpha)")
cosalpha = Scalar("cos(alpha)")
N2       = Scalar("N2")
dxi      = Scalar("cgi*dxi")
dyi      = Scalar("cgi*dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

bu_shear = Field("bu_shear", uloc)
bu_turb  = Field("bu_turb" , uloc)
bu_buoy  = Field("bu_buoy" , uloc)
bu_visc  = Field("bu_visc" , uloc)
bu_diss  = Field("bu_diss" , uloc)

rhs_shear = interpxz(w)*interpx(b-bmean)*gradz(interpz(umean))*dzi4 \
          + (u-umean)*interpxz(w) * gradz(interpxz(bmean))*dzi4 \
          + N2*((u-umean)**2 * sinalpha + (u-umean)*interpxz(w)*cosalpha)

rhs_turb = gradz( interpx(w) * interpz(u-umean) * interpxz(b-bmean) ) * dzi4

rhs_buoy = interpx( b-bmean )**2 * sinalpha 

rhs_visc = visc * gradz(gradz((u-umean) * interpx(b-bmean) ) * dzhi4) * dzi4

rhs_diss = 2.*visc * ( gradx( interpx(u-umean) ) * dxi  * gradx( (b-bmean) )*dxi \
                     + grady( interpy(u-umean) ) * dyi  * grady( interpxy(b-bmean) )*dyi \
                     + gradz( interpz(u-umean) ) * dzi4 * gradz( interpxz(b-bmean) )*dzi4)

printStencil(bu_shear, rhs_shear, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(bu_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bu_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bu_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bu_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bu_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(bu_buoy, rhs_buoy, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(bu_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(bu_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bu_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(bu_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(bu_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(bu_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bu_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(bu_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bu_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bu_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bu_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bu_diss, rhs_diss, "-=", "top", "[k]")