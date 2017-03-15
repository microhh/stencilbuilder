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
cosalpha = Scalar("cos(alpha)")
N2       = Scalar("N2")
dxi      = Scalar("cgi*dxi")
dyi      = Scalar("cgi*dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

bw_shear = Field("bw_shear", sloc)
bw_turb  = Field("bw_turb" , sloc)
bw_visc  = Field("bw_visc" , sloc)
bw_diss  = Field("bw_diss" , sloc)

rhs_shear = 2.*interpz(w)*(b-bmean) * gradz(interpz(bmean))*dzi4 \
          + 2.*N2*(interpx(u-umean)*(b-bmean) * sinalpha + interpz(w)*(b-bmean) * cosalpha ) 

rhs_turb = gradz( w * interpz(b-bmean)**2 ) * dzi4

rhs_visc = visc * gradz(gradz((b-bmean)**2 ) * dzhi4 ) * dzi4

rhs_diss = 2.*visc * ( gradx( interpx(b-bmean) ) * dxi  * gradx( interpx(b-bmean) ) * dxi \
                     + grady( interpy(b-bmean) ) * dyi  * grady( interpy(b-bmean) ) * dyi \
                     + gradz( interpz(b-bmean) ) * dzi4 * gradz( interpz(b-bmean) ) * dzi4 )

printStencil(bw_shear, rhs_shear, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(bw_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bw_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bw_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bw_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bw_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(bw_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(bw_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bw_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(bw_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(bw_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(bw_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bw_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(bw_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "top", "[k]")