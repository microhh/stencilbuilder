#!/usr/bin/python

from StencilBuilder import *
import numpy as np

p        = Field("p" , sloc )
b        = Field("b" , sloc )
u        = Field("u" , uloc )
w        = Field("w" , wloc )

umean    = Vector("umean", zloc)
pmean    = Vector("pmean", zloc)
bmean    = Vector("bmean", zloc)

visc     = Scalar("visc")
sinalpha = Scalar("sin(alpha)")
cosalpha = Scalar("cos(alpha)")
dxi      = Scalar("dxi")
dyi      = Scalar("dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

uw_shear = Field("uw_shear", uwloc)
uw_turb  = Field("uw_turb" , uwloc)
uw_buoy  = Field("uw_buoy" , uwloc)
uw_pres  = Field("uw_pres" , uwloc)
uw_rdstr = Field("uw_rdstr", uwloc)
uw_visc  = Field("uw_visc" , uwloc)
uw_diss  = Field("uw_diss" , uwloc)

rhs_shear = interpx(w)**2 * gradz(umean)*dzhi4

rhs_turb = gradz( interpxz(w)**2 * (u-umean)) * dzhi4

rhs_buoy = interpx(w)*interpxz(b-bmean) * sinalpha + interpz(u-umean)* interpxz(b-bmean) * cosalpha

rhs_pres = gradz( interpx(p-pmean) * (u-umean) ) * dzhi4

rhs_rdstr = interpxz( p-pmean ) * gradz( u-umean ) * dzhi4

rhs_visc = visc * gradz(gradz(interpz(u-umean)*interpx(w) ) * dzi4) * dzhi4

rhs_diss = 2.*visc * ( gradx( interpxz(u-umean) ) * dxi   * gradx( w ) * dxi \
                     + grady( interpyz(u-umean) ) * dyi   * grady( interpxy(w) )*dyi \
                     + gradz(         (u-umean) ) * dzhi4 * gradz( interpxz(w) )*dzhi4)

printStencil(uw_shear, rhs_shear, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(uw_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(uw_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(uw_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(uw_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(uw_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(uw_buoy, rhs_buoy, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(uw_pres, rhs_pres, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(uw_rdstr, rhs_rdstr, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(uw_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(uw_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(uw_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(uw_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(uw_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(uw_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(uw_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(uw_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "top", "[k]")