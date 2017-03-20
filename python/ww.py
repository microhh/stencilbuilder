#!/usr/bin/python

from StencilBuilder import *
import numpy as np

b        = Field("b", sloc )
w        = Field("w", wloc )
p        = Field("p", sloc)

umean    = Vector("umean", zloc)
bmean    = Vector("bmean", zloc)
pmean    = Vector("pmean", zloc)

visc     = Scalar("visc")
cosalpha = Scalar("cos(alpha)")
dxi      = Scalar("cgi*dxi")
dyi      = Scalar("cgi*dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

w2_turb  = Field("w2_turb" , wloc)
w2_buoy  = Field("w2_buoy" , wloc)
w2_pres  = Field("w2_pres" , wloc)
w2_rdstr = Field("w2_rdstr", wloc)
w2_visc  = Field("w2_visc" , wloc)
w2_diss  = Field("w2_diss" , wloc)

rhs_turb = gradz( interpz(w)**3) * dzhi4

rhs_buoy = 2.* w * interpz( b-bmean ) * cosalpha 

rhs_pres = 2. * gradz( interpz(w) * (p-pmean) ) * dzhi4

rhs_rdstr = 2.*interpz( p-pmean ) * gradz(interpz(w)) * dzhi4

rhs_visc = visc * gradz(gradz(w**2) * dzi4) * dzhi4

rhs_diss = 2.*visc * ( gradx( interpx(w) ) * dxi   * gradx( interpx(w) )*dxi \
                     + grady( interpy(w) ) * dyi   * grady( interpy(w) )*dyi \
                     + gradz( interpz(w) ) * dzhi4 * gradz( interpz(w) )*dzhi4)

printStencil(w2_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(w2_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(w2_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(w2_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(w2_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(w2_buoy, rhs_buoy, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(w2_pres, rhs_pres, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(w2_rdstr, rhs_rdstr, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(w2_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(w2_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(w2_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(w2_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(w2_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(w2_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(w2_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(w2_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(w2_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(w2_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(w2_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(w2_diss, rhs_diss, "-=", "top", "[k]")