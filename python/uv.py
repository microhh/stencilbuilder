#!/usr/bin/python

from StencilBuilder import *
import numpy as np

p        = Field("p" , sloc )
b        = Field("b" , sloc )
v        = Field("v" , vloc )
w        = Field("w" , wloc )

vmean    = Vector("vmean", zloc)
pmean    = Vector("pmean", zloc)
bmean    = Vector("bmean", zloc)

visc     = Scalar("visc")
sinalpha = Scalar("sin(alpha)")
cosalpha = Scalar("cos(alpha)")
dxi      = Scalar("cgi*dxi")
dyi      = Scalar("cgi*dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

vw_shear = Field("vw_shear", vwloc)
vw_turb  = Field("vw_turb" , vwloc)
vw_buoy  = Field("vw_buoy" , vwloc)
vw_pres  = Field("vw_pres" , vwloc)
vw_rdstr = Field("vw_rdstr", vwloc)
vw_visc  = Field("vw_visc" , vwloc)
vw_diss  = Field("vw_diss" , vwloc)

rhs_shear = interpy(w)**2 * gradz(vmean)*dzhi4

rhs_turb = gradz( interpyz(w)**2 * (v-vmean)) * dzhi4

rhs_buoy = interpz(v-vmean)* interpyz(b-bmean) * cosalpha

rhs_pres = gradz( interpy(p-pmean) * (v-vmean) ) * dzhi4

rhs_rdstr = interpyz( p-pmean ) * gradz( v-vmean ) * dzhi4

rhs_visc = visc * gradz(gradz(interpz(v-vmean)*interpy(w) ) * dzi4) * dzhi4

rhs_diss = 2.*visc * ( gradx( interpxz(v-vmean) ) * dxi   * gradx( interpxy(w) ) * dxi \
                     + grady( interpyz(v-vmean) ) * dyi   * grady( w )*dyi \
                     + gradz(         (v-vmean) ) * dzhi4 * gradz( interpyz(w) )*dzhi4)

printStencil(vw_shear, rhs_shear, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(vw_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(vw_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(vw_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(vw_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(vw_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(vw_buoy, rhs_buoy, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(vw_pres, rhs_pres, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(vw_rdstr, rhs_rdstr, "+=", "int", "[k]")

printEmptyLine(6)

printStencil(vw_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(vw_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(vw_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(vw_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(vw_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(vw_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(vw_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(vw_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(vw_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(vw_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(vw_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(vw_diss, rhs_diss, "-=", "top", "[k]")