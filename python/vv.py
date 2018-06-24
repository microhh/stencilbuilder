#!/usr/bin/python

from StencilBuilder import *
import numpy as np

b        = Field("b", sloc )
v        = Field("v", vloc )
w        = Field("w", wloc )

vmean    = Vector("vmean", zloc)
bmean    = Vector("bmean", zloc)

visc     = Scalar("visc")
sinalpha = Scalar("sin(alpha)")
dxi      = Scalar("dxi")
dyi      = Scalar("dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

v2_shear = Field("v2_shear", vloc)
v2_turb  = Field("v2_turb" , vloc)
v2_visc  = Field("v2_visc" , vloc)
v2_diss  = Field("v2_diss" , vloc)

rhs_shear = 2.*interpyz(w)*(v-vmean) * gradz(interpz(vmean))*dzi4

rhs_turb = gradz( interpy(w) * interpz(v-vmean)**2) * dzi4

rhs_visc = visc * gradz(gradz((v-vmean)**2 ) * dzhi4) * dzi4

rhs_diss = 2.*visc * ( gradx( interpx(v-vmean) ) * dxi  * gradx( interpx(v-vmean) )*dxi \
                     + grady( interpy(v-vmean) ) * dyi  * grady( interpy(v-vmean) )*dyi \
                     + gradz( interpz(v-vmean) ) * dzi4 * gradz( interpz(v-vmean) )*dzi4)

printStencil(v2_shear, rhs_shear, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(v2_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(v2_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(v2_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(v2_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(v2_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(v2_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(v2_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(v2_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(v2_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(v2_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(v2_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(v2_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(v2_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(v2_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(v2_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(v2_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(v2_diss, rhs_diss, "-=", "top", "[k]")