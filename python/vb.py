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

bv_shear = Field("bv_shear", vloc)
bv_turb  = Field("bv_turb" , vloc)
bv_visc  = Field("bv_visc" , vloc)
bv_diss  = Field("bv_diss" , vloc)

rhs_shear = interpyz(w)*interpy(b-bmean)*gradz(interpz(vmean))*dzi4 \
          + interpyz(w)*(v-umean) * gradz(interpyz(bmean))*dzi4 \
          + N2*(interpxy(u-umean)*(v-vmean)*sinalpha + interpyz(w)*(v-vmean)*cosalpha )

rhs_turb = gradz( interpz(v-vmean) * interpy(w) * interpyz(b-bmean) ) * dzi4

rhs_visc = visc * gradz(gradz((v-umean) * interpy(b-bmean) ) * dzhi4) * dzi4

rhs_diss = 2.*visc * ( gradx( interpx(v-vmean) ) * dxi  * gradx( interpxy(b-bmean) )*dxi \
                     + grady( interpy(v-vmean) ) * dyi  * grady( (b-bmean) )*dyi \
                     + gradz( interpz(v-vmean) ) * dzi4 * gradz( interpyz(b-bmean) )*dzi4)

printStencil(bv_shear, rhs_shear, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(bv_turb, rhs_turb, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bv_turb, rhs_turb, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bv_turb, rhs_turb, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bv_turb, rhs_turb, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bv_turb, rhs_turb, "-=", "top", "[k]")

printEmptyLine(6)

printStencil(bv_visc, rhs_visc, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(bv_visc, rhs_visc, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bv_visc, rhs_visc, "+=", "bot+2", "[k]")
printEmptyLine(3)
printStencil(bv_visc, rhs_visc, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(bv_visc, rhs_visc, "+=", "top-2", "[k]")
printEmptyLine(3)
printStencil(bv_visc, rhs_visc, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bv_visc, rhs_visc, "+=", "top", "[k]")

printEmptyLine(6)

printStencil(bv_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bv_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bv_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bv_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bv_diss, rhs_diss, "-=", "top", "[k]")