#!/usr/bin/python

from StencilBuilder import *
import numpy as np

b        = Field("b" , sloc )
u        = Field("v" , uloc )
v        = Field("v" , vloc )
w        = Field("w" , wloc )

umean    = Vector("umean", zloc)
vmean    = Vector("vmean", zloc)
bmean    = Vector("bmean", zloc)

visc     = Scalar("visc")
sinalpha = Scalar("sin(alpha)")
cosalpha = Scalar("cos(alpha)")
dxi      = Scalar("dxi")
dyi      = Scalar("dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

vw_shear = Field("uv_shear", uvloc)
vw_turb  = Field("uv_turb" , uvloc)
vw_buoy  = Field("uv_buoy" , uvloc)
vw_visc  = Field("uv_visc" , uvloc)
vw_diss  = Field("uv_diss" , uvloc)

rhs_shear = interpxyz(w)*interpx(v-vmean) * gradz(interpy(umean))*dzi4 \
          + interpxyz(w)*interpy(u-umean) * gradz(interpx(vmean))*dzi4

rhs_turb = gradz( interpxy(w) * interpyz(u-umean) * interpxz(v-vmean)) * dzi4

rhs_buoy = interpx(v-vmean)* interpxy(b-bmean) * sinalpha

rhs_visc = visc * gradz(gradz(interpy(u-umean)*interpx(v-vmean) ) * dzhi4) * dzi4

rhs_diss = 2.*visc * ( gradx( interpxy(u-umean) ) * dxi  * gradx(         (v-vmean) ) * dxi \
                     + grady(         (u-umean) ) * dyi  * grady( interpxy(v-vmean) ) * dyi \
                     + gradz( interpyz(u-umean) ) * dzi4 * gradz( interpxz(v-vmean) ) * dzi4)

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