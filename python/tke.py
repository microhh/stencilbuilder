#!/usr/bin/python

from StencilBuilder import *
import numpy as np

p        = Field("p" , sloc )
b        = Field("b" , sloc )
u        = Field("u" , uloc )
v        = Field("u" , vloc )
w        = Field("w" , wloc )

umean    = Vector("umean", zloc)
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

uw_shear = Field("tke_shear", sloc)
uw_turb  = Field("tke_turb" , sloc)
uw_buoy  = Field("tke_buoy" , sloc)
uw_pres  = Field("tke_pres" , sloc)
uw_visc  = Field("tke_visc" , sloc)
uw_diss  = Field("tke_diss" , sloc)

tke      = 0.5*( interpx((u-umean)**2) + interpy((v-vmean)**2) + interpz(w**2) )

rhs_shear = interpz(w)*interpx(u-umean) * gradz(interpz(umean))*dzi4 \
          + interpz(w)*interpy(v-vmean) * gradz(interpz(vmean))*dzi4

rhs_turb = gradz( w * interpz(tke) ) * dzi4

rhs_buoy = interpx(u-umean)*(b-bmean) * sinalpha + interpz(w)*(b-bmean) * cosalpha

rhs_pres = gradz( w * interpz(p-pmean) ) * dzi4

rhs_visc = visc * gradz(gradz( tke ) * dzhi4) * dzi4

rhs_diss = 2.*visc * ( (gradx(         (u-umean) ) * dxi )**2\
                     + (grady( interpxy(u-umean) ) * dyi )**2\
                     + (gradz( interpxz(u-umean) ) * dzi4)**2\
                     + (gradx( interpxy(v-vmean) ) * dxi )**2\
                     + (grady(         (v-vmean) ) * dyi )**2\
                     + (gradz( interpyz(v-vmean) ) * dzi4)**2\
                     + (gradx(       interpxz(w) ) * dxi )**2\
                     + (grady(       interpyz(w) ) * dyi )**2\
                     + (gradz(                w  ) * dzi4)**2)

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