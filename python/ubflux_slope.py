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
n2       = Scalar("n2")
dxi      = Scalar("cgi*dxi")
dyi      = Scalar("cgi*dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

bu_advec = Field("bu_advec", uloc)
bu_shear = Field("bu_shear", uloc)
bu_turb  = Field("bu_turb" , uloc)
bu_buoy  = Field("bu_buoy" , uloc)
bu_visc  = Field("bu_visc" , uloc)
bu_pres  = Field("bu_pres" , uloc)

rhs_advec = umean * gradx( interpx(u-umean) * (b-bmean)) * dxi \
          + vmean * grady( interpy(u-umean) * interpxy(b-bmean)) * dyi# \
          #+ w * gradz( (u-umean)*(b-bmean)) * dzih4

rhs_shear = (u-umean)*interpx(b-bmean)*gradx(umean)*dxi \
          + interpxy(v-vmean)*interpx(b-bmean)*grady(interpy(umean))*dyi \
          + interpxz(w)*interpx(b-bmean)*gradz(interpz(umean))*dzi4 \
          + (u-umean)**2 * gradx(bmean)*dxi \
          + (u-umean)*interpxy(v-vmean) * grady(interpxy(bmean))*dyi \
          + (u-umean)*interpxz(w) * gradz(interpxz(bmean))*dzi4

rhs_turb = gradx( interpx(u-umean)**2 * (b-bmean) ) * dxi \
         + grady( interpx(v-vmean) * interpy(u-umean) * interpxy(b-bmean) ) * dyi \
         + gradz( interpz(u-umean) * interpx(w) * interpxz(b-bmean) ) * dzi4

rhs_buoy = interpx( b-bmean )**2 * sinalpha \
         - n2*((u-umean)**2 * sinalpha \
             + (u-umean)*interpxy(v-vmean) \
             + (u-umean)*interpxz(w)*cosalpha) 
 
rhs_visc = visc * ( interpx(b-bmean)*gradx(gradx((u-umean)) * dxi) * dxi \
                  + interpx(b-bmean)*grady(grady((u-umean)) * dyi) * dyi \
                  + interpx(b-bmean)*gradz(gradz((u-umean)) * dzhi4) * dzi4 \
                  + (u-umean)*gradx(gradx(interpx(b-bmean))*dxi)*dxi \
                  + (u-umean)*grady(grady(interpx(b-bmean))*dyi)*dyi \
                  + (u-umean)*gradz(gradz(interpx(b-bmean))*dzhi4)*dzi4)
 
rhs_pres = interpx(b-bmean) * gradx( (p-pmean) ) * dxi

printStencil(bu_advec, rhs_advec, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bu_advec, rhs_advec, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bu_advec, rhs_advec, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bu_advec, rhs_advec, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bu_advec, rhs_advec, "-=", "top", "[k]")

printEmptyLine(6)

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

printStencil(bu_pres, rhs_pres, "-=", "int", "[k]")