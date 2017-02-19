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

bw_advec = Field("bw_advec", wloc)
bw_shear = Field("bw_shear", wloc)
bw_turb  = Field("bw_turb" , wloc)
bw_buoy  = Field("bw_buoy" , wloc)
bw_visc  = Field("bw_visc" , wloc)
bw_pres  = Field("bw_pres" , wloc)

rhs_advec = interpz(umean) * gradx( interpx(w) * interpxz(b-bmean)) * dxi \
          + interpz(vmean) * grady( interpy(w) * interpyz(b-bmean)) * dyi# \
          #+ w * gradz( (u-umean)*(b-bmean)) * dzih4

rhs_shear = w*interpxz(u-umean) * gradx(interpxz(bmean))*dxi \
          + w*interpyz(v-vmean) * grady(interpyz(bmean))*dyi \
          + w**2 * gradz((bmean))*dzhi4

rhs_turb = gradx( interpx(w) * interpz(u-umean) * interpxz(b-bmean) ) * dxi \
         + grady( interpy(w) * interpz(v-vmean) * interpyz(b-bmean) ) * dyi \
         + gradz( interpx(u-umean) * interpz(w) * (b-bmean) ) * dzhi4

rhs_buoy = interpz( b-bmean )**2 * sinalpha \
         - n2*( w*interpxz(u-umean) * sinalpha \
              + w*interpyz(v-vmean) \
              + w**2 * cosalpha ) 
 
rhs_visc = visc * ( interpz(b-bmean)*gradx(gradx(w) * dxi) * dxi \
                  + interpz(b-bmean)*grady(grady(w) * dyi) * dyi \
                  + interpz(b-bmean)*gradz(gradz(w) * dzi4) * dzhi4 \
                  + w*gradx(gradx(interpz(b-bmean))*dxi)*dxi \
                  + w*grady(grady(interpz(b-bmean))*dyi)*dyi \
                  + w*gradz(gradz(interpz(b-bmean))*dzi4)*dzhi4)
 
rhs_pres = interpz(b-bmean) * gradz( (p-pmean) ) * dzhi4

printStencil(bw_advec, rhs_advec, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bw_advec, rhs_advec, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bw_advec, rhs_advec, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bw_advec, rhs_advec, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bw_advec, rhs_advec, "-=", "top", "[k]")

printEmptyLine(6)

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

printStencil(bw_buoy, rhs_buoy, "+=", "int", "[k]")

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

printStencil(bw_pres, rhs_pres, "-=", "int", "[k]")