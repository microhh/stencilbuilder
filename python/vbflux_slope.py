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

bv_advec = Field("bv_advec", vloc)
bv_shear = Field("bv_shear", vloc)
bv_turb  = Field("bv_turb" , vloc)
bv_buoy  = Field("bv_buoy" , vloc)
bv_visc  = Field("bv_visc" , vloc)
bv_pres  = Field("bv_pres" , vloc)

rhs_advec = umean * gradx( interpx(v-vmean) * interpxy(b-bmean)) * dxi \
          + vmean * grady( interpy(v-vmean) * (b-bmean)) * dyi# \
          #+ w * gradz( (u-umean)*(b-bmean)) * dzih4

rhs_shear = interpxy(u-umean)*interpy(b-bmean)*gradx(interpx(vmean))*dxi \
          + (v-vmean)*interpy(b-bmean)*grady(vmean)*dyi \
          + interpyz(w)*interpy(b-bmean)*gradz(interpz(vmean))*dzi4 \
          + (v-vmean)*interpxy(u-umean) * gradx(interpxy(bmean))*dxi \
           + (v-umean)**2 * grady(bmean)*dyi \
           + (v-umean)*interpyz(w) * gradz(interpyz(bmean))*dzi4

rhs_turb = gradx( interpx(v-vmean)*interpy(u-umean)*interpxy(b-bmean) ) * dxi \
         + grady( interpy(v-vmean)**2  * (b-bmean) ) * dyi \
         + gradz( interpz(v-vmean) * interpy(w) * interpyz(b-bmean) ) * dzi4
 
rhs_buoy = n2*( (v-vmean)*interpxy(u-umean) * sinalpha \
              + (v-vmean)**2 \
              + (v-vmean)*interpyz(w)*cosalpha ) 
 
rhs_visc = visc * ( interpy(b-bmean)*gradx(gradx((v-vmean)) * dxi) * dxi \
                  + interpy(b-bmean)*grady(grady((v-vmean)) * dyi) * dyi \
                  + interpy(b-bmean)*gradz(gradz((v-vmean)) * dzhi4) * dzi4 \
                  + (v-vmean)*gradx(gradx(interpy(b-bmean))*dxi)*dxi \
                  + (v-vmean)*grady(grady(interpy(b-bmean))*dyi)*dyi \
                  + (v-vmean)*gradz(gradz(interpy(b-bmean))*dzhi4)*dzi4)
  
rhs_pres = interpy(b-bmean) * grady( (p-pmean) ) * dyi

printStencil(bv_advec, rhs_advec, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bv_advec, rhs_advec, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bv_advec, rhs_advec, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bv_advec, rhs_advec, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bv_advec, rhs_advec, "-=", "top", "[k]")

printEmptyLine(6)

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

printStencil(bv_buoy, rhs_buoy, "s=", "int", "[k]")

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

printStencil(bv_pres, rhs_pres, "-=", "int", "[k]")