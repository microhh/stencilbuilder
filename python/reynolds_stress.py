#!/usr/bin/python

from StencilBuilder import *

u  = Field("u" , uloc )
v  = Field("v" , vloc )
w  = Field("w" , wloc )
wx = Field("wx", uwloc)

umean = Vector("umean", zloc)
vmean = Vector("vmean", zloc)

uw_diss = Field("uw_diss", uwloc)

dxi = Scalar("cgi*dxi")
dyi = Scalar("cgi*dyi")

visc = Scalar("visc")
two  = Scalar("2.")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

rhs_turb  = gradz( interpz(wx) * (u-umean) ) * dzhi4
rhs_shear = wx * gradz(umean) * dzhi4

rhs_diss_pseudo_x = two * visc * ( gradx(interpxz(u-umean)) * dxi   * gradx(interpyz(v-vmean)) * dxi   )
rhs_diss_pseudo_y = two * visc * ( grady(interpyz(u-umean)) * dyi   * grady(interpxz(v-vmean)) * dyi   )
rhs_diss_pseudo_z = two * visc * ( gradz(        (u-umean)) * dzhi4 * gradz(interpxy(v-vmean)) * dzhi4 )

rhs_diss_x1 = gradx( interpxz( u-umean ) ) * dxi \
            * ( gradx( interpyz( v-vmean ) ) * dxi + grady( interpyz( u-umean ) ) * dyi )
rhs_diss_y1 = grady( interpyz( u-umean ) ) * dyi \
            * two * grady( interpxz( v-vmean ) ) * dyi
rhs_diss_z1 = gradz( u-umean ) * dzhi4 \
            * ( gradz( interpxy( v-vmean ) ) * dzhi4 + grady( interpxy( w ) ) * dyi )

rhs_diss_x2 = gradx( interpyz( v-vmean ) ) * dxi \
            * two * gradx( interpxz( u-umean ) ) * dxi
rhs_diss_y2 = grady( interpxz( v-vmean ) ) * dyi \
            * ( grady( interpyz( u-umean ) ) * dyi + gradx( interpyz( v-vmean ) ) * dxi )
rhs_diss_z2 = gradz( interpxy( v-vmean ) ) * dzhi4 \
            * ( gradz( u-umean ) * dzhi4 + gradx( w ) * dxi )

printStencil(uw_diss, rhs_diss_x1, "-=", "int", "[k]")
printEmptyLine(1)
printStencil(uw_diss, rhs_diss_y1, "-=", "int", "[k]")
printEmptyLine(1)
printStencil(uw_diss, rhs_diss_z1, "-=", "int", "[k]")
printEmptyLine(1)
printStencil(uw_diss, rhs_diss_x2, "-=", "int", "[k]")
printEmptyLine(1)
printStencil(uw_diss, rhs_diss_y2, "-=", "int", "[k]")
printEmptyLine(1)
printStencil(uw_diss, rhs_diss_z2, "-=", "int", "[k]")
