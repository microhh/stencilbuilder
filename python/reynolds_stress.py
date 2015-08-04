#!/usr/bin/python

from StencilBuilder import *

p  = Field("p" , sloc )
u  = Field("u" , uloc )
v  = Field("v" , vloc )
w  = Field("w" , wloc )
wx = Field("wx", uwloc)

umean = Vector("umean", zloc)
vmean = Vector("vmean", zloc)

uw_pres = Field("uw_pres", uwloc)

dxi = Scalar("cgi*dxi")
dyi = Scalar("cgi*dyi")

visc = Scalar("visc")
two  = Scalar("2.")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

rhs_turb  = gradz( interpz(wx) * (u-umean) ) * dzhi4
rhs_shear = wx * gradz(umean) * dzhi4

rhs_diss_pseudo_x = ( gradx(interpxz(u-umean)) * dxi   * gradx(w) * dxi   )
rhs_diss_pseudo_y = ( grady(interpyz(u-umean)) * dyi   * grady(interpxy(w)) * dyi   )
rhs_diss_pseudo_z = ( gradz(        (u-umean)) * dzhi4 * gradz(interpxz(w)) * dzhi4 )

rhs_rdstr = interpxz(p) * ( gradz(u-umean) * dzhi4 + gradx(w) * dxi )

rhs_pres = gradz( (u-umean) * interpx(p) ) * dzhi4 + gradx( w * interpz(p) ) * dxi

#rhs_diss_x1 = gradx( interpxz( u-umean ) ) * dxi \
#            * ( gradx( interpyz( v-vmean ) ) * dxi + grady( interpyz( u-umean ) ) * dyi )
#rhs_diss_y1 = grady( interpyz( u-umean ) ) * dyi \
#            * two * grady( interpxz( v-vmean ) ) * dyi
#rhs_diss_z1 = gradz( u-umean ) * dzhi4 \
#            * ( gradz( interpxy( v-vmean ) ) * dzhi4 + grady( interpxy( w ) ) * dyi )
#
#rhs_diss_x2 = gradx( interpyz( v-vmean ) ) * dxi \
#            * two * gradx( interpxz( u-umean ) ) * dxi
#rhs_diss_y2 = grady( interpxz( v-vmean ) ) * dyi \
#            * ( grady( interpyz( u-umean ) ) * dyi + gradx( interpyz( v-vmean ) ) * dxi )
#rhs_diss_z2 = gradz( interpxy( v-vmean ) ) * dzhi4 \
#            * ( gradz( u-umean ) * dzhi4 + gradx( w ) * dxi )

printStencil(uw_pres, rhs_pres, "-=", "int", "[k]")
