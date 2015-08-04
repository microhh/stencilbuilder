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

visc = Scalar("2.*visc")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

rhs_turb  = gradz( interpz(wx) * (u-umean) ) * dzhi4
rhs_shear = wx * gradz(umean) * dzhi4
#rhs_diss  = visc * ( gradx(interpz(interpx(u-umean))) * dxi   * gradx(interpz(interpy(v-vmean))) * dxi \
#                   + grady(interpz(interpy(u-umean))) * dyi   * grady(interpz(interpx(v-vmean))) * dyi \
#                   + gradz(               (u-umean) ) * dzhi4 * gradz(interpx(interpy(v-vmean))) * dzhi4 )
rhs_diss = gradz(interpx(interpy(v))) * dzhi4

uw_rhs = rhs_diss

printStencil(uw_diss, uw_rhs, "-=", "int", "[k]")
