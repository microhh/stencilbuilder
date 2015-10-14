#!/usr/bin/python

from StencilBuilder import *

p  = Field("p" , sloc )
b  = Field("b" , sloc )
u  = Field("u" , uloc )
v  = Field("v" , vloc )
w  = Field("w" , wloc )

wx = Field("wx", uwloc)
uz = Field("uz", uwloc)

umean = Vector("umean", zloc)
vmean = Vector("vmean", zloc)
bmean = Vector("bmean", zloc)

uw_turb  = Field("uw_turb" , uwloc)
uw_buoy  = Field("uw_buoy" , uwloc)
uw_shear = Field("uw_shear", uwloc)
uw_diss  = Field("uw_diss" , uwloc)
uw_rdstr = Field("uw_rdstr", uwloc)

dxi = Scalar("cgi*dxi")
dyi = Scalar("cgi*dyi")

visc = Scalar("visc")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

rhs_turb  = gradz( interpz(wx)**2 * (u-umean) ) * dzhi4
rhs_shear = wx * gradz(umean) * dzhi4

rhs_diss_x = ( gradx(interpxz(u-umean)) * dxi   * gradx(w) * dxi   )
rhs_diss_y = ( grady(interpyz(u-umean)) * dyi   * grady(interpxy(w)) * dyi   )
rhs_diss_z = ( gradz(        (u-umean)) * dzhi4 * gradz(interpxz(w)) * dzhi4 )

rhs_diss = 2. * visc * (rhs_diss_x + rhs_diss_y + rhs_diss_z)

rhs_rdstr = interpxz(p) * ( gradz(u-umean) * dzhi4 + gradx(w) * dxi )

rhs_pres = gradz( (u-umean) * interpx(p) ) * dzhi4 + gradx( w * interpz(p) ) * dxi

rhs_visc = visc * gradz( gradz( uz * interpx(w) ) * dzi4 ) * dzhi4

rhs_buoy = interpz( u-umean ) * interpxz( b-bmean )

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

#printStencil(uw_turb , rhs_turb  , "-=", "int"  , "[k]")
#printEmptyLine(6)
#printStencil(uw_buoy , rhs_buoy  , "+=", "int"  , "[k]")
#printEmptyLine(6)
#printStencil(uw_shear, rhs_shear , "-=", "int"  , "[k]")
#printEmptyLine(6)
#printStencil(uw_diss , rhs_diss_x, "-=", "int"  , "[k]")
#printEmptyLine(6)
#printStencil(uw_rdstr, rhs_rdstr , "+=", "int"  , "[k]")

printStencil(uw_diss, rhs_diss, "-=", "bot"  , "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "int"  , "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(uw_diss, rhs_diss, "-=", "top"  , "[k]")
