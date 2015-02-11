#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)
evisc = Field("evisc", sloc)

ut = Field("ut", uloc)

dxi = Scalar("dxi")
dyi = Scalar("dyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

eviscy = interpy( interpx(evisc) )
eviscz = interpz( interpx(evisc) )

utrhs = gradx( evisc  * ( gradx(u)*dxi   + gradx(u)*dxi ) ) \
      + grady( eviscy * ( grady(u)*dyi   + gradx(v)*dxi ) ) \
      + gradz( eviscz * ( gradz(u)*dzhi4 + gradx(w)*dxi ) )

printStencil(ut, utrhs, "=")
