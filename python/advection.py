#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)
wt = Field("wt", wloc)

dxi = Scalar("dxi")
dyi = Scalar("dyi")

dzi4  = Vector("dzi4" , zloc )
dzhi4 = Vector("dzhi4", zhloc)

utrhs = gradx( interpx(u) * interpx(u) ) * dxi \
      + grady( interpx(v) * interpy(u) ) * dyi \
      + gradz( interpx(w) * interpz(u) ) * dzi4

wtrhs = gradx( interpz(u) * interpx(w) ) * dxi \
      + grady( interpz(v) * interpy(w) ) * dyi \
      + gradz( interpz(w) * interpz(w) ) * dzhi4

printStencil(ut, utrhs, "=", "bot"   )
#printStencil(ut, utrhs, "=", "interior")
#printStencil(ut, utrhs, "=", "end"     )

printStencil(wt, wtrhs, "=", "bot"  )

#printStencil(wt, wtrhs, "=", "interior")
#printStencil(wt, wtrhs, "=", "end"     )
#printLoop(ut, utrhs, "=", kstart="grid->kstart+1", kend="grid->kend-1"  )
