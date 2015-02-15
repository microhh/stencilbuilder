#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)

ut = Field("ut", uloc)

dxi = Scalar("dxi")
dyi = Scalar("dyi")
dzi4 = Vector("dzi4", zloc)

utrhsb = gradx( interpx(u) * interpx (u) ) * dxi \
       + grady( interpx(v) * interpy (u) ) * dyi \
       + gradz( interpx(w) * interpzb(u) ) * dzi4

utrhs = gradx( interpx(u) * interpx(u) ) * dxi \
      + grady( interpx(v) * interpy(u) ) * dyi \
      + gradz( interpx(w) * interpz(u) ) * dzi4

utrhst = gradx( interpx(u) * interpx (u) ) * dxi \
       + grady( interpx(v) * interpy (u) ) * dyi \
       + gradz( interpx(w) * interpzt(u) ) * dzi4

printStencil(ut, utrhsb, "=")
printStencil(ut, utrhs , "=")
printStencil(ut, utrhst, "=")

printLoop(ut, utrhsb, "=", kstart="grid->kstart"  , kend="grid->kstart+1")
printLoop(ut, utrhs , "=", kstart="grid->kstart+1", kend="grid->kend-1"  )
printLoop(ut, utrhst, "=", kstart="grid->kend-1"  , kend="grid->kend"    )
