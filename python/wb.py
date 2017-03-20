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
N2       = Scalar("N2")
dxi      = Scalar("cgi*dxi")
dyi      = Scalar("cgi*dyi")

dzi4     = Vector("dzi4" , zloc )
dzhi4    = Vector("dzhi4", zhloc)

bw_shear = Field("bw_shear", wloc)
bw_turb  = Field("bw_turb" , wloc)
bw_pres  = Field("bw_pres" , wloc)
bw_rdstr = Field("bw_rdstr", wloc)
bw_buoy  = Field("bw_buoy" , wloc)
bw_visc  = Field("bw_visc" , wloc)
bw_diss  = Field("bw_diss" , wloc)

rhs_shear = w**2 * gradz((bmean))*dzhi4 \
          + N2*( interpxz(u-umean)*w * sinalpha + w**2 * cosalpha ) 

rhs_turb = gradz( interpz(w)**2 * (b-bmean) ) * dzhi4

rhs_buoy = interpz( b-bmean )**2 * cosalpha

rhs_pres = gradz( (p-pmean) * (b-bmean) ) * dzhi4

rhs_rdstr = interpz( p-pmean ) * gradz( b-bmean ) * dzhi4

rhs_visc = visc * gradz(gradz(w * interpz(b-bmean) ) * dzi4 ) * dzhi4

rhs_diss = 2.*visc * ( gradx( interpx(w) ) * dxi   * gradx( interpxz( b-bmean ) ) * dxi \
                     + grady( interpy(w) ) * dyi   * grady( interpyz( b-bmean ) ) * dyi \
                     + gradz( interpz(w) ) * dzhi4 * gradz( b-bmean ) * dzhi4 )

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

printStencil(bw_pres, rhs_pres, "-=", "int", "[k]")

printEmptyLine(6)

printStencil(bw_rdstr, rhs_rdstr, "+=", "int", "[k]")

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

printStencil(bw_diss, rhs_diss, "-=", "bot", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "int", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bw_diss, rhs_diss, "-=", "top", "[k]")