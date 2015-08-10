#!/usr/bin/python

from StencilBuilder import *

b = Field("b", sloc)
w = Field("w", wloc)

bt = Field("btest", sloc)
wt = Field("wtest", wloc)

brhs = gradz( gradz( b ) )
wrhs = gradz( gradz( w ) )

printStencil(bt, brhs, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(bt, brhs, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(bt, brhs, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(bt, brhs, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(bt, brhs, "+=", "top", "[k]")

printEmptyLine(3)
print("BLA")
printEmptyLine(3)

printStencil(wt, wrhs, "+=", "bot", "[k]")
printEmptyLine(3)
printStencil(wt, wrhs, "+=", "bot+1", "[k]")
printEmptyLine(3)
printStencil(wt, wrhs, "+=", "int", "[k]")
printEmptyLine(3)
printStencil(wt, wrhs, "+=", "top-1", "[k]")
printEmptyLine(3)
printStencil(wt, wrhs, "+=", "top", "[k]")
