#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
s = Field("s", sloc)

a = Field("a", sloc)

printStencil(a, s * interpx(u), "=")
