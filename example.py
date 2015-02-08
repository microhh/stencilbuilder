#!/usr/bin/python

from StencilBuilder import *

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)
s = Field("s", sloc)

a = s * interpx(u)
print("a[i,j,k] = {0};\n".format(a.getString(0,0,0,11)))

a = grady(s * interpx(u)) + grady(s * interpx(u)) * grady(s * interpx(u))
print("a[i,j,k] = {0};\n".format(a.getString(0,0,0,11)))

