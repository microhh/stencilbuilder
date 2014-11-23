#!/bin/python

import numpy as np
from StencilBuilder import *

uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
sloc = np.array([0,0,0])

u = Field("u", uloc)
v = Field("v", vloc)
w = Field("w", wloc)
s = Field("s", sloc)

a = s * interpx(u)
print("a[i,j,k] = {0};\n".format(a.getString(0,0,0,11)))

a = grady(s * interpx(u)) + grady(s * interpx(u)) * grady(s * interpx(u))
print("a[i,j,k] = {0};\n".format(a.getString(0,0,0,11)))

