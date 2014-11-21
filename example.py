#!/bin/python

import numpy as np
from StencilBuilder import *

uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
sloc = np.array([0,0,0])

u = Scalar("u", uloc)
v = Scalar("v", vloc)
w = Scalar("w", wloc)
s = Scalar("s", sloc)

a = Scalar("a", sloc)

a = s * interpx(u)
print("a[i] = {0};\n".format(a.getString(0,0,0,7,0)))

a = grady(s * interpx(u)) + grady(s * interpx(u)) * grady(s * interpx(u))
print("a[i] = {0};\n".format(a.getString(0,0,0,7,0)))

