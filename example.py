#!/bin/python

import numpy as np
from StencilBuilder import *

u  = Scalar("u")
v  = Scalar("v")
w  = Scalar("w")
ut = Scalar("ut")

print( "duu/dx")
ut = grad( interp(u) * interp(u) ) + grad( interp(v) * interp(u) ) + grad( interp(w) * interp(u) )
print("ut[i] = {0};".format(ut.getString(0,8)))
