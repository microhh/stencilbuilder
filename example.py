#!/bin/python

import numpy as np
from StencilBuilder import *

u  = Scalar("a")
v  = Scalar("b")
w  = Scalar("c")
ut = Scalar("d")

print( "duu/dx")
ut = grad( interp(u) * interp(u) )
print("ut[i] = {0};".format(ut.getString(0, 8)))
