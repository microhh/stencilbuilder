#!/bin/python

import numpy as np
from StencilBuilder import *

uloc = 1
sloc = 0

u  = Scalar("u" , uloc)
v  = Scalar("v" , uloc)
w  = Scalar("w" , uloc)
ut = Scalar("ut", uloc)

ut = grad( interp(u) * interp(u) ) \
   + grad( interp(v) * interp(u) ) \
   + grad( interp(w) * interp(u) )

print("ut[i] = {0};\n".format(ut.getString(0,8,0)))
