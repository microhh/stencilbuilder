#!/bin/python

import numpy as np
from StencilBuilder_advanced import *

u  = Scalar("u")
v  = Scalar("v")
w  = Scalar("w")
ut = Scalar("ut")

ut = grad( interp(u) * interp(u) ) \
   + grad( interp(v) * interp(u) ) \
   + grad( interp(w) * interp(u) )

print("ut[i] = {0};\n".format(ut.getString(0,8,0)))
