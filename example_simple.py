#!/bin/python

import numpy as np
from StencilBuilder_simple import *

a_data = np.random.uniform(0., 1., 7)
b_data = np.random.uniform(0., 1., 7)
c_data = np.random.uniform(0., 1., 7)
d_data = np.zeros(7)

a = Scalar(a_data, "a")
b = Scalar(b_data, "b")
c = Scalar(c_data, "c")
d = Scalar(d_data, "d")

print("interp( interp(a)")
d = interp( interp(a))
print("d[i] = {0};\n".format(d.getString(3,7,0)))

print("interp( a * b * c )")
d = interp( a * b * c )
print("d[i] = {0};\n".format(d.getString(3,7,0)))

print("interp( interp(a) + interp(b) ) * c")
d = interp( interp(a) + interp(b) ) * c
print("d[i] = {0};\n".format(d.getString(3,7,0)))

print("a * interp( interp(b) + interp(c) )")
d = a * interp( interp(b) + interp(c) )
print("d[i] = {0};\n".format(d.getString(3,7,0)))

print("interp( c * interp(a) + interp(b) )")
d = interp( c * interp(a) + interp(b) )
print("d[i] = {0};\n".format(d.getString(3,7,0)))

print("interp( interp( interp(a) + interp(b) ) )")
d = interp( interp( interp(a) + interp(b) ) )
print("d[i] = {0};\n".format(d.getString(3,7,0)))

print("interp( interp( interp( interp(a) ) ) )")
d = interp( interp( interp( interp(a) ) ) )
print("d[i] = {0};\n".format(d.getString(3,7,0)))

print("interp( interp(a) * interp(b) ) + interp( interp(a) * interp(c)")
d = interp( interp(a) * interp(b) ) + interp( interp(a) * interp(c) )
print("d[i] = {0};\n".format(d.getString(3,7,0)))
