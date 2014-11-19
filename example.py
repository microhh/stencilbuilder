#!/bin/python

import numpy as np
from StencilBuilder import *

a_data = np.random.uniform(0., 1., 7)
b_data = np.random.uniform(0., 1., 7)
c_data = np.random.uniform(0., 1., 7)
d_data = np.zeros(7)

a = Scalar(a_data, "a")
b = Scalar(b_data, "b")
c = Scalar(c_data, "c")
d = Scalar(d_data, "d")

print( "interp( interp(a) + interp(b) ) * c")
d = interp( interp(a) + interp(b) ) * c
print("d[i] = {0};".format(d.getString(3, 7)))

print( "interp( c * interp(a) + interp(b) )")
d = interp( c * interp(a) + interp(b) )
print("d[i] = {0};".format(d.getString(3, 7)))

print( "interp( interp( interp(a) + interp(b) ) )")
d = interp( interp( interp(a) + interp(b) ) )
print("d[i] = {0};".format(d.getString(3, 7)))

print( "interp( interp( interp( interp(a) ) ) )" )
d = interp( interp( interp( interp(a) ) ) )
print("d[i] = {0};".format(d.getString(3, 7)))
