#!/bin/python
import numpy as np

# Base Node class
class Node:
  def __add__(self, right):
    return NodeAdd(self, right)

  def __mul__(self, right):
    return NodeMult(self, right)

class NodeAdd(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right
    self.depth = max(left.depth, right.depth)
    self.pad   = 2

  def __getitem__(self, i):
    return self.left[i] + self.right[i]

  def getString(self, i, pad):
    pad += self.pad
    return "( {0} + {1} )".format(self.left.getString(i, pad), self.right.getString(i, pad))

class NodeMult(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right
    self.depth = max(left.depth, right.depth)
    self.pad   = 2

  def __getitem__(self, i):
    return self.left[i] * self.right[i]

  def getString(self, i, pad):
    pad += self.pad
    return "( {0} * {1} )".format(self.left.getString(i, pad), self.right.getString(i, pad))

class NodeStencil(Node):
  def __init__(self, inner):
    self.inner = inner
    self.depth = inner.depth + 1
    self.pad = 6

  def __getitem__(self, i):
    ci0 = 1.
    ci1 = 1.
    ci2 = 1.
    ci3 = 1.
    return ci0*self.inner[i-2] + ci1*self.inner[i-1] + ci2*self.inner[i] + ci3*self.inner[i+1]

  def getString(self, i, pad):
    if(self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad
      return "( ci0*{0}\n{ws}+ ci1*{1}\n{ws}+ ci2*{2}\n{ws}+ ci3*{3} )".format(self.inner.getString(i-2, pad), self.inner.getString(i-1, pad), self.inner.getString(i, pad), self.inner.getString(i+1, pad), ws=ws)
    else:
      return "( ci0*{0} + ci1*{1} + ci2*{2} + ci3*{3} )".format(self.inner.getString(i-2, pad), self.inner.getString(i-1, pad), self.inner.getString(i, pad), self.inner.getString(i+1, pad))

# Scalar class representing one grid cell
class Scalar(Node):
  def __init__(self, data, name):
    self.data  = data
    self.name  = name
    self.depth = 0

  def __getitem__(self, i):
    return self.data[i]

  def getString(self, i, pad):
    nn = i-3
    if(nn > 0):
      return "{0}[i+{1}]".format(self.name, nn)
    elif(nn < 0):
      return "{0}[i-{1}]".format(self.name, abs(nn))
    else:
      return "{0}[i  ]".format(self.name, abs(nn))

# Define functions.
def interp(inner):
  return NodeStencil(inner)

# test
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
