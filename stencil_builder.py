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
    self.nest  = left.nest + 1

  def __getitem__(self, i):
    return self.left[i] + self.right[i]

  def getString(self, i):
    return "( {0} + {1} )".format(self.left.getString(i), self.right.getString(i))

class NodeMult(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right
    self.depth = max(left.depth, right.depth)
    self.nest  = left.nest + 1

  def __getitem__(self, i):
    return self.left[i] * self.right[i]

  def getString(self, i):
    return "( {0} * {1} )".format(self.left.getString(i), self.right.getString(i))

class NodeStencil(Node):
  def __init__(self, inner):
    self.inner = inner
    self.depth = inner.depth + 1
    self.nest  = inner.nest + 1

  def __getitem__(self, i):
    ci0 = 1.
    ci1 = 1.
    ci2 = 1.
    ci3 = 1.
    return ci0*self.inner[i-2] + ci1*self.inner[i-1] + ci2*self.inner[i] + ci3*self.inner[i+1]

  def getString(self, i):
    if(self.depth > 1):
      pad = (self.depth-1)*2
      return "( ci0*{0}\n{4:{pad}}+ ci1*{1}\n{4:{pad}}+ ci2*{2}\n{4:{pad}}+ ci3*{3} )".format(self.inner.getString(i-2), self.inner.getString(i-1), self.inner.getString(i), self.inner.getString(i+1), ' ', pad = pad)
    else:
      return "( ci0*{0} + ci1*{1} + ci2*{2} + ci3*{3} )".format(self.inner.getString(i-2), self.inner.getString(i-1), self.inner.getString(i), self.inner.getString(i+1))

# Scalar class representing one grid cell
class Scalar(Node):
  def __init__(self, data, name):
    self.data  = data
    self.name  = name
    self.depth = 0
    self.nest  = 0

  def __getitem__(self, i):
    return self.data[i]

  def getString(self, i):
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
print("{0};".format(d.getString(3)))


print( "interp( c * interp(a) + interp(b) )")
d = interp( c * interp(a) + interp(b) )
print("{0};".format(d.getString(3)))

print( "interp( interp( interp(a) + interp(b) ) )")
d = interp( interp( interp(a) + interp(b) ) )
print("{0};".format(d.getString(3)))

