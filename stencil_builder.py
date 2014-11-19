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

  def __getitem__(self, i):
    return self.left[i] + self.right[i]

  def getString(self, i):
    return "( {0} + {1} )".format(self.left.getString(i), self.right.getString(i))

class NodeMult(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right

  def __getitem__(self, i):
    return self.left[i] * self.right[i]

  def getString(self, i):
    return "( {0} * {1} )".format(self.left.getString(i), self.right.getString(i))

class NodeStencil(Node):
  def __init__(self, inner):
    self.inner = inner

  def __getitem__(self, i):
    ci0 = 1.
    ci1 = 1.
    ci2 = 1.
    ci3 = 1.
    return ci0*self.inner[i-2] + ci1*self.inner[i-1] + ci2*self.inner[i] + ci3*self.inner[i+1]

  def getString(self, i):
    return "( ci0*{0} + ci1*{1} + ci2*{2} + ci3*{3} )".format(self.inner.getString(i-2), self.inner.getString(i-1), self.inner.getString(i), self.inner.getString(i+1))

# Scalar class representing one grid cell
class Scalar(Node):
  def __init__(self, data, name):
    self.data = data
    self.name = name
    self.nest = 0

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

# d = ScalarAdd( ScalarMult(a,b), c )
d = interp( interp(a) + interp(b) ) * c
print("d[i] = {0};".format(d.getString(3)))
print(d[3])
