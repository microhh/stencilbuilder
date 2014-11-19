#!/bin/python
import numpy as np

# Node class
class NodeAdd:
  def __init__(self, left, right):
    self.left  = left
    self.right = right

  def __getitem__(self, i):
    return self.left[i] + self.right[i]

  def getString(self):
    return "( {0} + {1} )".format(self.left.getString(), self.right.getString())

  def __add__(self, right):
    return NodeAdd(self, right)

  def __mul__(self, right):
    return NodeMult(self, right)

class NodeMult:
  def __init__(self, left, right):
    self.left  = left
    self.right = right
  def __getitem__(self, i):
    return self.left[i] * self.right[i]

  def getString(self):
    return "( {0} * {1} )".format(self.left.getString(), self.right.getString())

  def __add__(self, right):
    return NodeAdd(self, right)

  def __mul__(self, right):
    return NodeMult(self, right)

# Scalar class representing one grid cell
class Scalar:
  def __init__(self, data, name):
    self.data = data
    self.name = name
    self.nest = 0

  def __getitem__(self, i):
    return self.data[i]

  def getString(self):
    return "{0}[i]".format(self.name)

  def __add__(self, right):
    return NodeAdd(self, right)

  def __mul__(self, right):
    return NodeMult(self, right)

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
d = ( a + b ) * c

print(d.getString())

# Check
print(a_data[3] * b_data[3] + c_data[3])
print(d[3])
