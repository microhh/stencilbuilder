#!/bin/python

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
    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(1, self.depth):
        lb = lb + '\n'

      return "( {0}\n{lb}{ws}+ {1} )".format(self.left.getString(i, pad), self.right.getString(i, pad), ws=ws, lb=lb)

    else:
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
    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(1, self.depth):
        lb = lb + '\n'

      return "( {0}\n{lb}{ws}* {1} )".format(self.left.getString(i, pad), self.right.getString(i, pad), ws=ws, lb=lb)

    else:
      return "( {0} * {1} )".format(self.left.getString(i, pad), self.right.getString(i, pad))

class NodeStencil(Node):
  def __init__(self, inner):
    self.inner = inner
    self.depth = inner.depth + 1
    if (type(inner) == Scalar):
      self.pad = 6
    else:
      self.pad = 8

  def __getitem__(self, i):
    ci0 = 1.
    ci1 = 1.
    ci2 = 1.
    ci3 = 1.
    return ci0*self.inner[i-2] + ci1*self.inner[i-1] + ci2*self.inner[i] + ci3*self.inner[i+1]

  def getString(self, i, pad):
    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "( ci0 * {0}\n{lb}{ws}+ ci1 * {1}\n{lb}{ws}+ ci2 * {2}\n{lb}{ws}+ ci3 * {3} )".format(self.inner.getString(i-2, pad), self.inner.getString(i-1, pad), self.inner.getString(i, pad), self.inner.getString(i+1, pad), ws=ws, lb=lb)
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
    if (nn > 0):
      return "{0}[i+{1}]".format(self.name, nn)
    elif (nn < 0):
      return "{0}[i-{1}]".format(self.name, abs(nn))
    else:
      return "{0}[i  ]".format(self.name, abs(nn))

# Define functions.
def interp(inner):
  return NodeStencil(inner)
