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

  def getString(self, i, pad):
    pad += self.pad
    return "( {0} + {1} )".format(self.left.getString(i, pad), self.right.getString(i, pad))

class NodeMult(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right
    self.depth = max(left.depth, right.depth)
    self.pad   = 2

  def getString(self, i, pad):
    pad += self.pad
    return "( {0} * {1} )".format(self.left.getString(i, pad), self.right.getString(i, pad))

class NodeStencil(Node):
  def __init__(self, inner):
    self.inner = inner
    self.depth = inner.depth + 1
    if (type(inner) == Scalar):
      self.pad = 6
    else:
      self.pad = 8

class NodeStencilInterp(NodeStencil):
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

class NodeStencilGrad(NodeStencil):
  def getString(self, i, pad):
    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "( cg0 * {0}\n{lb}{ws}+ cg1 * {1}\n{lb}{ws}+ cg2 * {2}\n{lb}{ws}+ cg3 * {3} )".format(self.inner.getString(i-2, pad), self.inner.getString(i-1, pad), self.inner.getString(i, pad), self.inner.getString(i+1, pad), ws=ws, lb=lb)
    else:
      return "( cg0*{0} + cg1*{1} + cg2*{2} + cg3*{3} )".format(self.inner.getString(i-2, pad), self.inner.getString(i-1, pad), self.inner.getString(i, pad), self.inner.getString(i+1, pad))

# Scalar class representing one grid cell
class Scalar(Node):
  def __init__(self, name):
    self.name  = name
    self.depth = 0

  def getString(self, i, pad):
    if (i > 0):
      return "{0}[i+{1}]".format(self.name, i)
    elif (i < 0):
      return "{0}[i-{1}]".format(self.name, abs(i))
    else:
      return "{0}[i  ]".format(self.name, abs(i))

  def printStencil(self, pad):
    self.getString(0, pad)

# Define functions.
def interp(inner):
  return NodeStencilInterp(inner)

def grad(inner):
  return NodeStencilGrad(inner)

