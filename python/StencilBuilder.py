#!/usr/bin/python

import numpy as np
import copy

# Define the location arrays
uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
sloc = np.array([0,0,0])
zloc  = 0
zhloc = 1

# Base Node class
class Node(object):
  def __add__(self, right):
    return NodeOperator(self, right, "+")

  def __sub__(self, right):
    return NodeOperator(self, right, "-")

  def __mul__(self, right):
    return NodeOperator(self, right, "*")

class NodeOperator(Node):
  def __init__(self, left, right, operatorString):
    self.left  = left
    self.right = right
    self.operatorString = operatorString
    self.depth = max(left.depth, right.depth)
    if (self.depth > 1):
      self.pad = 2
    else:
      self.pad = 0

    if (type(left) == Scalar):
      self.loc = right.loc
    elif (type(right) == Scalar):
      self.loc = left.loc

    # In case of Vector, only check the k-location
    # CvH types of Vector*Vector will go wrong...
    elif (type(left) == Vector):
      if (left.loc[2] == right.loc[2]):
        self.loc = right.loc
      else:
        raise (RuntimeError)
    elif (type(right) == Vector):
      if (left.loc[2] == right.loc[2]):
        self.loc = left.loc
      else:
        raise (RuntimeError)
    elif (np.array_equal(left.loc, right.loc)):
      self.loc = copy.deepcopy(left.loc)
    else:
      raise (RuntimeError)

  def getString(self, i, j, k, pad):
    ob = '( '
    cb = ' )'

    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(1, self.depth):
        lb = lb + '\n'

      return "{ob}{0}\n{lb}{ws}{os} {1}{cb}".format(self.left.getString(i, j, k, pad),
                                                 self.right.getString(i, j, k, pad),
                                                 ws=ws, lb=lb, ob=ob, os=self.operatorString, cb=cb)
    else:
      return "{ob}{0} {os} {1}{cb}".format(self.left.getString(i, j, k, pad),
                                        self.right.getString(i, j, k, pad),
                                        ob=ob, os=self.operatorString, cb=cb)

class NodeStencilFour(Node):
  def __init__(self, inner, dim, c0, c1, c2, c3):
    self.inner = inner
    self.depth = inner.depth + 1
    self.pad = 8

    self.dim = dim
    self.loc = copy.deepcopy(inner.loc)
    self.loc[dim] = not self.loc[dim]

    self.c0 = c0
    self.c1 = c1
    self.c2 = c2
    self.c3 = c3

  def getString(self, i, j, k, pad):
    i0 = i1 = i2 = i3 = i
    j0 = j1 = j2 = j3 = j
    k0 = k1 = k2 = k3 = k
    if (self.dim == 0):
      i0 += -1-self.loc[0]
      i1 +=   -self.loc[0]
      i2 += +1-self.loc[0]
      i3 += +2-self.loc[0]
    elif (self.dim == 1):
      j0 += -1-self.loc[1]
      j1 +=   -self.loc[1]
      j2 += +1-self.loc[1]
      j3 += +2-self.loc[1]
    elif (self.dim == 2):
      k0 += -1-self.loc[2]
      k1 +=   -self.loc[2]
      k2 += +1-self.loc[2]
      k3 += +2-self.loc[2]

    ob = '( '
    cb = ' )'

    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "{ob}{c0}*{0}\n{lb}{ws}+ {c1}*{1}\n{lb}{ws}+ {c2}*{2}\n{lb}{ws}+ {c3}*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad),
          self.inner.getString(i1, j1, k1, pad),
          self.inner.getString(i2, j2, k2, pad),
          self.inner.getString(i3, j3, k3, pad),
          ws=ws, lb=lb, ob=ob, cb=cb,
          c0=self.c0, c1=self.c1, c2=self.c2, c3=self.c3)
    else:
      return "{ob}{c0}*{0} + {c1}*{1} + {c2}*{2} + {c3}*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad),
          self.inner.getString(i1, j1, k1, pad),
          self.inner.getString(i2, j2, k2, pad),
          self.inner.getString(i3, j3, k3, pad),
          ob=ob, cb=cb,
          c0=self.c0, c1=self.c1, c2=self.c2, c3=self.c3)

def formatIndex(n, nstr):
  if (n > 0):
    nn = "+{0}{1}".format(nstr, abs(n))
  elif (n < 0):
    nn = "-{0}{1}".format(nstr, abs(n))
  else:
    nn = " " * (len(nstr) + 2)
  return nn

class Field(Node):
  def __init__(self, name, loc):
    self.name  = name
    self.depth = 0
    self.loc   = loc

  def getString(self, i, j, k, pad):
    ii = formatIndex(i, "ii")
    jj = formatIndex(j, "jj")
    kk = formatIndex(k, "kk")
    return "{0}[ijk{1}{2}{3}]".format(self.name, ii, jj, kk)

# Vector class representing a profile
class Vector(Node):
  def __init__(self, name, loc):
    self.name  = name
    self.depth = 0
    self.loc   = np.array([0, 0, loc])

  def getString(self, i, j, k, pad):
    kk = formatIndex(k, "")
    return "{0}[k{1}]".format(self.name, kk)

class Scalar(Node):
  def __init__(self, name):
    self.name  = name
    self.depth = 0

  def getString(self, i, j, k, pad):
    return "{0}".format(self.name)

# Define functions.
def interpx(inner):
  return NodeStencilFour(inner, 0, "ci0", "ci1", "ci2", "ci3")
def interpy(inner):
  return NodeStencilFour(inner, 1, "ci0", "ci1", "ci2", "ci3")
def interpz(inner):
  return NodeStencilFour(inner, 2, "ci0", "ci1", "ci2", "ci3")

def gradx(inner):
  return NodeStencilFour(inner, 0, "cg0", "cg1", "cg2", "cg3")
def grady(inner):
  return NodeStencilFour(inner, 1, "cg0", "cg1", "cg2", "cg3")
def gradz(inner):
  return NodeStencilFour(inner, 2, "cg0", "cg1", "cg2", "cg3")

def printStencil(lhs, rhs, operator):
  if (not np.array_equal(lhs.loc, rhs.loc)):
    raise (RuntimeError)

  index = "[ijk]"
  indent = len(lhs.name) + len(index) + len(operator) + 2
  print("{0}{1} {2} {3};".format(lhs.name, index, operator, rhs.getString(0, 0, 0, indent)))

