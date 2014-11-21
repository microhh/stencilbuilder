#!/bin/python
import numpy as np
import copy

# Base Node class
class Node(object):
  def __add__(self, right):
    return NodeAdd(self, right)

  def __mul__(self, right):
    return NodeMult(self, right)

class NodeAdd(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right
    self.depth = max(left.depth, right.depth)
    if (self.depth > 1):
      self.pad = 2
    else:
      self.pad = 0

    if (np.array_equal(left.loc, right.loc)):
      self.loc = copy.deepcopy(left.loc)
    else:
      raise (RuntimeError)

  def getString(self, i, j, k, pad, maxDepth):
    maxDepth = max(self.depth, maxDepth)
    if (self.depth == maxDepth):
      pad = pad - 2
      ob = ''
      cb = ''
    else:
      ob = '( '
      cb = ' )'

    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(1, self.depth):
        lb = lb + '\n'

      return "{ob}{0}\n{lb}{ws}+ {1}{cb}".format(self.left.getString(i, j, k, pad, maxDepth),
                                                 self.right.getString(i, j, k, pad, maxDepth),
                                                 ws=ws, lb=lb, ob=ob, cb=cb)
    elif (type(self.left) == Scalar and type(self.right) == Scalar):
      return "{ob}{0}+{1}{cb}".format(self.left.getString(i, j, k, pad, maxDepth),
                                      self.right.getString(i, j, k, pad, maxDepth),
                                      ob=ob, cb=cb)

    else:
      return "{ob}{0} + {1}{cb}".format(self.left.getString(i, j, k, pad, maxDepth),
                                        self.right.getString(i, j, k, pad, maxDepth),
                                        ob=ob, cb=cb)

class NodeMult(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right
    self.depth = max(left.depth, right.depth)
    if (self.depth > 1):
      self.pad = 2
    else:
      self.pad = 0

    if (np.array_equal(left.loc, right.loc)):
      self.loc = copy.deepcopy(left.loc)
    else:
      raise (RuntimeError)

  def getString(self, i, j, k, pad, maxDepth):
    maxDepth = max(self.depth, maxDepth)
    if (self.depth == maxDepth):
      pad = pad - 2
      ob = ''
      cb = ''
    else:
      ob = '( '
      cb = ' )'

    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(1, self.depth, maxDepth):
        lb = lb + '\n'

      return "{ob}{0}\n{lb}{ws}* {1}{cb}".format(self.left.getString(i, j, k, pad, maxDepth),
                                                 self.right.getString(i, j, k, pad, maxDepth),
                                                 ws=ws, lb=lb, ob=ob, cb=cb)

    elif (type(self.left) == Scalar and type(self.right) == Scalar):
      return "{ob}{0}*{1}{cb}".format(self.left.getString(i, j, k, pad, maxDepth),
                                      self.right.getString(i, j, k, pad, maxDepth),
                                      ob=ob, cb=cb)
    else:
      return "{ob}{0} * {1}{cb}".format(self.left.getString(i, j, k, pad, maxDepth),
                                        self.right.getString(i, j, k, pad, maxDepth),
                                        ob=ob, cb=cb)

class NodeStencilInterp(Node):
  def __init__(self, inner, dim):
    self.inner = inner
    self.depth = inner.depth + 1
    self.pad = 6

    self.dim = dim
    self.loc = copy.deepcopy(inner.loc)
    self.loc[dim] = not self.loc[dim]

  def getString(self, i, j, k, pad, maxDepth):
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

    maxDepth = max(self.depth, maxDepth)
    if (self.depth == maxDepth):
      pad = pad - 2
      ob = ''
      cb = ''
    else:
      ob = '('
      cb = ')'

    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "{ob}ci0*{0}\n{lb}{ws}+ ci1*{1}\n{lb}{ws}+ ci2*{2}\n{lb}{ws}+ ci3*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad, maxDepth),
          self.inner.getString(i1, j1, k1, pad, maxDepth),
          self.inner.getString(i2, j2, k2, pad, maxDepth),
          self.inner.getString(i3, j3, k3, pad, maxDepth),
          ws=ws, lb=lb, ob=ob, cb=cb)
    else:
      return "{ob}ci0*{0} + ci1*{1} + ci2*{2} + ci3*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad, maxDepth),
          self.inner.getString(i1, j1, k1, pad, maxDepth),
          self.inner.getString(i2, j2, k2, pad, maxDepth),
          self.inner.getString(i3, j3, k3, pad, maxDepth),
          ob=ob, cb=cb)

class NodeStencilGrad(Node):
  def __init__(self, inner, dim):
    self.inner = inner
    self.depth = inner.depth + 1
    self.pad = 6

    self.dim = dim
    self.loc = copy.deepcopy(inner.loc)
    self.loc[dim] = not self.loc[dim]

  def getString(self, i, j, k, pad, maxDepth):
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

    maxDepth = max(self.depth, maxDepth)
    if (self.depth == maxDepth):
      pad = pad - 2
      ob = ''
      cb = ''
    else:
      ob = '('
      cb = ')'

    if (self.depth > 1):
      ws = ''.rjust(pad)
      pad += self.pad

      lb = ''
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "{ob}cg0*{0}\n{lb}{ws}+ cg1*{1}\n{lb}{ws}+ cg2*{2}\n{lb}{ws}+ cg3*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad, maxDepth),
          self.inner.getString(i1, j1, k1, pad, maxDepth),
          self.inner.getString(i2, j2, k2, pad, maxDepth),
          self.inner.getString(i3, j3, k3, pad, maxDepth),
          ws=ws, lb=lb, ob=ob, cb=cb)
    else:
      return "{ob}cg0*{0} + cg1*{1} + cg2*{2} + cg3*{3}{cb}".format(
          self.inner.getString(i0, j0, k0, pad, maxDepth),
          self.inner.getString(i1, j1, k1, pad, maxDepth),
          self.inner.getString(i2, j2, k2, pad, maxDepth),
          self.inner.getString(i3, j3, k3, pad, maxDepth),
          ob=ob, cb=cb)

# Scalar class representing one grid cell
class Scalar(Node):
  def __init__(self, name, loc):
    self.name  = name
    self.depth = 0
    self.loc   = loc

  def getString(self, i, j, k, pad, maxDepth):
    if (i > 0):
      ii = "+{0}".format(i)
    elif (i < 0):
      ii = "{0}".format(i)
    else:
      ii = "  "
    if (j > 0):
      jj = "+{0}".format(j)
    elif (j < 0):
      jj = "{0}".format(j)
    else:
      jj = "  "
    if (k > 0):
      kk = "+{0}".format(k)
    elif (k < 0):
      kk = "{0}".format(k)
    else:
      kk = "  "

    return "{0}[i{1},j{2},k{3}]".format(self.name, ii, jj, kk)

# Define functions.
def interpx(inner):
  return NodeStencilInterp(inner, 0)
def interpy(inner):
  return NodeStencilInterp(inner, 1)
def interpz(inner):
  return NodeStencilInterp(inner, 2)

def gradx(inner):
  return NodeStencilGrad(inner, 0)
def grady(inner):
  return NodeStencilGrad(inner, 1)
def gradz(inner):
  return NodeStencilGrad(inner, 2)
