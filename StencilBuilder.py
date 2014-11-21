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
    self.pad = 8

    self.loc = copy.deepcopy(inner.loc)
    self.loc[dim] = not self.loc[dim]

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
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "{ob}ci0 * {0}\n{lb}{ws}+ ci1 * {1}\n{lb}{ws}+ ci2 * {2}\n{lb}{ws}+ ci3 * {3}{cb}".format(
          self.inner.getString(i-1 - self.loc[0], j, k, pad, maxDepth),
          self.inner.getString(i   - self.loc[0], j, k, pad, maxDepth),
          self.inner.getString(i+1 - self.loc[0], j, k, pad, maxDepth),
          self.inner.getString(i+2 - self.loc[0], j, k, pad, maxDepth),
          ws=ws, lb=lb, ob=ob, cb=cb)
    elif (type(self.inner) == Scalar):
      return "{ob}ci0*{0} + ci1*{1} + ci2*{2} + ci3*{3}{cb}".format(
          self.inner.getString(i-1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i   - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+2 - self.loc[0], j, k,  pad, maxDepth),
          ob=ob, cb=cb)
    else:
      return "{ob}ci0 * {0} + ci1 * {1} + ci2 * {2} + ci3 * {3}{cb}".format(
          self.inner.getString(i-1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i   - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+2 - self.loc[0], j, k,  pad, maxDepth),
          ob=ob, cb=cb)

class NodeStencilGrad(Node):
  def __init__(self, inner, dim):
    self.inner = inner
    self.depth = inner.depth + 1
    self.pad = 8

    self.loc = copy.deepcopy(inner.loc)
    self.loc[dim] = not self.loc[dim]

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
      for n in range(2, self.depth):
        lb = lb + '\n'
      return "{ob}cg0 * {0}\n{lb}{ws}+ cg1 * {1}\n{lb}{ws}+ cg2 * {2}\n{lb}{ws}+ cg3 * {3}{cb}".format(
          self.inner.getString(i-1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i   - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+2 - self.loc[0], j, k,  pad, maxDepth),
          ws=ws, lb=lb, ob=ob, cb=cb)
    elif (type(self.inner) == Scalar):
      return "{ob}cg0*{0} + cg1*{1} + cg2*{2} + cg3*{3}{cb}".format(
          self.inner.getString(i-1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i   - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+2 - self.loc[0], j, k,  pad, maxDepth),
          ob=ob, cb=cb)
    else:
      return "{ob}cg0 * {0} + cg1 * {1} + cg2 * {2} + cg3 * {3}{cb}".format(
          self.inner.getString(i-2 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i-1 - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i   - self.loc[0], j, k,  pad, maxDepth),
          self.inner.getString(i+1 - self.loc[0], j, k,  pad, maxDepth),
          ob=ob, cb=cb)


# Scalar class representing one grid cell
class Scalar(Node):
  def __init__(self, name, loc):
    self.name  = name
    self.depth = 0
    self.loc   = loc

  def getString(self, i, j, k, pad, maxDepth):
    if (i > 0):
      return "{0}[i+{1}]".format(self.name, i)
    elif (i < 0):
      return "{0}[i-{1}]".format(self.name, abs(i))
    else:
      return "{0}[i  ]".format(self.name, abs(i))

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
