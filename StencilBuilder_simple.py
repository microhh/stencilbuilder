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
    if(self.depth > 1):
      self.pad = 2
    else:
      self.pad = 0

  def getString(self, i, pad, maxDepth):
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

      return "{ob}{0}\n{lb}{ws}+ {1}{cb}".format(self.left.getString(i, pad, maxDepth),
                                                 self.right.getString(i, pad, maxDepth),
                                                 ws=ws, lb=lb, ob=ob, cb=cb)

    else:
      return "{ob}{0} + {1}{cb}".format(self.left.getString(i, pad, maxDepth),
                                        self.right.getString(i, pad, maxDepth),
                                        ob=ob, cb=cb)

class NodeMult(Node):
  def __init__(self, left, right):
    self.left  = left
    self.right = right
    self.depth = max(left.depth, right.depth)
    if(self.depth > 1):
      self.pad = 2
    else:
      self.pad = 0

  def getString(self, i, pad, maxDepth):
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

      return "{ob}{0}\n{lb}{ws}* {1}{cb}".format(self.left.getString(i, pad, maxDepth),
                                                 self.right.getString(i, pad, maxDepth),
                                                 ws=ws, lb=lb, ob=ob, cb=cb)

    else:
      return "{ob}{0} * {1}{cb}".format(self.left.getString(i, pad, maxDepth),
                                        self.right.getString(i, pad, maxDepth),
                                        ob=ob, cb=cb)

class NodeStencil(Node):
  def __init__(self, inner):
    self.inner = inner
    self.depth = inner.depth + 1
    if (type(inner) == Scalar):
      self.pad = 6
    else:
      self.pad = 8

  def getString(self, i, pad, maxDepth):
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
          self.inner.getString(i-2, pad, maxDepth),
          self.inner.getString(i-1, pad, maxDepth),
          self.inner.getString(i  , pad, maxDepth),
          self.inner.getString(i+1, pad, maxDepth),
          ws=ws, lb=lb, ob=ob, cb=cb)
    else:
      return "{ob}ci0*{0} + ci1*{1} + ci2*{2} + ci3*{3}{cb}".format(
          self.inner.getString(i-2, pad, maxDepth),
          self.inner.getString(i-1, pad, maxDepth),
          self.inner.getString(i  , pad, maxDepth),
          self.inner.getString(i+1, pad, maxDepth),
          ob=ob, cb=cb)

# Scalar class representing one grid cell
class Scalar(Node):
  def __init__(self, data, name):
    self.data  = data
    self.name  = name
    self.depth = 0

  def getString(self, i, pad, maxDepth):
    if (i > 0):
      return "{0}[i+{1}]".format(self.name, i)
    elif (i < 0):
      return "{0}[i-{1}]".format(self.name, abs(i))
    else:
      return "{0}[i  ]".format(self.name, abs(i))

# Define functions.
def interp(inner):
  return NodeStencil(inner)
