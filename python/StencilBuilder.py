#!/usr/bin/python

import numpy as np

# Define the location arrays
uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
sloc = np.array([0,0,0])

zloc  = np.array([None, None, 0])
zhloc = np.array([None, None, 1])

# Check locations
def checkLocs(left, right):
    loc = np.array([ None, None, None ])
    for i in range(len(loc)):
        if (left.loc[i] == None and right.loc[i] == None):
            continue
        elif (left.loc[i] == None):
            loc[i] = right.loc[i]
        elif (right.loc[i] == None):
            loc[i] = left.loc[i]
        elif (left.loc[i] != right.loc[i]):
            raise RuntimeError("Left ({0}) and right ({0}) are not at the same grid location".format(
                type(left).__name__, type(right).__name__))
        else:
            loc[i] = left.loc[i]

    return loc

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
        self.left = left
        self.right = right
        self.operatorString = operatorString
        self.depth = max(left.depth, right.depth)
        self.depthk = max(left.depthk, right.depthk)

        if (self.depth > 1):
            self.pad = 2
        else:
            self.pad = 0

        self.loc = checkLocs(left, right)

    def getString(self, i, j, k, pad, plane, loc):
        ob = '( '
        cb = ' )'

        if (self.depth > 1):
            ws = ''.rjust(pad)
            pad += self.pad

            lb = ''
            for n in range(1, self.depth):
                lb = lb + '\n'

            return "{ob}{0}\n{lb}{ws}{os} {1}{cb}".format(self.left.getString(i, j, k, pad, plane, loc),
                                                          self.right.getString(i, j, k, pad, plane, loc),
                                                          ws=ws, lb=lb, ob=ob,
                                                          os=self.operatorString, cb=cb)
        else:
            return "{ob}{0} {os} {1}{cb}".format(self.left.getString(i, j, k, pad, plane, loc),
                                                 self.right.getString(i, j, k, pad, plane, loc),
                                                 ob=ob, os=self.operatorString, cb=cb)

class NodeStencilFour(Node):
    def __init__(self, inner, dim, c0, c1, c2, c3):
        self.inner = inner
        self.depth = inner.depth + 1
        self.depthk = inner.depthk + 1 if dim == 2 else inner.depth
        if (self.depthk > 2):
            raise RuntimeError("Type ({0}) exceeds maximum depth of 2".format(type(inner).__name__))
        self.pad = 6

        self.dim = dim
        self.loc = np.copy(inner.loc)
        self.loc[dim] = not self.loc[dim]

        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3

    def getString(self, i, j, k, pad, plane, loc):
        i0 = i1 = i2 = i3 = i
        j0 = j1 = j2 = j3 = j
        k0 = k1 = k2 = k3 = k

        bias = 0
        c0 = 'c' + self.c0[1:]
        c1 = 'c' + self.c1[1:]
        c2 = 'c' + self.c2[1:]
        c3 = 'c' + self.c3[1:]

        # Check in which cells biased schemes need to be applied.
        if (self.dim == 2):
            if ( ( loc == "bot"  and k == -1 ) or
                 ( loc == "both" and ( (self.loc[2] == 0 and k == -1) or 
                                       (self.loc[2] == 1 and k ==  0 and self.depthk == 2) ) ) or
                 ( loc == "bot+1h" and ( (self.loc[2] == 0 and k == -2) ) ) ):
                bias = 1
                c0 = 'b' + self.c0[1:]
                c1 = 'b' + self.c1[1:]
                c2 = 'b' + self.c2[1:]
                c3 = 'b' + self.c3[1:]
            elif ( ( loc == "top"  and k == 2 ) or
                   ( loc == "toph" and ( (self.loc[2] == 0 and k == 0) or
                                         (self.loc[2] == 1 and k == 0) ) ) or
                   ( loc == "top-1h" and ( (self.loc[2] == 0 and k == 1) ) ) ):
                bias = -1
                c0 = 't' + self.c0[1:]
                c1 = 't' + self.c1[1:]
                c2 = 't' + self.c2[1:]
                c3 = 't' + self.c3[1:]

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
            k0 += -1-self.loc[2] + bias
            k1 +=   -self.loc[2] + bias
            k2 += +1-self.loc[2] + bias
            k3 += +2-self.loc[2] + bias

        ob = '( '
        cb = ' )'

        newplane = np.copy(plane)
        if (self.depth-1 < 2):
            newplane[self.dim] = 1

        if (self.depth > 1):

            ws = ''.rjust(pad)
            pad += self.pad

            lb = ''
            for n in range(2, self.depth):
                lb = lb + '\n'
            return "{ob}{c0}*{0}\n{lb}{ws}+ {c1}*{1}\n{lb}{ws}+ {c2}*{2}\n{lb}{ws}+ {c3}*{3}{cb}".format(
                self.inner.getString(i0, j0, k0, pad, newplane, loc),
                self.inner.getString(i1, j1, k1, pad, newplane, loc),
                self.inner.getString(i2, j2, k2, pad, newplane, loc),
                self.inner.getString(i3, j3, k3, pad, newplane, loc),
                ws=ws, lb=lb, ob=ob, cb=cb,
                c0=c0, c1=c1, c2=c2, c3=c3)
        else:
            return "{ob}{c0}*{0} + {c1}*{1} + {c2}*{2} + {c3}*{3}{cb}".format(
                self.inner.getString(i0, j0, k0, pad, newplane, loc),
                self.inner.getString(i1, j1, k1, pad, newplane, loc),
                self.inner.getString(i2, j2, k2, pad, newplane, loc),
                self.inner.getString(i3, j3, k3, pad, newplane, loc),
                ob=ob, cb=cb,
                c0=c0, c1=c1, c2=c2, c3=c3)

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
        self.name = name
        self.depth = 0
        self.depthk = 0
        self.loc = np.copy(loc)

    def getString(self, i, j, k, pad, plane, loc):
        compact = True
        ii = jj = kk = ""
        if (plane[0] or not compact):
            ii = formatIndex(i, "ii")
        if (plane[1] or not compact):
            jj = formatIndex(j, "jj")
        if (plane[2] or not compact):
            kk = formatIndex(k, "kk")
        return "{0}[ijk{1}{2}{3}]".format(self.name, ii, jj, kk)

# Vector class representing a vertical profile
class Vector(Node):
    def __init__(self, name, loc):
        self.name = name
        self.depth = 0
        self.depthk = 0
        self.loc = np.copy(loc)

    def getString(self, i, j, k, pad, plane, loc):
        kk = formatIndex(k, "")
        # This is an unelegant solution but threats the double biased spatial operators
        if ( loc == "both" and (self.loc[2] == 1 and k == 0) ):
            return "{0}bot".format(self.name, kk)
        elif ( loc == "toph" and (self.loc[2] == 1 and k == 0) ):
            return "{0}top".format(self.name, kk)

        else:
            return "{0}[k{1}]".format(self.name, kk)

class Scalar(Node):
    def __init__(self, name):
        self.name = name
        self.depth = 0
        self.depthk = 0
        self.loc = np.array([ None, None, None ])

    def getString(self, i, j, k, pad, plane, loc):
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

def printStencil(lhs, rhs, operator, loc):
    checkLocs(lhs, rhs)

    # If the location is on the half level, add suffix h
    if (lhs.loc[2]):
        loc = loc + 'h'

    index = "[ijk]"
    indent = len(lhs.name) + len(index) + len(operator) + 2

    plane = np.array([0,0,0])
    print("{0}{1} {2} {3};".format(lhs.name, index, operator, rhs.getString(0, 0, 0, indent, plane, loc)))

def printLoop(lhs, rhs, operator, istart="grid->istart", iend="grid->iend",
                                  jstart="grid->jstart", jend="grid->jend",
                                  kstart="grid->kstart", kend="grid->kend"):
    checkLocs(lhs, rhs)

    index = "[ijk]"
    indexIndent = len(lhs.name) + len(index) + len(operator) + 2

    plane = np.array([0,0,0])

    # Set one indentation unit in C++ code
    indent = "    "

    print("{0}for (int k={1}; k<{2}; ++k)".format(indent * 0, kstart, kend))
    print("{0}for (int j={1}; j<{2}; ++j)".format(indent * 1, jstart, jend))
    print("{0}for (int i={1}; i<{2}; ++i)".format(indent * 2, istart, iend))
    print("{0}{{".format(indent * 2))
    print("{0}const int ijk = i + j*jj + k*kk;".format(indent * 3))
    print("{0}{1}{2} {3} {4};".format(indent * 3,
        lhs.name, index, operator, rhs.getString(0, 0, 0, indexIndent + len(indent * 3), plane, 0)))
    print("{0}}}".format(indent * 2))

def printEmptyLine(n=1):
    for i in range(n):
        print("")
