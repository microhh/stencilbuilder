import numpy as np

# Define the location arrays
uloc = np.array([1,0,0])
vloc = np.array([0,1,0])
wloc = np.array([0,0,1])
sloc = np.array([0,0,0])

uwloc = np.array([1,0,1])
uvloc = np.array([1,1,0])
vwloc = np.array([0,1,1])

zloc  = np.array([None, None, 0])
zhloc = np.array([None, None, 1])

# Check whether float or int
def is_int_or_float(var):
    return ( isinstance(var, int) or isinstance(var, float) )

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

    def __pow__(self, right):
        return NodeOperatorPower(self, right)

    # Cover the case of multiplication with a scalar. Self does not exist and arguments are swapped
    __radd__ = __add__
    __rmul__ = __mul__

class NodeOperator(Node):
    def __init__(self, left, right, operatorString):
        # Swap the terms in case right is a scalar multiplication
        self.left  = left  if ( not is_int_or_float(right) ) else Scalar('{0}'.format(right))
        self.right = right if ( not is_int_or_float(right) ) else left
        self.operatorString = operatorString
        self.depth = max(self.left.depth, self.right.depth)
        self.depthk = max(self.left.depthk, self.right.depthk)

        if (self.depth > 1):
            self.pad = 2
        else:
            self.pad = 0

        self.loc = checkLocs(self.left, self.right)

    def getString(self, i, j, k, pad, plane, label, max_depthk, outer_depthk):
        max_depthk = max(max_depthk, self.depthk)

        ob = '( '
        cb = ' )'

        if (self.depth > 1):
            ws = ''.rjust(pad)
            pad += self.pad

            lb = ''
            for n in range(1, self.depth):
                lb = lb + '\n'

            return "{ob}{0}\n{lb}{ws}{os} {1}{cb}".format(self.left.getString(i, j, k, pad, plane, label, max_depthk, self.depthk),
                                                          self.right.getString(i, j, k, pad, plane, label, max_depthk, self.depthk),
                                                          ws=ws, lb=lb, ob=ob,
                                                          os=self.operatorString, cb=cb)
        else:
            return "{ob}{0} {os} {1}{cb}".format(self.left.getString(i, j, k, pad, plane, label, max_depthk, self.depthk),
                                                 self.right.getString(i, j, k, pad, plane, label, max_depthk, self.depthk),
                                                 ob=ob, os=self.operatorString, cb=cb)

class NodeOperatorPower(Node):
    def __init__(self, inner, power):
        self.inner = inner
        self.power = power
        self.depth = inner.depth
        self.depthk = inner.depthk

        if ( not is_int_or_float(self.power) ):
            raise RuntimeError("Only integer and float powers are supported")

        if (self.depth > 1):
            self.pad = 10
        else:
            self.pad = 0

        self.loc = inner.loc

    def getString(self, i, j, k, pad, plane, label, max_depthk, outer_depthk):
        max_depthk = max(max_depthk, self.depthk)

        ob = 'std::pow( '
        cb = ' )'

        if (self.depth > 1):
            ws = ''.rjust(pad)
            pad += self.pad

            lb = ''
            for n in range(1, self.depth):
                lb = lb + '\n'

            return "{ob}{0}\n{lb}{ws}, {1}{cb}".format(self.inner.getString(i, j, k, pad, plane, label, max_depthk, self.depthk),
                                                       self.power, ws=ws, lb=lb, ob=ob, cb=cb)
        else:
            return "{ob}{0}, {1}{cb}".format(self.inner.getString(i, j, k, pad, plane, label, max_depthk, self.depthk),
                                             self.power, ob=ob, cb=cb)

class NodeStencilFour(Node):
    def __init__(self, inner, dim, c0, c1, c2, c3):
        self.inner = inner
        self.depth = inner.depth + 1
        self.depthk = inner.depthk + 1 if dim == 2 else inner.depthk

        self.pad = 6

        self.dim = dim
        self.loc = np.copy(inner.loc)
        self.loc[dim] = not self.loc[dim]

        self.c0 = c0
        self.c1 = c1
        self.c2 = c2
        self.c3 = c3

    def getString(self, i, j, k, pad, plane, label, max_depthk, outer_depthk):
        max_depthk = max(max_depthk, self.depthk)

        i0 = i1 = i2 = i3 = i
        j0 = j1 = j2 = j3 = j
        k0 = k1 = k2 = k3 = k

        bias = 0
        c0 = self.c0[:] + "c"
        c1 = self.c1[:]
        c2 = self.c2[:]
        c3 = self.c3[:]

        # Check in which cells biased schemes need to be applied.
        if (self.dim == 2):
            # Check if the stencil has an odd depth, to take into account the staggering
            # If the stencil is evaluated at the cell faces, add one grid point, otherwise remove one
            if (max_depthk % 2 == 1):
                top_shift = 1 if ( (self.depthk % 2) != self.loc[2] ) else -1
            else:
                top_shift = 0

            # RULES:
            if ( ( label == "bot" and self.depthk == 1 and self.loc[2] == 0 and k == -1 ) or
                 ( label == "bot" and self.depthk == 1 and self.loc[2] == 1 and k == -1 ) or
                 ( label == "bot" and self.depthk == 2 and self.loc[2] == 0 and k == -1 ) or
                 ( label == "bot" and self.depthk == 2 and self.loc[2] == 1 and k ==  0 ) or
                 ( label == "bot" and self.depthk == 3 and self.loc[2] == 1 and k ==  0 ) or

                 ( label == "bot+1" and self.depthk == 1 and self.loc[2] == 0 and k == -2 ) or
                 ( label == "bot+1" and self.depthk == 1 and self.loc[2] == 1 and k == -2 ) or
                 ( label == "bot+1" and self.depthk == 2 and self.loc[2] == 0 and k == -2 ) or
                 ( label == "bot+1" and self.depthk == 2 and self.loc[2] == 1 and k == -1 ) or
                 ( label == "bot+1" and self.depthk == 3 and self.loc[2] == 1 and k == -1 ) or

                 ( label == "bot+2" and self.depthk == 1 and self.loc[2] == 0 and k == -3 ) or
                 ( label == "bot+2" and self.depthk == 1 and self.loc[2] == 1 and k == -3 ) or
                 ( label == "bot+2" and self.depthk == 2 and self.loc[2] == 0 and k == -3 ) or
                 ( label == "bot+2" and self.depthk == 2 and self.loc[2] == 1 and k == -2 ) or
                 ( label == "bot+2" and self.depthk == 3 and self.loc[2] == 1 and k == -2 ) ):

                bias = 1
                c0 = self.c0[:] + "b"
                c1 = self.c1[:]
                c2 = self.c2[:]
                c3 = self.c3[:]

            # RULES (More complex than bot, because of grid indexing):
            elif ( ( label == "top" and self.depthk == 1 and self.loc[2] == 0 and k == 0 + top_shift ) or
                   ( label == "top" and self.depthk == 2 and self.loc[2] == 1 and k == 0 + top_shift ) or

                   ( label == "top" and self.depthk == 1 and self.loc[2] == 1 and k == 2 + top_shift ) or
                   ( label == "top" and self.depthk == 2 and self.loc[2] == 0 and k == 1 + top_shift ) or
                   ( label == "top" and self.depthk == 3 and self.loc[2] == 1 and k == 1 + top_shift ) or

                   ( label == "top-1" and self.depthk == 1 and self.loc[2] == 0 and k == 1 + top_shift ) or
                   ( label == "top-1" and self.depthk == 2 and self.loc[2] == 1 and k == 1 + top_shift ) or

                   ( label == "top-1" and self.depthk == 1 and self.loc[2] == 1 and k == 3 + top_shift ) or
                   ( label == "top-1" and self.depthk == 2 and self.loc[2] == 0 and k == 2 + top_shift ) or
                   ( label == "top-1" and self.depthk == 3 and self.loc[2] == 1 and k == 2 + top_shift ) or

                   ( label == "top-2" and self.depthk == 1 and self.loc[2] == 0 and k == 2 + top_shift ) or
                   ( label == "top-2" and self.depthk == 2 and self.loc[2] == 1 and k == 2 + top_shift ) or

                   ( label == "top-2" and self.depthk == 1 and self.loc[2] == 1 and k == 4 + top_shift ) or
                   ( label == "top-2" and self.depthk == 2 and self.loc[2] == 0 and k == 3 + top_shift ) or
                   ( label == "top-2" and self.depthk == 3 and self.loc[2] == 1 and k == 3 + top_shift ) ):

                bias = -1
                c0 = self.c0[:] + "t"
                c1 = self.c1[:]
                c2 = self.c2[:]
                c3 = self.c3[:]

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
        if (self.depth-1 < 3):
            newplane[self.dim] = 1

        if (self.depth > 1):

            ws = ''.rjust(pad)
            pad += self.pad

            lb = ''
            for n in range(2, self.depth):
                lb = lb + '\n'
            return "{c0}{ob}{0} ,\n{lb}{ws}{1} ,\n{lb}{ws}{2} ,\n{lb}{ws}{3}{cb}".format(
                self.inner.getString(i0, j0, k0, pad, newplane, label, max_depthk, self.depthk),
                self.inner.getString(i1, j1, k1, pad, newplane, label, max_depthk, self.depthk),
                self.inner.getString(i2, j2, k2, pad, newplane, label, max_depthk, self.depthk),
                self.inner.getString(i3, j3, k3, pad, newplane, label, max_depthk, self.depthk),
                ws=ws, lb=lb, ob=ob, cb=cb,
                c0=c0, c1=c1, c2=c2, c3=c3)
        else:
            return "{c0}{ob}{0}, {1}, {2}, {3}{cb}".format(
                self.inner.getString(i0, j0, k0, pad, newplane, label, max_depthk, self.depthk),
                self.inner.getString(i1, j1, k1, pad, newplane, label, max_depthk, self.depthk),
                self.inner.getString(i2, j2, k2, pad, newplane, label, max_depthk, self.depthk),
                self.inner.getString(i3, j3, k3, pad, newplane, label, max_depthk, self.depthk),
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

    def getString(self, i, j, k, pad, plane, label, max_depthk, outer_depthk):
        max_depthk = max(max_depthk, self.depthk)

        compact = True
        ii = formatIndex(i, "ii") if ( plane[0] or not compact ) else ""
        jj = formatIndex(j, "jj") if ( plane[1] or not compact ) else ""
        kk = formatIndex(k, "kk") if ( plane[2] or not compact ) else ""
        return "{0}[ijk{1}{2}{3}]".format(self.name, ii, jj, kk)

# Vector class representing a vertical profile
class Vector(Node):
    def __init__(self, name, loc):
        self.name = name
        self.depth = 0
        self.depthk = 0
        self.loc = np.copy(loc)

    def getString(self, i, j, k, pad, plane, label, max_depthk, outer_depthk):
        compact = True
        kk = formatIndex(k, "") if ( plane[2] or not compact ) else ""

        # In case the vector is requested at the boundary in a stencil of a depth 
        if ( label == "bot" and self.loc[2] == 1 and k == 0 and outer_depthk >= 2):
            return "{0}bot".format(self.name)
        elif ( label == "top" and self.loc[2] == 1 and k == 0 and outer_depthk >= 2):
            return "{0}top".format(self.name)
        else:
            return "{0}[k{1}]".format(self.name, kk)

class Scalar(Node):
    def __init__(self, name):
        self.name = name
        self.depth = 0
        self.depthk = 0
        self.loc = np.array([ None, None, None ])

    def getString(self, i, j, k, pad, plane, label, max_depthk, outer_depthk):
        return "{0}".format(self.name)

# Define functions.
def interpx(inner):
    return NodeStencilFour(inner, 0, "interp4", "", "", "")
def interpy(inner):
    return NodeStencilFour(inner, 1, "interp4", "", "", "")
def interpz(inner):
    return NodeStencilFour(inner, 2, "interp4", "", "", "")

# Shortcuts for double interpolation
def interpxy(inner):
    return interpy(interpx(inner))
def interpxz(inner):
    return interpz(interpx(inner))
def interpyz(inner):
    return interpz(interpy(inner))
def interpxyz(inner):
    return interpz(interpy(interpx(inner)))

def gradx(inner):
    return NodeStencilFour(inner, 0, "grad4", "", "", "")
def grady(inner):
    return NodeStencilFour(inner, 1, "grad4", "", "", "")
def gradz(inner):
    return NodeStencilFour(inner, 2, "grad4", "", "", "")

def printEmptyLine(n=1):
    for i in range(n):
        print("")

def printStencil(lhs, rhs, operator, label, index="[ijk]"):
    checkLocs(lhs, rhs)

    indent = len(lhs.name) + len(index) + len(operator) + 2

    plane = np.array([0,0,0])
    max_depthk = 0

    print("Printing stencil {0} with vertical depth {1}, location {2}".format(lhs.name, rhs.depthk, label))
    printEmptyLine()

    printString = rhs.getString(0, 0, 0, indent, plane, label, max_depthk, 0)

    print("{0}{1} {2} {3};".format(lhs.name, index, operator, printString))

"""
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
"""
