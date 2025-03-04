"""Module provides functions to generate sparse operators for the chebyshev polynomials with a linear map y = ax + b."""

from __future__ import division
from __future__ import unicode_literals

import numpy as np
import scipy.sparse as spsp

import quicc.base.utils as utils

def linear_map(lower, upper):
    """Map lower and upper bound to y = ax + b"""

    a = (upper - lower)/2.0
    b = (upper + lower)/2.0

    return (a,b)

def d1(ny, lower, upper, zr = 1):
    """Create operator for 1st derivative"""

    a, b = linear_map(lower, upper)

    row = [2*j for j in range(0,ny)]
    mat = spsp.lil_matrix((ny,ny))
    for i in range(0,ny-1):
        mat[i,i+1:ny:2] = row[i+1:ny:2]
    mat[-zr:,:] = 0

    return (1.0/a)*mat

def d2(ny, lower, upper, zr = 2):
    """Create operator for 2nd derivative"""

    a, b = linear_map(lower, upper)

    mat = spsp.lil_matrix((ny,ny))
    for i in range(0,ny-2):
        mat[i,i+2:ny:2] = [j*(j**2 - i**2) for j in range(0,ny)][i+2:ny:2]
    mat[-zr:,:] = 0

    return (1.0/a**2)*mat

def d4(ny, lower, upper, zr = 4):
    """Create operator for 4th derivative"""

    mat_d2 = d2(ny + 4)
    mat = mat_d2*mat_d2
    mat = mat[0:-4, 0:-4]
    mat[-zr:,:] = 0

    return mat

def y1_diags(ny, lower, upper):
    """Create operator for y multiplication"""

    ns = np.arange(0, ny)
    offsets = np.arange(-1,2)
    nzrow = -1

    a, b = linear_map(lower, upper)

    # Generate 1st subdiagonal
    def d_1(n):
        return np.full(n.shape, a/2.0)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, b)

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n)

    ds = [d_1, d0, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def y2_diags(ny, lower, upper):
    """Create operator for y^2 multiplication"""

    ns = np.arange(0, ny)
    offsets = np.arange(-2,3)
    nzrow = -1

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return np.full(n.shape, a**2/4.0)

    # Generate 1st subdiagonal
    def d_1(n):
        return np.full(n.shape, a*b)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, (a**2 + 2.0*b**2)/2.0)

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n)

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def y4_diags(ny, lower, upper):
    """Create operator for y^4 multiplication"""

    ns = np.arange(0, ny)
    offsets = np.arange(-4,5)
    nzrow = -1

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return np.full(n.shape, a**4/16.0)

    # Generate 3rd subdiagonal
    def d_3(n):
        return np.full(n.shape, a**3*b/2.0)

    # Generate 2nd subdiagonal
    def d_2(n):
        return np.full(n.shape, a**2*(a**2 + 6.0*b**2)/4.0)

    # Generate 1st subdiagonal
    def d_1(n):
        return np.full(n.shape, a*b*(3.0*a**2 + 4.0*b**2)/2.0)

    # Generate diagonal
    def d0(n):
        return np.full(n.shape, (3.0*a**4 + 24.0*a**2*b**2 + 8.0*b**4)/8.0)

    # Generate 1st superdiagonal
    def d1(n):
        return  d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n)

    # Generate 3rd superdiagonal
    def d3(n):
        return d_3(n)

    # Generate 4th superdiagonal
    def d4(n):
        return d_4(n)

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i1_diags(ny, lower, upper):
    """Create operator for 1st integral f dy."""

    ns = np.arange(0, ny)
    offsets = np.arange(-1,2,2)
    nzrow = 0

    a, b = linear_map(lower, upper)

    # Generate 1st subdiagonal
    def d_1(n):
        return a/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    ds = [d_1, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i1y1_diags(ny, lower, upper):
    """Create operator for I y f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-2,3)
    nzrow = 0

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2/(4.0*n)

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b/(2.0*n)

    # Generate diagonal
    def d0(n):
        return 0

    # Generate 1st superdiagonal
    def d1(n):
        return -d_1(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return -d_2(n)

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2_diags(ny, lower, upper):
    """Create operator for I^2 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2/(4.0*n*(n - 1.0))

    # Generate diagonal
    def d0(n):
        return -a**2/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n+1.0)

    ds = [d_2, d0, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y2_diags(ny, lower, upper):
    """Create operator for I^2 y^2 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-4,5)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4/(16.0*n*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return (a**3*b)/(4.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (a**2*(2.0*b**2*n + a**2 + 2.0*b**2))/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(a**3*b)/(4.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -(a**2*(a**2 + 4.0*b**2))/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return d_1(n - 1.0)

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(a**2 - 2.0*b**2*n + 2.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return d_3(n + 1.0)

    # Generate 4th superdiagonal
    def d4(n):
        return d_4(n + 1.0)

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y3_diags(ny, lower, upper):
    """Create operator for I^2 y^3 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-5,6)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5/(32.0*n*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return 3.0*a**4*b/(16.0*n*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(a**2*n + 3.0*a**2 + 12.0*b**2*n + 12.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(3.0*a**2 + 2.0*b**2*n + 2.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3*(a**2 + 6.0*b**2)/(16.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*b*(3.0*a**2 + 4.0*b**2)/(8.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a**3*(a**2 + 6.0*b**2)/(16.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(3.0*a**2 - 2.0*b**2*n + 2.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*(a**2*n - 3.0*a**2 + 12.0*b**2*n - 12.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return 3.0*a**4*b/(16.0*n*(n + 1.0))

    # Generate 5th superdiagonal
    def d5(n):
        return a**5/(32.0*n*(n + 1.0))

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2lapl_3d_diags(ny, lower, upper, k, l):
    """Create operator for I^2 Lapl3D f(y) dy for 1 Chebyshev direction"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(k**2 + l**2)/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (a**2*(k**2 + l**2) + 2.0*n**2 - 2.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(k**2 + l**2)/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2lapl_2d_diags(ny, lower, upper, k):
    """Create operator for I2 Lapl2D f(y) dy for 1 chebyshev direction"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-2,3,2)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*k**2/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return (a**2*k**2 + 2.0*n**2 - 2.0)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*k**2/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y2sphlapl_diags(ny, lower, upper, l):
    """Create operator for I^2 y^2 SphLapl(f(y)) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-2,3)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(a**2*(l - n + 2.0)*(l + n - 1.0))/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (a*b*(n - 1.0))/n

    # Generate main diagonal
    def d0(n):
        return (a**2*l**2 + a**2*l + a**2*n**2 - a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(n + 1.0)/n

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(l - n - 1.0)*(l + n + 2.0)/(4.0*n*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y3sphlapl_diags(ny, lower, upper, l):
    """Create operator for I^2 y^3 SphLapl(f(y)) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-3,4)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 3rd subdiagonal
    def d_3(n):
        return -a**3*(l - n + 3.0)*(l + n - 2.0)/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*b*(l**2 + l - 3.0*n**2 + 11.0*n - 10.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**2*l**2 + a**2*l + 3.0*a**2*n**2 - a**2*n - 6.0*a**2 + 12.0*b**2*n**2 - 4.0*b**2*n - 16.0*b**2)/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return b*(a**2*l**2 + a**2*l + 3.0*a**2*n**2 - 5.0*a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*(a**2*l**2 + a**2*l + 3.0*a**2*n**2 + a**2*n - 6.0*a**2 + 12.0*b**2*n**2 + 4.0*b**2*n - 16.0*b**2)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(l**2 + l - 3.0*n**2 - 11.0*n - 10.0)/(4.0*n*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*(l - n - 2.0)*(l + n + 3.0)/(8.0*n*(n + 1.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i3_diags(ny, lower, upper):
    """Create operator for I^3 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-3,4,2)
    nzrow = 2

    a, b = linear_map(lower, upper)

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3/(8.0*n*(n - 2.0)*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0*a**3/(8.0*n*(n - 2.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*a**3/(8.0*n*(n - 1.0)*(n + 2.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3/(8.0*n*(n + 1.0)*(n + 2.0))

    ds = [d_3, d_1, d1, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4_diags(ny, lower, upper):
    """Create operator for I^4 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate diagonal
    def d0(n):
        return 3.0*a**4/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**4/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4d1_diags(ny, lower, upper):
    """Create operator for I^4 D f(y) dy"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-3,4,2)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 3rd subdiagonal
    def d_2(n):
        return a**3/(8.0*n*(n - 2.0)*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3*3.0/(8.0*n*(n - 2.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a**3*3.0/(8.0*n*(n - 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**3/(8.0*n*(n + 1.0)*(n + 2.0))

    ds = [d_2, d_1, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4d2_diags(ny, lower, upper):
    """Create operator for I^4 D^2 f(y) dy"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-2,3,2)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2/(4.0*n*(n - 1.0))

    # Generate main diagonal
    def d0(n):
        return -a**2/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2/(4.0*n*(n + 1.0))

    ds = [d_2, d0, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y1_diags(ny, lower, upper):
    """Create operator for I^4 y f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-5,6)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*b/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -3.0*a**5/(32.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*b/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a**5*(n - 8.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return 3.0*a**4*b/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a**5*(n + 8.0)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**4*b/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -3.0*a**5/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*b/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return a**5/(32.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y1d1y1_diags(ny, lower, upper):
    """Create operator for I^4 y D(f(y) y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-5,6)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5*(n - 4.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*b*(2.0*n - 7.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -a**3*(a**2*n - 8.0*a**2 - 4.0*b**2*n - 4.0*b**2)/(32.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*b*(n - 4.0)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3*(a**2*n**2 + 2.0*a**2*n - 20.0*a**2 + 6.0*b**2*n**2 - 6.0*b**2*n - 36.0*b**2)/(16.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return -9.0*a**4*b/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a**3*(a**2*n**2 - 2.0*a**2*n - 20.0*a**2 + 6.0*b**2*n**2 + 6.0*b**2*n - 36.0*b**2)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**4*b*(n + 4.0)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*(a**2*n + 8.0*a**2 - 4.0*b**2*n + 4.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*b*(2.0*n + 7.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*(n + 4.0)/(32.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y3d2_diags(ny, lower, upper):
    """Create operator for I^4 y^3 D^2 f(y) dy."""

    ns = np.arange(0, ny)
    offsets = np.arange(-5,6)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5*(n - 6.0)*(n - 5.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return 3.0*a**4*b*(n - 5.0)*(n - 4.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(a**2*n**2 + 7.0*a**2*n - 54.0*a**2 + 12.0*b**2*n**2 - 36.0*b**2*n - 48.0*b**2)/(32.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(15.0*a**2*n - 57.0*a**2 + 2.0*b**2*n**2 - 4.0*b**2*n - 6.0*b**2)/(8.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a**3*(a**2*n**3 - 9.0*a**2*n**2 - 16.0*a**2*n + 132.0*a**2 + 6.0*b**2*n**3 - 54.0*b**2*n**2 + 12.0*b**2*n + 288.0*b**2)/(16.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*b*(3.0*a**2*n**2 - 66.0*a**2 + 4.0*b**2*n**2 - 16.0*b**2)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a**3*(a**2*n**3 + 9.0*a**2*n**2 - 16.0*a**2*n - 132.0*a**2 + 6.0*b**2*n**3 + 54.0*b**2*n**2 + 12.0*b**2*n - 288.0*b**2)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(15.0*a**2*n + 57.0*a**2 - 2.0*b**2*n**2 - 4.0*b**2*n + 6.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*(a**2*n**2 - 7.0*a**2*n - 54.0*a**2 + 12.0*b**2*n**2 + 36.0*b**2*n - 48.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return 3.0*a**4*b*(n + 4.0)*(n + 5.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return a**5*(n + 5.0)*(n + 6.0)/(32.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y4_diags(ny, lower, upper):
    """Create operator for I^4 y^4 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-8,9)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 8th subdiagonal
    def d_8(n):
        return a**8/(256.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0))

    # Generate 7th subdiagonal
    def d_7(n):
        return (a**7*b)/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return (3.0*a**6*(2.0*b**2*n + a**2 + 2.0*b**2))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -(a**5*b*(a**2*n - 4.0*b**2*n - 11.0*a**2 - 4.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -(a**4*(a**4*n**2 - 19.0*a**4 + 12.0*a**2*b**2*n**2 - 36.0*a**2*b**2*n - 120.0*a**2*b**2 - 4.0*b**4*n**2 - 12.0*b**4*n - 8.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(3.0*a**5*b*(a**2*n + 4.0*b**2*n + 6.0*a**2 + 8.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(a**4*(9.0*a**4*n + 33.0*a**4 + 6.0*a**2*b**2*n**2 + 120*a**2*b**2*n + 306.0*a**2*b**2 + 16.0*b**4*n**2 + 80.0*b**4*n + 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (a**5*b*(3.0*a**2*n**2 - 15.0*a**2*n - 102.0*a**2 + 8.0*b**2*n**2 - 40.0*b**2*n - 192.0*b**2))/(32.0*n*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*(a**4*n**2 - 29.0*a**4 + 16.0*a**2*b**2*n**2 - 304.0*a**2*b**2 + 16.0*b**4*n**2 - 144.0*b**4))/(128.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (a**5*b*(3.0*a**2*n**2 + 15.0*a**2*n - 102.0*a**2 + 8.0*b**2*n**2 + 40.0*b**2*n - 192.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (a**4*(9.0*a**4*n - 33.0*a**4 - 6.0*a**2*b**2*n**2 + 120.0*a**2*b**2*n - 306.0*a**2*b**2 - 16.0*b**4*n**2 + 80.0*b**4*n - 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(3.0*a**5*b*(a**2*n + 4.0*b**2*n - 6.0*a**2 - 8.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -(a**4*(a**4*n**2 - 19.0*a**4 + 12.0*a**2*b**2*n**2 + 36.0*a**2*b**2*n - 120.0*a**2*b**2 - 4.0*b**4*n**2 + 12.0*b**4*n - 8.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -(a**5*b*(a**2*n - 4.0*b**2*n + 11.0*a**2 + 4.0*b**2))/(32.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -(3.0*a**6*(a**2 - 2.0*b**2*n + 2.0*b**2))/(64.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 7th superdiagonal
    def d7(n):
        return d_7(n + 3.0)

    # Generate 8th superdiagonal
    def d8(n):
        return d_8(n + 3.0)

    ds = [d_8, d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7, d8]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y4sphlaply1d1y1_diags(ny, lower, upper):
    """Create operator for I^4 y^4 SphLapl(y D(f(y) y)) dy"""

    raise RuntimeError("This doesn't look right!")

    ns = np.arange(0, ny)
    offsets = np.arange(-5,6)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5*(n - 6.0)*(n - 5.0)*(n - 4.0)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*b*(n - 5.0)*(n - 4.0)*(4.0*n - 15.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return 3.0*a**3*(n - 4.0)*(a**2*n**2 - a**2*n - 14.0*a**2 + 8.0*b**2*n**2 - 20.0*b**2*n - 28.0*b**2)/(32.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(4.0*a**2*n**3 - 12.0*a**2*n**2 - 67.0*a**2*n + 213.0*a**2 + 8.0*b**2*n**3 - 42.0*b**2*n**2 + 28.0*b**2*n + 78.0*b**2)/(8.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**4*n**4 + 7.0*a**4*n**3 - 52.0*a**4*n**2 - 52.0*a**4*n + 384.0*a**4 + 12.0*a**2*b**2*n**4 + 30.0*a**2*b**2*n**3 - 354.0*a**2*b**2*n**2 - 12.0*a**2*b**2*n + 1440.0*a**2*b**2 + 8.0*b**4*n**4 - 16.0*b**4*n**3 - 56.0*b**4*n**2 + 64.0*b**4*n + 96.0*b**4)/(16.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return 9.0*a**2*b*(3.0*a**2*n**2 - 26.0*a**2 + 4.0*b**2*n**2 - 16.0*b**2)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*(a**4*n**4 - 7.0*a**4*n**3 - 52.0*a**4*n**2 + 52.0*a**4*n + 384.0*a**4 + 12.0*a**2*b**2*n**4 - 30.0*a**2*b**2*n**3 - 354.0*a**2*b**2*n**2 + 12.0*a**2*b**2*n + 1440.0*a**2*b**2 + 8.0*b**4*n**4 + 16.0*b**4*n**3 - 56.0*b**4*n**2 - 64.0*b**4*n + 96.0*b**4)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(4.0*a**2*n**3 + 12.0*a**2*n**2 - 67.0*a**2*n - 213.0*a**2 + 8.0*b**2*n**3 + 42.0*b**2*n**2 + 28.0*b**2*n - 78.0*b**2)/(8.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -3.0*a**3*(n + 4.0)*(a**2*n**2 + a**2*n - 14.0*a**2 + 8.0*b**2*n**2 + 20.0*b**2*n - 28.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*b*(n + 4.0)*(n + 5.0)*(4.0*n + 15.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*(n + 4.0)*(n + 5.0)*(n + 6.0)/(32.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y4sphlapl_diags(ny, lower, upper, l):
    """Create operator for I^4 y^4 SphLapl f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-6,7)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 6th subdiagonal
    def d_6(n):
        return -(a**6*(l - n + 6.0)*(l + n - 5.0))/(64.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -(a**5*b*(l**2 + l - 2.0*n**2 + 19.0*n - 45.0))/(16.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return (a**4*(a**2*l**2*n - 5.0*a**2*l**2 + a**2*l*n - 5.0*a**2*l + a**2*n**3 - 3.0*a**2*n**2 - 28.0*a**2*n + 96.0*a**2 - 2.0*b**2*l**2*n - 2.0*b**2*l**2 - 2.0*b**2*l*n - 2.0*b**2*l + 12.0*b**2*n**3 - 84.0*b**2*n**2 + 96.0*b**2*n + 192.0*b**2))/(32.0*n*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return (a**3*b*(3.0*a**2*l**2 + 3.0*a**2*l + 2.0*a**2*n**2 + 11.0*a**2*n - 75.0*a**2 + 8.0*b**2*n**2 - 20.0*b**2*n - 28.0*b**2))/(16.0*n*(n - 1.0)*(n - 2.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return (a**2*(a**4*l**2*n + 17.0*a**4*l**2 + a**4*l*n + 17.0*a**4*l - a**4*n**3 + 24.0*a**4*n**2 - 5.0*a**4*n - 294.0*a**4 + 16.0*a**2*b**2*l**2*n + 32.0*a**2*b**2*l**2 + 16.0*a**2*b**2*l*n + 32.0*a**2*b**2*l + 192.0*a**2*b**2*n**2 - 288.0*a**2*b**2*n - 1344.0*a**2*b**2 + 16.0*b**4*n**3 - 112.0*b**4*n - 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -(a**3*b*(a**2*l**2*n - 8.0*a**2*l**2 + a**2*l*n - 8.0*a**2*l + 2.0*a**2*n**3 - 15.0*a**2*n**2 - 23.0*a**2*n + 180.0*a**2 + 4.0*b**2*n**3 - 30.0*b**2*n**2 + 2.0*b**2*n + 156.0*b**2))/(8.0*n*(n - 2.0)*(n - 3.0)*(n + 2.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -(a**2*(a**4*l**2*n**2 - 19.0*a**4*l**2 + a**4*l*n**2 - 19.0*a**4*l + a**4*n**4 - 37.0*a**4*n**2 + 312.0*a**4 + 6.0*a**2*b**2*l**2*n**2 - 54.0*a**2*b**2*l**2 + 6.0*a**2*b**2*l*n**2 - 54.0*a**2*b**2*l + 12.0*a**2*b**2*n**4 - 300.0*a**2*b**2*n**2 + 1728.0*a**2*b**2 + 8.0*b**4*n**4 - 104.0*b**4*n**2 + 288.0*b**4))/(16.0*(n - 1.0)*(n - 2.0)*(n - 3.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -(a**3*b*(a**2*l**2*n + 8.0*a**2*l**2 + a**2*l*n + 8.0*a**2*l + 2.0*a**2*n**3 + 15.0*a**2*n**2 - 23.0*a**2*n - 180.0*a**2 + 4.0*b**2*n**3 + 30.0*b**2*n**2 + 2.0*b**2*n - 156.0*b**2))/(8.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return (a**2*(a**4*l**2*n - 17.0*a**4*l**2 + a**4*l*n - 17.0*a**4*l - a**4*n**3 - 24.0*a**4*n**2 - 5.0*a**4*n + 294.0*a**4 + 16.0*a**2*b**2*l**2*n - 32.0*a**2*b**2*l**2 + 16.0*a**2*b**2*l*n - 32.0*a**2*b**2*l - 192.0*a**2*b**2*n**2 - 288.0*a**2*b**2*n + 1344.0*a**2*b**2 + 16.0*b**4*n**3 - 112.0*b**4*n + 96.0*b**4))/(64.0*n*(n - 1.0)*(n - 2.0)*(n + 3.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return (a**3*b*(3.0*a**2*l**2 + 3.0*a**2*l + 2.0*a**2*n**2 - 11.0*a**2*n - 75.0*a**2 + 8.0*b**2*n**2 + 20.0*b**2*n - 28.0*b**2))/(16.0*n*(n - 1.0)*(n + 2.0)*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return (a**4*(a**2*l**2*n + 5.0*a**2*l**2 + a**2*l*n + 5.0*a**2*l + a**2*n**3 + 3.0*a**2*n**2 - 28.0*a**2*n - 96.0*a**2 - 2.0*b**2*l**2*n + 2.0*b**2*l**2 - 2.0*b**2*l*n + 2.0*b**2*l + 12.0*b**2*n**3 + 84.0*b**2*n**2 + 96.0*b**2*n - 192.0*b**2))/(32.0*n*(n - 1.0)*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -(a**5*b*(l**2 + l - 2.0*n**2 - 19.0*n - 45.0))/(16.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -(a**6*(l - n - 5.0)*(l + n + 6.0))/(64.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4lapl_3d_diags(ny, lower, upper, k, l):
    """Create operator for I^4 Lapl3D f(y) dy with 1 chebyshev direction"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return -a**4*(k**2 + l**2)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(a**2*(k**2 + l**2) + n**2 - 2.0*n - 3.0)/(4*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*(3.0*a**2*(k**2 + l**2) + 4.0*n**2 - 16.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2*(a**2*(k**2 + l**2) + n**2 + 2.0*n - 3.0)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*(k**2 + l**2)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4lapl_2d_diags(ny, lower, upper, k):
    """Create operator for I^4 Lapl2D f(y) dy with 1 Chebyshev direction"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return -a**4*k**2/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(a**2*k**2 + n**2 - 2.0*n - 3.0)/(4*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return -a**2*(3.0*a**2*k**2 + 4.0*n**2 - 16.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2*(a**2*k**2 + n**2 + 2.0*n - 3.0)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*k**2/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4lapl2_3d_diags(ny, lower, upper, k, l):
    """Create operator for I^4 Lapl3D^2 f(y) dy with 1 chebyshev direction"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*(k**2 + l**2)**2/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*(k**2 + l**2)*(a**2*(k**2 + l**2) + 2.0*n**2 - 4.0*n - 6.0)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*(k**4 + 2.0*k**2*l**2 + l**4) + 8.0*a**2*(k**2*n**2 - 4.0*k**2 + l**2*n**2 - 4.0*l**2) + 8.0*n**4 - 40.0*n**2 + 32.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(k**2 + l**2)*(a**2*(k**2 + l**2) + 2.0*n**2 + 4.0*n - 6.0)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*(k**2 + l**2)**2/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4lapl2_2d_diags(ny, lower, upper, k):
    """Create operator for I^4 Lapl2D^2 f(y) dy  with 1 Chebyshev direction"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-4,5,2)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*k**4/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**2*k**2*(a**2*k**2 + 2.0*n**2 - 4.0*n - 6.0)/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*k**4 + 8.0*a**2*(k**2*n**2 - 4.0*k**2) + 8.0*n**4 - 40.0*n**2 + 32.0)/(8.0*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*k**2*(a**2*k**2 + 2.0*n**2 + 4.0*n - 6.0)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*k**4/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_4, d_2, d0, d2, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y2sphlapl2_l1_diags(ny, lower, upper):
    """Create operator for I^4 y^2 SphLapl^2 f(y) dy for l = 1."""

    ns = np.arange(0, ny)
    offsets = np.arange(-2,3)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(n - 5.0)/(4.0*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*b*(n - 2.0)/n

    # Generate main diagonal
    def d0(n):
        return (a**2*n**2 + 3.0*a**2 + 2.0*b**2*n**2 - 2.0*b**2)/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a*b*(n + 2.0)/n

    # Generate 2nd superdiagonal
    def d2(n):
        return a**2*(n + 5.0)/(4.0*(n + 1.0))

    ds = [d_2, d_1, d0, d1, d2]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y4sphlapl2_diags(ny, lower, upper, l):
    """Create operator for I^4 y^4 SphLapl^2 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-4,5)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 4th subdiagonal
    def d_4(n):
        return (a**4*(l - n + 6.0)*(l + n - 5.0)*(l - n + 4.0)*(l + n - 3.0))/(16.0*n*(n - 3.0)*(n - 1.0)*(n - 2.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -(a**3*b*(n - 4.0)*(l**2 + l - n**2 + 8.0*n - 15.0))/(2.0*n*(n - 1.0)*(n - 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -(a**2*(a**2*l**4 + 2.0*a**2*l**3 + 5.0*a**2*l**2*n - 20.0*a**2*l**2 + 5.0*a**2*l*n - 21.0*a**2*l - a**2*n**4 + 9.0*a**2*n**3 - 17.0*a**2*n**2 - 39.0*a**2*n + 108.0*a**2 + 2.0*b**2*l**2*n**2 - 4.0*b**2*l**2*n - 6.0*b**2*l**2 + 2.0*b**2*l*n**2 - 4.0*b**2*l*n - 6.0*b**2*l - 6.0*b**2*n**4 + 54.0*b**2*n**3 - 138.0*b**2*n**2 + 18.0*b**2*n + 216.0*b**2))/(4.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return (a*b*(a**2*l**2*n - 8.0*a**2*l**2 + a**2*l*n - 8.0*a**2*l + 3.0*a**2*n**3 - 12.0*a**2*n**2 - 15.0*a**2*n + 72.0*a**2 + 4.0*b**2*n**3 - 16.0*b**2*n**2 + 4.0*b**2*n + 24.0*b**2))/(2.0*n*(n + 1.0)*(n - 2.0))

    # Generate main diagonal
    def d0(n):
        return (3.0*a**4*l**4 + 6.0*a**4*l**3 + 2.0*a**4*l**2*n**2 - 47.0*a**4*l**2 + 2.0*a**4*l*n**2 - 50.0*a**4*l + 3.0*a**4*n**4 - 51.0*a**4*n**2 + 228.0*a**4 + 8.0*a**2*b**2*l**2*n**2 - 32.0*a**2*b**2*l**2 + 8.0*a**2*b**2*l*n**2 - 32.0*a**2*b**2*l + 24.0*a**2*b**2*n**4 - 264.0*a**2*b**2*n**2 + 672.0*a**2*b**2 + 8.0*b**4*n**4 - 40.0*b**4*n**2 + 32.0*b**4)/(8.0*(n - 1.0)*(n - 2.0)*(n + 2.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return (a*b*(a**2*l**2*n + 8.0*a**2*l**2 + a**2*l*n + 8.0*a**2*l + 3.0*a**2*n**3 + 12.0*a**2*n**2 - 15.0*a**2*n - 72.0*a**2 + 4.0*b**2*n**3 + 16.0*b**2*n**2 + 4.0*b**2*n - 24.0*b**2))/(2.0*n*(n + 2.0)*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -(a**2*(a**2*l**4 + 2.0*a**2*l**3 - 5.0*a**2*l**2*n - 20.0*a**2*l**2 - 5.0*a**2*l*n - 21.0*a**2*l - a**2*n**4 - 9.0*a**2*n**3 - 17.0*a**2*n**2 + 39.0*a**2*n + 108.0*a**2 + 2.0*b**2*l**2*n**2 + 4.0*b**2*l**2*n - 6.0*b**2*l**2 + 2.0*b**2*l*n**2 + 4.0*b**2*l*n - 6.0*b**2*l - 6.0*b**2*n**4 - 54.0*b**2*n**3 - 138.0*b**2*n**2 - 18.0*b**2*n + 216.0*b**2))/(4.0*n*(n + 3.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -(a**3*b*(n + 4.0)*(l**2 + l - n**2 - 8.0*n - 15.0))/(2.0*n*(n + 2.0)*(n + 1.0))

    # Generate 4th superdiagonal
    def d4(n):
        return (a**4*(l + n + 6.0)*(l - n - 3.0)*(l + n + 4.0)*(l - n - 5.0))/(16.0*n*(n + 3.0)*(n + 2.0)*(n + 1.0))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2d1_diags(ny, lower, upper):
    """Create operator for I2 D f(y) dy"""

    ns = np.arange(0, ny, 1)
    offsets = np.arange(-1,2,2)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 1st subdiagonal
    def d_1(n):
        return a/(2.0*n)

    # Generate 1st superdiagonal
    def d1(n):
        return -a/(2.0*n)

    ds = [d_1, d1]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y1_diags(ny, lower, upper):
    """Create operator for I^2 y f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-3,4)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -d_3(n + 1.0)

    # Generate main diagonal
    def d0(n):
        return -a**2*b/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -d_3(n)

    # Generate 2nd superdiagonal
    def d2(n):
        return d_2(n + 1.0)

    # Generate 3rd superdiagonal
    def d3(n):
        return d_3(n + 1.0)

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y1d1y1_diags(ny, lower, upper):
    """Create operator for I^2 y D(y f(y)) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-3,4)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(n - 2.0)/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(2.0*n - 3.0)/(4.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**2*n + 2.0*a**2 + 4.0*b**2*n + 4.0*b**2)/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return a**2*b/(2.0*(n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*(a**2*n - 2.0*a**2 + 4.0*b**2*n - 4.0*b**2)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(2.0*n + 3.0)/(4.0*n*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*(n + 2.0)/(8.0*n*(n + 1.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y2d1y1_diags(ny, lower, upper):
    """Create operator for I^2 y^2 D(y f(y)) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-3,4)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 4rd subdiagonal
    def d_4(n):
        return a**4*(n - 3)/(16*n*(n - 1))

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*b*(3*n - 7)/(8*n*(n - 1))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*(a**2*n**2 - 3*a**2 + 6*b**2*n**2 - 4*b**2*n - 10*b**2)/(8*n*(n - 1)*(n + 1))

    # Generate 1st subdiagonal
    def d_1(n):
        return -a*b*(3*a**2*n - 7*a**2 + 4*b**2*n - 4*b**2)/(8*n*(n - 1))

    # Generate main diagonal
    def d0(n):
        return a**2*(a**2 + 4*b**2)/(4*(n - 1)*(n + 1))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*b*(3*a**2*n - 7*a**2 + 4*b**2*n - 4*b**2)/(8*n*(n - 1))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*(a**2*n**2 - 3*a**2 + 6*b**2*n**2 + 4*b**2*n - 10*b**2)/(8*n*(n - 1)*(n + 1))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*b*(3*n + 7)/(8*n*(n + 1))
    
    # Generate 4rd superdiagonal
    def d4(n):
        return -a**4*(n + 3)/(16*n*(n + 1))

    ds = [d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i2y2d1_diags(ny, lower, upper):
    """Create operator for I^2 y^2 D f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-3,4)
    nzrow = 1

    a, b = linear_map(lower, upper)

    # Generate 3rd subdiagonal
    def d_3(n):
        return a**3*(n - 3.0)/(8.0*n*(n - 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return a**2*b*(n - 2.0)/(2.0*n*(n - 1.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a*(a**2*n + 3.0*a**2 + 4.0*b**2*n + 4.0*b**2)/(8.0*n*(n + 1.0))

    # Generate main diagonal
    def d0(n):
        return a**2*b/((n - 1.0)*(n + 1.0))

    # Generate 1st superdiagonal
    def d1(n):
        return -a*(a**2*n - 3.0*a**2 + 4.0*b**2*n - 4.0*b**2)/(8.0*n*(n - 1.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**2*b*(n + 2.0)/(2.0*n*(n + 1.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -a**3*(n + 3.0)/(8.0*n*(n + 1.0))

    ds = [d_3, d_2, d_1, d0, d1, d2, d3]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y2_diags(ny, lower, upper):
    """Create operator for I^4 y^2 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-6,7)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 6th subdiagonal
    def d_6(n):
        return a**6/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5*b/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -a**4*(a**2*n - 5.0*a**2 - 2.0*b**2*n - 2.0*b**2)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -3.0*a**5*b/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*(a**2*n + 17.0*a**2 + 16.0*b**2*n + 32.0*b**2)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return a**5*b*(n - 8.0)/(8.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0))

    # Generate main diagonal
    def d0(n):
        return a**4*(a**2*n**2 - 19.0*a**2 + 6.0*b**2*n**2 - 54.0*b**2)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return a**5*b*(n + 8.0)/(8.0*n*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**4*(a**2*n - 17.0*a**2 + 16.0*b**2*n - 32.0*b**2)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -3*a**5*b/(16.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*(a**2*n + 5.0*a**2 - 2.0*b**2*n + 2.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return a**5*b/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return a**6/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y3_diags(ny, lower, upper):
    """Create operator for I^4 y^3 f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-7,8)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 7th subdiagonal
    def d_7(n):
        return a**7/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return 3.0*a**6*b/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return -a**5*(a**2*n - 11.0*a**2 - 12.0*b**2*n - 12.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return -a**4*b*(3.0*a**2*n - 15.0*a**2 - 2.0*b**2*n - 2.0*b**2)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -3.0*a**5*(a**2*n + 6.0*a**2 + 12.0*b**2*n + 24.0*b**2)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*b*(3.0*a**2*n + 51.0*a**2 + 16.0*b**2*n + 32.0*b**2)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return 3.0*a**5*(a**2*n**2 - 5.0*a**2*n - 34.0*a**2 + 8.0*b**2*n**2 - 40.0*b**2*n - 192.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate main diagonal
    def d0(n):
        return 3.0*a**4*b*(a**2*n**2 - 19.0*a**2 + 2.0*b**2*n**2 - 18.0*b**2)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*a**5*(a**2*n**2 + 5.0*a**2*n - 34.0*a**2 + 8.0*b**2*n**2 + 40.0*b**2*n - 192.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return -a**4*b*(3.0*a**2*n - 51.0*a**2 + 16.0*b**2*n - 32.0*b**2)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return -3.0*a**5*(a**2*n - 6.0*a**2 + 12.0*b**2*n - 24.0*b**2)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return -a**4*b*(3.0*a**2*n + 15.0*a**2 - 2.0*b**2*n + 2.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*(a**2*n + 11.0*a**2 - 12.0*b**2*n + 12.0*b**2)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return 3.0*a**6*b/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 7th superdiagonal
    def d7(n):
        return a**7/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y3d1y1_diags(ny, lower, upper):
    """Create operator for I^4 y^3 D(f(y) y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-7,8)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 7th subdiagonal
    def d_7(n):
        return a**7*(n - 6.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return a**6*b*(4.0*n - 21.0)/(64.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5*(a**2*n**2 + 7.0*a**2*n - 54.0*a**2 + 24.0*b**2*n**2 - 84.0*b**2*n - 108.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*b*(21.0*a**2*n - 81.0*a**2 + 8.0*b**2*n**2 - 22.0*b**2*n - 30.0*b**2)/(32.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -a**3*(3.0*a**4*n**2 - 12.0*a**4*n - 84.0*a**4 + 24.0*a**2*b**2*n**2 - 180.0*a**2*b**2*n - 456.0*a**2*b**2 - 16.0*b**4*n**2 - 48.0*b**4*n - 32.0*b**4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*b*(12.0*a**2*n**2 - 9.0*a**2*n - 261.0*a**2 + 32.0*b**2*n**2 - 80.0*b**2*n - 288.0*b**2)/(64.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0*a**3*(a**4*n**3 + 9.0*a**4*n**2 - 24.0*a**4*n - 156.0*a**4 + 16.0*a**2*b**2*n**3 + 88.0*a**2*b**2*n**2 - 264.0*a**2*b**2*n - 1152.0*a**2*b**2 + 16.0*b**4*n**3 + 32.0*b**4*n**2 - 144.0*b**4*n - 288.0*b**4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate main diagonal
    def d0(n):
        return -3.0*a**4*b*(7.0*a**2*n**2 - 93.0*a**2 + 14.0*b**2*n**2 - 126.0*b**2)/(16.0*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*a**3*(a**4*n**3 - 9.0*a**4*n**2 - 24.0*a**4*n + 156.0*a**4 + 16.0*a**2*b**2*n**3 - 88.0*a**2*b**2*n**2 - 264.0*a**2*b**2*n + 1152.0*a**2*b**2 + 16.0*b**4*n**3 - 32.0*b**4*n**2 - 144.0*b**4*n + 288.0*b**4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**4*b*(12.0*a**2*n**2 + 9.0*a**2*n - 261.0*a**2 + 32.0*b**2*n**2 + 80.0*b**2*n - 288.0*b**2)/(64.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*(3.0*a**4*n**2 + 12.0*a**4*n - 84.0*a**4 + 24.0*a**2*b**2*n**2 + 180.0*a**2*b**2*n - 456.0*a**2*b**2 - 16.0*b**4*n**2 + 48.0*b**4*n - 32.0*b**4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*b*(21.0*a**2*n + 81.0*a**2 - 8.0*b**2*n**2 - 22.0*b**2*n + 30.0*b**2)/(32.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*(a**2*n**2 - 7.0*a**2*n - 54.0*a**2 + 24.0*b**2*n**2 + 84.0*b**2*n - 108.0*b**2)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -a**6*b*(4.0*n + 21.0)/(64.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 7th superdiagonal
    def d7(n):
        return -a**7*(n + 6.0)/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def i4y4d1_diags(ny, lower, upper):
    """Create operator for I^4 y^4 D f(y) dy"""

    ns = np.arange(0, ny)
    offsets = np.arange(-7,8)
    nzrow = 3

    a, b = linear_map(lower, upper)

    # Generate 7th subdiagonal
    def d_7(n):
        return a**7*(n - 7.0)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 6th subdiagonal
    def d_6(n):
        return a**6*b*(n - 6.0)/(16.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0))

    # Generate 5th subdiagonal
    def d_5(n):
        return a**5*(n - 5.0)*(a**2*n + 13.0*a**2 + 24.0*b**2*n + 24.0*b**2)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 4th subdiagonal
    def d_4(n):
        return a**4*b*(n - 4.0)*(3.0*a**2 + b**2*n + b**2)/(4.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0))

    # Generate 3rd subdiagonal
    def d_3(n):
        return -a**3*(3.0*a**4*n**2 - 15.0*a**4*n - 102.0*a**4 + 24.0*a**2*b**2*n**2 - 216.0*a**2*b**2*n - 528.0*a**2*b**2 - 16.0*b**4*n**2 - 48.0*b**4*n - 32.0*b**4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 2nd subdiagonal
    def d_2(n):
        return -a**4*b*(3.0*a**2*n**2 - 3.0*a**2*n - 78.0*a**2 + 8.0*b**2*n**2 - 24.0*b**2*n - 80.0*b**2)/(16.0*n*(n - 3.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 1st subdiagonal
    def d_1(n):
        return -3.0*a**3*(a**4*n**3 + 10.0*a**4*n**2 - 29.0*a**4*n - 190.0*a**4 + 16.0*a**2*b**2*n**3 + 96.0*a**2*b**2*n**2 - 304.0*a**2*b**2*n - 1344.0*a**2*b**2 + 16.0*b**4*n**3 + 32.0*b**4*n**2 - 144.0*b**4*n - 288.0*b**4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate main diagonal
    def d0(n):
        return -3.0*a**4*b*(a**2*n**2 - 14.0*a**2 + 2.0*b**2*n**2 - 18.0*b**2)/(2*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 1st superdiagonal
    def d1(n):
        return 3.0*a**3*(a**4*n**3 - 10.0*a**4*n**2 - 29.0*a**4*n + 190.0*a**4 + 16.0*a**2*b**2*n**3 - 96.0*a**2*b**2*n**2 - 304.0*a**2*b**2*n + 1344.0*a**2*b**2 + 16.0*b**4*n**3 - 32.0*b**4*n**2 - 144.0*b**4*n + 288.0*b**4)/(128.0*n*(n - 3.0)*(n - 2.0)*(n - 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 2nd superdiagonal
    def d2(n):
        return a**4*b*(3.0*a**2*n**2 + 3.0*a**2*n - 78.0*a**2 + 8.0*b**2*n**2 + 24.0*b**2*n - 80.0*b**2)/(16.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 3.0))

    # Generate 3rd superdiagonal
    def d3(n):
        return a**3*(3.0*a**4*n**2 + 15.0*a**4*n - 102.0*a**4 + 24.0*a**2*b**2*n**2 + 216.0*a**2*b**2*n - 528.0*a**2*b**2 - 16.0*b**4*n**2 + 48.0*b**4*n - 32.0*b**4)/(128.0*n*(n - 2.0)*(n - 1.0)*(n + 1.0)*(n + 2.0))

    # Generate 4th superdiagonal
    def d4(n):
        return a**4*b*(n + 4.0)*(3.0*a**2 - b**2*n + b**2)/(4.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 5th superdiagonal
    def d5(n):
        return -a**5*(n + 5.0)*(a**2*n - 13.0*a**2 + 24.0*b**2*n - 24.0*b**2)/(128.0*n*(n - 1.0)*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 6th superdiagonal
    def d6(n):
        return -a**6*b*(n + 6.0)/(16.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    # Generate 7th superdiagonal
    def d7(n):
        return -a**7*(n + 7.0)/(128.0*n*(n + 1.0)*(n + 2.0)*(n + 3.0))

    ds = [d_7, d_6, d_5, d_4, d_3, d_2, d_1, d0, d1, d2, d3, d4, d5, d6, d7]
    diags, offsets = utils.build_diagonals(ns, nzrow, ds, offsets)
    return (diags, offsets)

def qid(ny, lower, upper, q, coeff = 1.0):
    """Create a quasi identity block of order q"""

    mat = spsp.coo_matrix((ny,ny))
    if coeff != 1.0:
        mat.data = coeff*np.ones((ny-q))
    else:
        mat.data = np.ones((ny-q))
    mat.row = np.arange(q,ny)
    mat.col = mat.row
    return mat

def sid(ny, lower, upper, s, coeff = 1.0):
    """Create a identity block with last s rows zeroed"""

    mat = spsp.coo_matrix((ny,ny))
    if coeff != 1.0:
        mat.data = coeff*np.ones((ny-s))
    else:
        mat.data = np.ones((ny-s))
    mat.row = np.arange(0,ny-s)
    mat.col = mat.row
    return mat

def integral(ny, lower, upper):
    """Compute the definite integral of the expansion"""

    a, b = linear_map(lower, upper)

    mat = spsp.lil_matrix((1,ny))
    mat[0,::2] = [4.0*a*(n/(n**2 - 1.0) - 1.0/(n - 1.0)) for n in np.arange(0,ny,2)]
    mat[0,0] = mat[0,0]/2.0

    return mat

def avg(ny, lower, upper):
    """Compute the average of the expansion"""

    mat = spsp.lil_matrix((1,ny))
    mat[0,::2] = [2.0*(n/(n**2 - 1.0) - 1.0/(n - 1.0)) for n in np.arange(0,ny,2)]
    mat[0,0] = mat[0,0]/2.0

    return mat
