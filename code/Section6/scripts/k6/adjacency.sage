# -*- coding: utf-8 -*-
# SageMath script
#
# Purpose
# -------
# Compute the primitive 1-PS normal vector r for the face spanned by a given
# index set of cubic monomials (here: the k = 6 case for cubic fivefolds).
# The index set consists of exponent vectors i = (i0,...,i6) with sum(i) = 3.
#
# Method
# ------
# If S = {i^(1), ..., i^(m)} is the index set, any normal vector r to the face
# satisfies r·(i^(a) - i^(1)) = 0 for a = 2,...,m. In addition, for SL(7) we
# impose sum_j r_j = 0. Solving these linear equations over QQ yields a
# 1-dimensional kernel; its primitive integer generator is the desired r.
#
# Output convention
# -----------------
# To match the Mathematica output used previously, we *only* print the final
# vector after multiplying by (-1) and sorting in descending order. (This loses
# coordinate positions intentionally, per the user’s request.)
#
# Notes
# -----
# - All comments are in English for public release.
# - No command-line options: running `sage adjacency.sage` prints the result.

from sage.all import matrix, vector, QQ, ZZ, lcm, gcd

def primitive_integer_vector(vv):
    """
    Convert a QQ-vector vv to a primitive ZZ-vector.
    Steps: clear denominators, then divide by gcd of entries; fix the sign so
    the first nonzero entry is positive.
    """
    den_lcm = 1
    for x in vv:
        # x is a rational in QQ; denominator() is positive
        den_lcm = lcm(den_lcm, x.denominator())
    w = [ZZ(x * den_lcm) for x in vv]
    g = 0
    for z in w:
        g = gcd(g, abs(z))
    if g != 0:
        w = [z // g for z in w]
    # Normalize sign: first nonzero should be positive
    for z in w:
        if z != 0:
            if z < 0:
                w = [-t for t in w]
            break
    return vector(ZZ, w)

def normal_vector_from_indexset(indexset):
    """
    Given an index set (list of 7-tuples summing to 3), build the linear system:
        r · (i^(a) - i^(1)) = 0 for a = 2..m
        sum_j r_j = 0
    Solve for r in the right kernel over QQ and return a primitive ZZ-vector.
    """
    if len(indexset) < 2:
        raise ValueError("indexset must contain at least two exponent vectors.")
    base = vector(QQ, indexset[0])
    rows = []
    for idx in indexset[1:]:
        rows.append(vector(QQ, idx) - base)
    # SL(7) trace-zero constraint: sum r_j = 0
    rows.append(vector(QQ, [1]*7))
    M = matrix(QQ, rows)
    ker = M.right_kernel()
    basis = ker.basis()
    if len(basis) != 1:
        raise RuntimeError("Kernel dimension is not 1; got %d" % len(basis))
    r = primitive_integer_vector(basis[0])
    return r

def mathematica_style(v):
    """
    Return the list obtained by multiplying v by -1 and sorting in descending order.
    This matches the requested display style (ignores original coordinate order).
    """
    lst = [-int(x) for x in v]
    lst.sort(reverse=True)
    return lst

# ----------------------------------------------------------------------
# k = 6 : index set provided by the user
# ----------------------------------------------------------------------
indexset = [
    [0, 0, 1, 0, 2, 0, 0],
    [0, 0, 2, 0, 0, 1, 0],
    [0, 1, 0, 1, 0, 1, 0],
    [1, 0, 0, 0, 1, 1, 0],
    [0, 2, 0, 0, 0, 0, 1],
    [1, 0, 0, 1, 0, 0, 1],
]

if __name__ == "__main__":
    r = normal_vector_from_indexset(indexset)
    # Print only the Mathematica-style vector (as a singleton list for consistency)
    print([mathematica_style(r)])
