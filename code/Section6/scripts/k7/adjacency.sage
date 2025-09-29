# adjacency_k7.sage
# ------------------------------------------------------------
# Purpose
#   Given a set of exponent vectors (indexset) that lie on a face of the
#   cubic simplex and contain the barycenter, compute a normal vector r
#   (1-PS weights) orthogonal to that face. This detects the T-face used
#   in the adjacency computation for cubic fivefolds.
#
# Method
#   Let I ⊂ Z^7_(>=0) be degree-3 exponent vectors. For a face passing
#   through the barycenter η = (3/7,...,3/7), the normal vector r must
#   satisfy
#       r · e_i = 0  for every e_i in the chosen indexset,
#       sum(r_j) = 0 (equivalently, r · η = 0).
#   Hence, stacking the exponent rows together with a row of ones and
#   taking the right kernel gives a 1-dimensional space spanned by r.
#
# Output (Mathematica-style)
#   We print a single list [[...]] which is:
#     (i) the primitive integer generator of the kernel,
#     (ii) multiplied by (-1),
#     (iii) sorted in non-increasing order (descending).
#   This matches the sign/order convention the user’s Mathematica code used.
#
# Notes
#   - All arithmetic is done over QQ; conversion to a primitive ZZ-vector
#     uses an lcm of denominators followed by division by gcd of entries.
#   - The script is self-contained and has no command-line options, as requested.
# ------------------------------------------------------------

from sage.all import *

# -----------------------------
# Utilities
# -----------------------------
def _gcd_list(int_list):
    """gcd over a (possibly signed) list of Python/Sage integers."""
    from functools import reduce
    if not int_list:
        return Integer(1)
    return reduce(lambda a, b: gcd(a, b), [abs(Integer(x)) for x in int_list])

def primitive_integer_vector(v):
    """
    Scale a rational vector v ∈ QQ^n to a primitive integer vector in ZZ^n
    (divide by the gcd of entries). If v = 0, return v.
    """
    vv = vector(QQ, v)
    if vv.is_zero():
        return vector(ZZ, [0]*len(vv))
    den_lcm = lcm([c.denominator() for c in vv])
    w = vector(ZZ, [ (c*den_lcm).numerator() for c in vv ])  # integral now
    g = _gcd_list(list(w))
    if g != 0:
        w = vector(ZZ, [zi // g for zi in w])
    return w

def mathematica_style(v):
    """
    Convert vector v to the 'Mathematica-style' the user requested:
      - multiply by (-1)
      - sort entries in non-increasing (descending) order
    Return as a plain Python list.
    """
    w = [-int(a) for a in v]
    w.sort(reverse=True)
    return w

def check_indexset(indexset):
    """
    Basic sanity checks:
      - each exponent vector has length 7,
      - degree is 3 (sum of entries = 3),
      - entries are integers ≥ 0.
    """
    for row in indexset:
        if len(row) != 7:
            raise ValueError("Each exponent vector must have length 7.")
        if sum(row) != 3:
            raise ValueError("Each exponent vector must have total degree 3.")
        if any((x < 0) for x in row):
            raise ValueError("Exponent vectors must be nonnegative.")

# -----------------------------
# Core computation
# -----------------------------
def normal_vector_from_indexset(indexset):
    """
    Build the (m+1) × 7 matrix whose rows are the m exponent vectors
    in 'indexset' together with the row (1,1,1,1,1,1,1). The right kernel
    is expected to be 1-dimensional. Return its primitive integer generator.
    """
    check_indexset(indexset)
    M = matrix(QQ, indexset)
    ones = matrix(QQ, [[1,1,1,1,1,1,1]])  # sum(r_i) = 0  <=>  (1,...,1)·r = 0
    A = M.stack(ones)
    K = A.right_kernel()
    if K.dimension() != 1:
        raise RuntimeError("Right kernel is not 1-dimensional. "
                           "Please re-check the indexset; got dim = {}".format(K.dimension()))
    r = primitive_integer_vector(K.basis()[0])
    return r

# -----------------------------
# k = 7 : indexset (given by the user)
# -----------------------------
indexset = [
    [0, 0, 0, 2, 1, 0, 0],
    [0, 1, 0, 0, 2, 0, 0],
    [0, 0, 1, 1, 0, 1, 0],
    [1, 0, 0, 0, 0, 2, 0],
    [0, 2, 0, 0, 0, 0, 1],
    [1, 0, 1, 0, 0, 0, 1],
]

# -----------------------------
# Run and print (Mathematica-style)
# -----------------------------
r = normal_vector_from_indexset(indexset)
out = mathematica_style(r)
print([out])
