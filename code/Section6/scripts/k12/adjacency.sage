# adjacency.sage
# Compute the 1-PS weight vector r from the given exponent set (Case k=12).
# The script finds r such that:
#   - r · i is constant for all i in the index set,
#   - sum(r_i) = 0,
# and outputs the fingerprint (−r sorted in descending order).

from sage.all import *

# Exponent set for Case k=12 (I(r_12)_{=0})
indexset = [
    (0, 1, 0, 0, 2, 0, 0),
    (0, 0, 2, 0, 0, 1, 0),
    (0, 0, 1, 1, 0, 1, 0),
    (0, 0, 0, 2, 0, 1, 0),
    (1, 0, 0, 0, 1, 1, 0),
    (0, 2, 0, 0, 0, 0, 1),
    (1, 0, 1, 0, 0, 0, 1),
    (1, 0, 0, 1, 0, 0, 1),
]

def primitive_integer_vector(v):
    """Convert a rational vector to a primitive integer vector."""
    L = lcm([c.denominator() for c in v])
    w = [(c * L).numerator() for c in v]
    g = gcd([abs(x) for x in w if x != 0]) or 1
    return vector(ZZ, [x // g for x in w])

# Construct matrix with difference vectors and the sum=0 constraint
E = [vector(QQ, e) for e in indexset]
i0 = E[0]
rows = [(E[j] - i0) for j in range(1, len(E))]
rows.append(vector(QQ, [1] * 7))
M = Matrix(QQ, rows)

# Compute the right kernel (should be one-dimensional)
ker = M.right_kernel()
if ker.dimension() != 1:
    raise ValueError("Right kernel dimension is not 1")
r = primitive_integer_vector(ker.basis()[0])

# Normalize orientation: enforce non-increasing order
if not all(r[i] >= r[i + 1] for i in range(6)):
    r = -r

# Output fingerprint (−r sorted in descending order)
fingerprint = sorted([-int(ri) for ri in r], reverse=True)
print([fingerprint])
