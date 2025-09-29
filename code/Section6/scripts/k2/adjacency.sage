# -*- coding: utf-8 -*-
# adjacency.sage
#
# Title:
#   Adjacency detector for k = 2 (cubic fivefolds, Section 6)
#
# Purpose:
#   Given six degree-3 exponent vectors (the monomials used for Case k = 2),
#   form every 5-subset together with the barycenter row η = (1,...,1) to make
#   a 6×7 matrix. The right kernel is one-dimensional; its primitive integer
#   generator r ∈ Z^7 is (up to sign) the normal of the corresponding face
#   and coincides with a 1‑PS weight. We then print a positionless “fingerprint”
#   of r in the paper’s convention:
#      (i) clear denominators → primitive integer vector,
#     (ii) fix sign so the first nonzero entry is positive,
#    (iii) for printing only: multiply by −1 and sort entries in decreasing order.
#
#   For k = 2, this procedure returns [[5, 3, 2, 1, -1, -4, -6]], which matches
#   the neighbor weight r6 in the adjacency pair {Φ2, Φ6} (Section 6).  See also
#   the weight identity (2) listing the six monomials of ϕ2.  [SageMath script
#   independent of external packages.]
#
# How to run:
#   sage adjacency.sage
#
# Output:
#   One or more fingerprints (deduplicated). For k = 2 with the exponent set
#   below, the expected output is:
#       [[5, 3, 2, 1, -1, -4, -6]]
#
# References:
#   Y. Shibata, *The Boundary of the Moduli Space of Stable Cubic Fivefolds*,
#   esp. Section 6 (adjacency/wall crossing) and the k = 2 data. :contentReference[oaicite:1]{index=1}

from itertools import combinations

# ===== Parameters (ambient) =====
dim = 7   # variables x0,...,x6
deg = 3   # total degree of monomials (cubic); kept for clarity

# ===== Exponent set for k = 2 (six monomials of φ₂) =====
#   x2^2*x4, x1*x4^2, x1*x3*x5, x0*x5^2, x1*x2*x6, x0*x3*x6
#   Each vector lists exponents of (x0, x1, x2, x3, x4, x5, x6).
indexset = [
    vector(ZZ, [0, 0, 2, 0, 1, 0, 0]),  # x2^2 * x4
    vector(ZZ, [0, 1, 0, 0, 2, 0, 0]),  # x1 * x4^2
    vector(ZZ, [0, 1, 0, 1, 0, 1, 0]),  # x1 * x3 * x5
    vector(ZZ, [1, 0, 0, 0, 0, 2, 0]),  # x0 * x5^2
    vector(ZZ, [0, 1, 1, 0, 0, 0, 1]),  # x1 * x2 * x6
    vector(ZZ, [1, 0, 0, 1, 0, 0, 1]),  # x0 * x3 * x6
]

# ===== Utilities =====
def primitive_integer_vector(v):
    """
    Convert a rational vector v ∈ QQ^n to a primitive integer vector in ZZ^n:
      (1) clear denominators (multiply by lcm of denominators and take numerators);
      (2) divide by gcd of nonzero entries;
      (3) fix sign so the first nonzero entry is positive.
    Returns: vector(ZZ, n)
    """
    vv = vector(QQ, v)
    if vv.is_zero():
        return vector(ZZ, [0]*vv.length())

    den_lcm = lcm([x.denominator() for x in vv])
    w = vector(ZZ, [(x * den_lcm).numerator() for x in vv])

    g = gcd([abs(a) for a in w if a != 0])
    if g != 0:
        w = vector(ZZ, [a // g for a in w])

    for a in w:
        if a != 0:
            if a < 0:
                w = -w
            break
    return w

def fingerprint_from_nullspace(M):
    """
    Take a 6×7 matrix M (top row = all-ones), compute a basis of the 1-dim
    right kernel over QQ, normalize to a primitive integer vector, and return a
    positionless “fingerprint”: the entries sorted in decreasing order.
    """
    K = Matrix(QQ, M).right_kernel()
    if K.dimension() != 1:
        return None
    v = primitive_integer_vector(K.basis()[0])
    return sorted(list(v), reverse=True)

def to_paper_convention(desc_sorted_list):
    """
    Apply the paper/Mathematica printing convention:
      multiply by −1 and sort in decreasing order.
    """
    flipped = [-x for x in desc_sorted_list]
    return sorted(flipped, reverse=True)

# ===== Main routine =====
def main(print_debug=False):
    """
    Build 6×7 matrices with top row (1,...,1) and the remaining five rows a
    5-subset of the six exponent vectors. For each rank-6 matrix, extract the
    1‑PS direction from the right kernel and print it in the paper’s convention.
    """
    all_ones = vector(ZZ, [1]*dim)
    result_set = set()

    for subset in combinations(indexset, dim - 2):  # choose 5 of the 6
        rows = [all_ones] + list(subset)
        M = Matrix(ZZ, rows)
        if print_debug:
            rk = M.rank()
            kd = Matrix(QQ, M).right_kernel().dimension()
            print("rank(M) =", rk, "kernel_dim =", kd)
        if M.rank() == dim - 1:  # rank 6 ⇒ 1-dimensional kernel
            fp = fingerprint_from_nullspace(M)  # descending, signless
            if fp is not None:
                paper_vec = to_paper_convention(fp)
                result_set.add(tuple(paper_vec))

    # Deduplicate and print in a stable order
    formatted = [list(v) for v in sorted(result_set, reverse=True)]
    print(formatted)

if __name__ == "__main__":
    # Set print_debug=True to see ranks and kernel dimensions for each subset
    main(print_debug=False)
