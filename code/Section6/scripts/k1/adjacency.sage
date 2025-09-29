# ================================================================
# Cubic fivefold (in P^6): 1-step specializations (loop over k)
# ---------------------------------------------------------------
# This SageMath script loops over k=1..21. For each k it builds
#   Indexset := I(r_k)_{=0}
# and computes candidate facet normals ("walls") in the 5D slice
#   A_k = { e ∈ R^7 : sum(e_i)=3,  r_k·e=0 }.
#
# Key points:
#   • We solve for an *affine* supporting hyperplane w·e = c with
#     constraints  w·1=0, w·r_k=0, and w·e_j=c for 5 chosen points e_j.
#   • Orientation is fixed so that w·e - c ≥ 0 for all e in Indexset.
#   • A facet is accepted iff its equality set has ≥5 points with
#     affine dimension ≥4 inside A_k.
#
# Usage:
#   sage adjacency_step1_loop.sage
#
# Public release: comments are in English; license can be set to MIT.
# ================================================================

from itertools import combinations

# ---------------------------
# Global configuration
# ---------------------------
DIM     = 7     # number of variables x0,...,x6
DEGREE  = 3     # cubic
K_RANGE = range(1, 22)  # k = 1..21 ; r22 ≃ r21 (unused by default)

# ---------------------------
# All exponent vectors of total degree = 3
# (weak compositions: nonnegative entries summing to DEGREE)
# ---------------------------
ALL_EXPONENTS = [list(v) for v in IntegerVectors(DEGREE, DIM)]  # 84 vectors

# ---------------------------
# 1-PS weight list r_k  (as in Prop. 2.3)
# ---------------------------
R_LIST = {
  1:  ( 8, 3, 2, -1, -2, -4, -6),
  2:  ( 6, 4, 1, -1, -2, -3, -5),
  3:  ( 4, 2, 1, -1, -1, -2, -3),
  4:  ( 3, 2, 1,  0, -1, -2, -3),
  5:  ( 4, 2, 1,  0, -1, -2, -4),
  6:  ( 5, 3, 2,  1, -1, -4, -6),
  7:  ( 6, 4, 2,  1, -2, -3, -8),
  8:  ( 4, 1, 1,  0, -2, -2, -2),
  9:  ( 2, 2, 0,  0, -1, -1, -2),
 10:  ( 2, 1, 0,  0, -1, -1, -1),
 11:  ( 2, 0, 0,  0,  0, -1, -1),
 12:  ( 3, 2, 1,  1, -1, -2, -4),
 13:  ( 2, 1, 1,  0, -1, -1, -2),
 14:  ( 2, 2, 0, -1, -1, -1, -1),
 15:  ( 2, 1, 1,  0,  0, -2, -2),
 16:  ( 2, 1, 0,  0,  0, -1, -2),
 17:  ( 1, 1, 1,  0,  0, -1, -2),
 18:  ( 1, 1, 0,  0,  0, -1, -1),
 19:  ( 2, 2, 2,  0, -1, -1, -4),
 20:  ( 1, 1, 1,  1,  0, -2, -2),
 21:  ( 1, 1, 0,  0,  0,  0, -2),
 22:  ( 1, 0, 0,  0,  0,  0, -1)  # ≃ r21; normally unused
}

def R(k):
    """Return r_k as a QQ^DIM vector."""
    return vector(QQ, R_LIST[k])

# ---------------------------
# Build Indexset = I(r)_{=0}
# ---------------------------
def I_eq0_exponents(r):
    """
    Return I(r)_{=0}: all exponent vectors e with sum(e)=DEGREE and r·e=0.
    Output: list of Python lists (each length DIM).
    """
    rv = vector(QQ, r)
    out = []
    for e in ALL_EXPONENTS:
        if rv * vector(QQ, e) == 0:
            out.append(e)
    return out

# ---------------------------
# Utilities: normalization, dot, affine dimension in the slice
# ---------------------------
def primitive_integer_vector(v):
    """
    Normalize a rational vector to a primitive integer vector
    with gcd=1 and first nonzero entry positive. Returns a tuple.
    """
    vv = vector(QQ, v)
    L  = Integer(1)
    for x in vv:
        L = lcm(L, Integer(x.denominator()))
    ints = [Integer(L * x) for x in vv]
    g = Integer(0)
    for x in ints:
        g = gcd(g, abs(x))
    if g != 0:
        ints = [x // g for x in ints]
    # fix orientation
    for x in ints:
        if x != 0:
            if x < 0:
                ints = [-y for y in ints]
            break
    return tuple(int(x) for x in ints)

def dot(v, e):
    """Dot product (QQ): v·e with v a QQ^DIM vector/iterable, e a list/tuple."""
    return vector(QQ, v) * vector(QQ, e)

def affine_dim_in_slice(points):
    """
    Affine dimension inside A_k (sum=3 and r_k·e=0) computed as rank of
    difference vectors (points[i]-points[0]) over QQ.
    (Differences automatically lie in the slice's tangent space.)
    """
    if len(points) <= 1:
        return 0
    base = vector(QQ, points[0])
    diffs = [vector(QQ, p) - base for p in points[1:]]
    if not diffs:
        return 0
    M = matrix(QQ, diffs)
    return M.rank()

# ---------------------------
# Core: compute facet normals from a given Indexset and r_k
# ---------------------------
def facet_normals_from_indexset(Indexset, rk, facet_min_equal=5):
    """
    Given Indexset = I(r_k)_{=0} and rk, compute facet normals in the 5D slice A_k
    by solving for (w,c) with constraints:
        w·1 = 0,  w·rk = 0,  and  w·e_j = c  for five chosen points e_j.
    Keep (w,c) if it supports all of Indexset (nonnegativity after orientation)
    and if the equality set has ≥ facet_min_equal points with affine dimension ≥4.
    Returns a sorted list of primitive integer directions (tuples).
    """
    ones = vector(QQ, [1]*DIM)
    if len(Indexset) < 5:
        return []

    normals = set()

    for combo in combinations(Indexset, 5):
        # Build the 7×8 system A * (w,c)^T = 0 in right-kernel form:
        # rows: [ones | 0], [rk | 0], and [e_j | -1] for j=1..5
        rows = []
        rows.append(list(ones) + [QQ(0)])
        rows.append(list(vector(QQ, rk)) + [QQ(0)])
        for e in combo:
            rows.append(list(vector(QQ, e)) + [QQ(-1)])

        A = matrix(QQ, rows)
        K = A.right_kernel()
        if K.dimension() != 1:
            continue

        u = K.basis()[0]                  # u = (w, c)
        w = vector(QQ, list(u)[:-1])
        c = QQ(u[-1])

        # Orient so that w·e - c ≥ 0 for all e in Indexset
        vals = [dot(w, e) - c for e in Indexset]
        if all(v >= 0 for v in vals):
            w_oriented, c_oriented = w, c
        elif all(v <= 0 for v in vals):
            w_oriented, c_oriented = -w, -c
        else:
            continue  # not a supporting functional

        # Equality set on Indexset
        zero_pts = [e for e in Indexset if dot(w_oriented, e) - c_oriented == 0]

        # Facet test (5D slice): need ≥5 points and affine dimension ≥4
        if len(zero_pts) >= facet_min_equal and affine_dim_in_slice(zero_pts) >= 4:
            normals.add(primitive_integer_vector(w_oriented))

    return sorted(normals)

# ---------------------------
# Main loop over k = 1..21
# ---------------------------
if __name__ == "__main__":
    normals_per_k = {}
    any_found = False

    for k in K_RANGE:
        rk = R(k)
        Indexset = I_eq0_exponents(rk)      # <-- THIS is your Indexset for r_k

        # >>> Hook <<< If your existing script expects a global `Indexset`
        # and then runs some downstream logic, you can insert it here, e.g.:
        #   result_k = your_existing_wall_logic(Indexset)
        # For a self-contained run, we compute walls via (w,c) as below.
        normals_k = facet_normals_from_indexset(Indexset, rk, facet_min_equal=5)
        normals_per_k[k] = normals_k
        if normals_k:
            any_found = True

    if not any_found:
        print("1-step specialization: none")
    else:
        print("1-step specialization: found")
        nonempty = {k: v for k, v in normals_per_k.items() if v}
        print(f"Number of k with detected walls: {len(nonempty)}")
        for k in sorted(nonempty):
            print(f"k={k} -> {nonempty[k]}")
