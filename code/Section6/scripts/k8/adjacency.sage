# adjacency_k8.sage
# -------------------------------------------
# Purpose:
#   Recover the primitive 1-PS weight vector r \in Z^7 for the
#   given exponent index set (a face through the barycenter),
#   and print it in the "Mathematica-style" convention:
#     multiply by -1 and sort in descending order.
#
# Method:
#   If the face F is given by exponent vectors v_i in Z_{\ge 0}^7 with sum(v_i)=3,
#   then its affine span is orthogonal to a unique (up to scale) normal r.
#   Adding the all-ones row enforces sum(r)=0 (barycenter condition).
#   We compute the right kernel over QQ and extract a primitive integer vector.
#
# Output:
#   A single line like:
#       [[4, 1, 1, 0, -2, -2, -2]]
#   (No adjacency pairs are printed.)
#
# Reference:
#   Y. Shibata, "The Boundary of the Moduli Space of Stable Cubic Fivefolds",
#   esp. §2 (Lemma 2.1, Algorithm 2.2) and §3 for the list of r_k.  See the
#   accompanying preprint for the k=8 case.  (This script is a minimal,
#   self-contained extractor of the 1-PS vector from a face.)
# -------------------------------------------

# ---- Utility: primitive integer vector from a QQ-vector
def primitive_integer_vector(vv):
    """
    Take vv \in QQ^7 (kernel generator) and return a primitive vector in ZZ^7.
    """
    qv = [QQ(x) for x in vv]
    den_lcm = lcm([x.denominator() for x in qv]) if qv else 1
    ints = [int(x * den_lcm) for x in qv]
    # Make primitive by dividing by gcd of entries (ignoring zeros)
    nz = [abs(z) for z in ints if z != 0]
    g = gcd(nz) if nz else 1
    ints = [z // g for z in ints]
    # Canonicalize sign minimally (first nonzero > 0)
    for z in ints:
        if z != 0:
            if z < 0:
                ints = [-w for w in ints]
            break
    return vector(ZZ, ints)

# ---- Core: compute a normal r from the affine-difference matrix + ones
def normal_from_indexset(indexset):
    """
    Given indexset (list of 7-tuples in ZZ), build a matrix whose rows are
    all affine differences (v_i - v_0) together with the all-ones row.
    Return a primitive integer generator of its right kernel.
    """
    # unique points (defensive)
    pts = [tuple(p) for p in indexset]
    pts = sorted(set(pts))
    if len(pts) < 2:
        raise ValueError("indexset must contain at least two distinct points.")

    v0 = vector(QQ, pts[0])
    diff_rows = [vector(QQ, tuple(p)) - v0 for p in pts[1:]]
    ones = vector(QQ, [1]*7)  # enforces sum(r) = 0 (barycenter condition)

    M = matrix(QQ, diff_rows + [ones])   # rows x 7
    K = M.right_kernel()
    basis = K.basis()
    if len(basis) == 0:
        raise RuntimeError("Right kernel is trivial; check indexset.")
    # If kernel has dimension > 1, we take the first generator; for faces used
    # in this context, it is generically 1-dimensional.
    r = primitive_integer_vector(basis[0])
    return r

# ---- Pretty-print in the requested convention (flip sign and sort desc)
def mathematica_style(r):
    """
    Multiply by -1 and sort in descending order (as requested).
    Return a Python list of ints.
    """
    arr = [-int(z) for z in r]       # flip sign
    arr_sorted = sorted(arr, reverse=True)  # descending
    return arr_sorted

# ------------------------------
# k = 8 : set the given indexset
# ------------------------------
indexset = [
    (0, 0, 0, 3, 0, 0, 0),
    (0, 2, 0, 0, 1, 0, 0),
    (0, 1, 1, 0, 1, 0, 0),
    (0, 0, 2, 0, 1, 0, 0),
    (1, 0, 0, 0, 2, 0, 0),
    (0, 2, 0, 0, 0, 1, 0),
    (0, 1, 1, 0, 0, 1, 0),
    (0, 0, 2, 0, 0, 1, 0),
    (1, 0, 0, 0, 1, 1, 0),
    (1, 0, 0, 0, 0, 2, 0),
    (0, 2, 0, 0, 0, 0, 1),
    (0, 1, 1, 0, 0, 0, 1),
    (0, 0, 2, 0, 0, 0, 1),
    (1, 0, 0, 0, 1, 0, 1),
    (1, 0, 0, 0, 0, 1, 1),
    (1, 0, 0, 0, 0, 0, 2)
]

# ---- Compute r and print only the Mathematica-style list (no adjacency pairs)
r = normal_from_indexset(indexset)
out = mathematica_style(r)

# Print in the same shape as earlier scripts: a list containing one vector
print([out])
