#!/usr/bin/env sage
# -*- coding: utf-8 -*-
# Combined Sage script (runs normalVectors first, then maximalList)
# - First half (normalVectors): reuse the same counter (progress_interval) as in the uploaded program
# - Second half (maximalList): read normalVectors.m and compute maximalList
# Generated: 2025-09-25T07:55:16

# ======= BEGIN: normalVectors.sage (with original counter kept) =======
# git_projective_hypersurfaces.sage
# -*- coding: utf-8 -*-

from itertools import combinations

# ======= Parameters (change as needed) =======
dim = 7
deg = 3
# ===========================================

def generate_index_set(dim, deg):
    """
    Mathematica: equivalent to selecting only those from Tuples[Range[0, deg], dim] whose sum equals deg.
    In Sage, IntegerVectors(deg, dim) is efficient.
    """
    return [list(v) for v in IntegerVectors(deg, dim)]

def lcm2(a, b):
    from math import gcd as _g
    return a // _g(a, b) * b

def lcm_list(L):
    L = [int(abs(x)) for x in L if int(abs(x)) != 0]
    if not L:
        return 1
    from functools import reduce
    return reduce(lcm2, L, 1)

def gcd_list(L):
    from math import gcd as _g
    g = 0
    for x in L:
        g = _g(g, int(abs(x)))
    return g

def primitive_sorted(v):
    """
    For a basis vector v of the null space (in QQ^dim),
      1) multiply by the LCM of denominators to clear fractions,
      2) divide by the gcd to make it primitive,
      3) sort its components in descending order (matches Mathematica's Sort[..., Greater]),
    and return the resulting tuple.
    """
    den_lcm = lcm_list([c.denominator() for c in v])
    ints = [ZZ(c * den_lcm) for c in v]
    g = gcd_list(ints)
    if g > 0:
        ints = [int(x // g) for x in ints]
    ints.sort(reverse=True)     # ← this matches the Mathematica convention
    return tuple(ints)

def compute_normal_vectors(dim, deg, out_path="normalVectors-sage.m",
                           progress_interval=5_000_000):
    indexset = generate_index_set(dim, deg)
    ones = [1]*dim

    seen = set()
    count = 0

    for block in combinations(indexset, dim-2):   # Mathematica: equivalent to NextSubset
        rows = [ones] + [list(u) for u in block]  # put {1,1,...,1} as the top row
        A = matrix(QQ, rows)

        if A.rank() == dim - 1:
            ker = A.right_kernel()
            if ker.dimension() == 1:              # use only when the nullspace has dimension 1
                v = ker.basis()[0]
                seen.add(primitive_sorted(v))

        count += 1
        if count % progress_interval == 0:
            # same idea as Mathematica's Print[Length[result]]
            print(len(seen))

    # Mathematica's Union performs deduplication and canonical (lexicographic ascending) ordering
    result = sorted(seen)

    # Save in Mathematica format (with braces)
    with open(out_path, "w") as f:
        f.write('{' + ', '.join('{' + ', '.join(str(x) for x in row) + '}' for row in result) + '}')

    print("Length of normalVectors is", len(result))
    return result

# Compute here when the script is executed directly
# ======= END: normalVectors.sage =======

# ======= BEGIN: maximalList.sage =======
# -*- coding: utf-8 -*-
# SageMath script to compute maximalList from normalVectors and indexset.
# Usage:
#   sage maximalList.sage
#
# Output:
#   maximalList.m  (Mathematica format)

# (When running as a .sage file, 'from sage.all_cmdline import *' is automatically inserted.)
# Python standard library
import os, ast, re
from typing import Any, Iterable, List, Tuple

# ------------------------------ (1) .m I/O ------------------------------

def load_mathematica_list(path: str) -> Any:
    """
    Read an .m file that contains only a Mathematica { ... } list and convert it into a Python list.
    Remove comments (* ... *), replace { and } with [ and ], then feed it to ast.literal_eval.
    """
    with open(path, 'r', encoding='utf-8') as f:
        s = f.read()
    # remove comments of the form (* ... *) (non-greedy, possibly multiline)
    s = re.sub(r'\(\*.*?\*\)', '', s, flags=re.DOTALL)
    # turn Mathematica braces into a Python list
    s = s.replace('{', '[').replace('}', ']')
    # tolerate a trailing semicolon
    s = s.strip()
    if s.endswith(';'):
        s = s[:-1]
    return ast.literal_eval(s)

def to_mathematica(obj: Any) -> str:
    """
    Convert a Python object (nested list/tuple, int) into a Mathematica-style { ... } string.
    """
    # Simple handling to treat Sage integers as ints
    try:
        from sage.all import Integer as SageInteger, Rational
    except Exception:
        class SageInteger: pass
        class Rational: pass

    if isinstance(obj, (list, tuple)):
        return '{' + ', '.join(to_mathematica(x) for x in obj) + '}'
    # Fractions should not appear here, but handle just in case
    if 'Rational' in globals() and isinstance(obj, Rational):
        return f"{obj.numerator()}/{obj.denominator()}"
    if isinstance(obj, (int,)) or obj.__class__.__name__ == 'Integer':
        return str(int(obj))
    return str(obj)

# --------------- (2) Reproducing Mathematica's Union[{x}, list] ----------------

def _lex_key(e: Any) -> Tuple:
    """
    Mathematica's Union sorts by the default (exact) lexicographic order.
    Here, if an element is a list/tuple we convert it to a tuple; otherwise we make it a 1-element tuple,
    so that we can stably sort using Python's tuple comparison.
    """
    if isinstance(e, (list, tuple)):
        return tuple(e)
    return (e,)

def _canon_set_of_vectors(S: Iterable[Iterable[int]]) -> Tuple[Tuple[int, ...], ...]:
    """
    Assuming x is a list of vectors (tuples), convert each element to a tuple
    and return a lexicographically sorted "normalized form" to use for set equality and as a sort key.
    """
    return tuple(sorted(tuple(v) for v in S))

def union_elem_and_list_like_mathematica(elem: List[List[int]],
                                         lst: Iterable[List[List[int]]]
                                         ) -> List[List[List[int]]]:
    """
    Equivalent to Mathematica's Union[{elem}, lst].
    - form {elem} ∪ lst
    - remove exact duplicates
    - sort lexicographically as a set
    """
    all_items = [elem] + list(lst)
    # Build a map from normalized form to the original object to deduplicate
    uniq = {}
    for it in all_items:
        uniq[_canon_set_of_vectors(it)] = it
    # Sort before returning (Mathematica's Union also sorts)
    return [uniq[k] for k in sorted(uniq.keys())]

# ---------------- (3) maximalFunction (extract only inclusion-maximal sets) ---------------

def is_subset_set_of_vectors(A: List[List[int]], B: List[List[int]]) -> bool:
    """
    A and B are both lists of vectors. Check whether A ⊆ B (element equality is by vectors).
    Equivalent to Mathematica's Intersection[A, B] == A.
    """
    setA = set(map(tuple, A))
    setB = set(map(tuple, B))
    return setA.issubset(setB)

def maximalFunction(m: List[List[List[int]]]) -> List[List[List[int]]]:
    """
    Mathematica version:
      For i != j with m[[i]] != m[[j]] and Intersection[m[[i]], m[[j]]] == m[[i]],
      exclude such cases and return only those that are inclusion-maximal.
    """
    out: List[List[List[int]]] = []
    for i in range(len(m)):
        cand = m[i]
        maximal = True
        for j in range(len(m)):
            if i == j:
                continue
            if m[i] != m[j] and is_subset_set_of_vectors(m[i], m[j]):
                maximal = False
                break
        if maximal:
            # Avoid duplicates to match the effect of Union
            if not any(_canon_set_of_vectors(cand) == _canon_set_of_vectors(z) for z in out):
                out.append(cand)
    # Sort lexicographically as in Union (optional but for stability)
    out = [z for _, z in sorted(((_canon_set_of_vectors(z)), z) for z in out)]
    return out

# ------------------------------ (4) positive -------------------------------

def make_positive_with_r(r: List[int]):
    """
    Implement positive(u_) := u . r >= 0 in Python; r is captured by the closure.
    """
    def positive(u: List[int]) -> bool:
        return sum(ui * ri for ui, ri in zip(u, r)) >= 0
    return positive

# ------------------------------ (5) main ----------------------------------

def weak_compositions(n: int, k: int):
    """
    Weak compositions (split n into k nonnegative integers). A generator.
    """
    if k == 1:
        yield (n,)
        return
    for i in range(n + 1):
        for rest in weak_compositions(n - i, k - 1):
            yield (i,) + rest

def main():
    # (a) Load normalVectors
    NV_PATH = 'normalVectors.m'
    if not os.path.exists(NV_PATH):
        raise FileNotFoundError(f"{NV_PATH} が見つかりません。先に normalVectors を作成してください。")
    normalVectors: List[List[int]] = load_mathematica_list(NV_PATH)

    # (b) Load or generate indexset
    IX_PATH = 'indexset.m'
    if os.path.exists(IX_PATH):
        indexset: List[List[int]] = load_mathematica_list(IX_PATH)
    else:
        # If indexset.m is missing, regenerate it (default cubic: deg=3)
        # Change DEGREE as needed.
        DEGREE = 3  # e.g., cubic; for quartic set 4
        DIM = len(normalVectors[0])
        indexset = [list(t) for t in weak_compositions(DEGREE, DIM)]
        print(f"[info] indexset.m が無いので再生成しました: dim={DIM}, deg={DEGREE}, |indexset|={len(indexset)}")

    # (c) Main loop (equivalent to Mathematica's For)
    mlist: List[List[List[int]]] = []  # list = {}
    for i in range(len(normalVectors)):
        r = normalVectors[i]
        positive = make_positive_with_r(r)
        # x = Select[indexset, positive]
        x = [u for u in indexset if positive(u)]
        # list = maximalFunction[Union[{x}, list]]
        merged = union_elem_and_list_like_mathematica(x, mlist)
        mlist = maximalFunction(merged)

    maximalList = mlist
    print("Length of maximalList is", len(maximalList))

    # (d) Write output
    with open('maximalList.m', 'w', encoding='utf-8') as f:
        f.write(to_mathematica(maximalList) + '\n')

# ======= END: maximalList.sage =======

# ======= RUN ORDER (entry point when run as a script) =======
if __name__ == "__main__":
    # Use the same counter behavior (progress_interval argument) as the uploaded normalVectors.sage.
    # To adjust the print frequency, change PROGRESS_INTERVAL below.
    PROGRESS_INTERVAL = 5_000_000   # e.g., set 100_000 for more frequent progress output
    
    # Match the output filename expected by maximalList
    NV_OUT = "normalVectors.m"
    
    # First: run normalVectors
    compute_normal_vectors(dim, deg, out_path=NV_OUT, progress_interval=PROGRESS_INTERVAL)
    
    # Second: run maximalList
    # (calls the main() function inside maximalList)
    main()