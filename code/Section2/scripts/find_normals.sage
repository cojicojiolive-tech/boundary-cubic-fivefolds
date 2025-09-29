#!/usr/bin/env sage -python
# -*- coding: utf-8 -*-
"""
Port of the Mathematica notebook "Program to find normal vectors giving maxima.nb" to SageMath/Python.

Prerequisites:
- The same directory contains Mathematica-style list files
  - normalVectors.m  (list of candidate normal vectors)
  - maximalList.m    (list of maximal sets of exponents)
  We assume these are present.
- We assume these .m files end with a Mathematica list ({ ... }) as the final expression.
  (In the notebook this was loaded via `<< file; var = %`.)

Output:
- Print the result `normal` to stdout (each element has the form [r0, r1, ..., r6, [j]]).
- Also save it as JSON to `normal.json`.

Usage example:
    sage -python find_normals.sage
It also runs with plain Python (no Sage features are used).
"""

from __future__ import annotations
import os
import re
import ast
import json
from typing import List, Sequence, Tuple

# ---- Utility: load a Mathematica-style { ... } list ----
def _strip_mathematica_comments(s: str) -> str:
    # remove (* ... *) comments (nesting not supported)
    return re.sub(r'\(\*.*?\*\)', '', s, flags=re.S)

def _extract_last_list_literal(s: str) -> str:
    """
    Extract the trailing list expression { ... } and return it as a string.
    Ignore any trailing semicolon or whitespace.
    """
    s = _strip_mathematica_comments(s).strip()
    # Skip trailing semicolons and whitespace
    i = len(s) - 1
    while i >= 0 and s[i] in ' \t\r\n;':
        i -= 1
    # Find the last '}', then search backwards for its matching '{'
    end = s.rfind('}', 0, i + 1)
    if end == -1:
        raise ValueError("Mathematica リストの '}' が見つかりません。")
    depth = 0
    start = None
    for j in range(end, -1, -1):
        c = s[j]
        if c == '}':
            depth += 1
        elif c == '{':
            depth -= 1
            if depth == 0:
                start = j
                break
    if start is None:
        raise ValueError("対応する '{' が見つかりません。")
    return s[start:end+1]

def load_mathematica_list(path: str):
    """
    Convert a Mathematica { ... } list (assumed to be a nested list of integers) into a Python list and return it.
    Example: {{1,2,3},{4,5,6}}  -> [[1,2,3],[4,5,6]]
    """
    with open(path, 'r', encoding='utf-8') as f:
        raw = f.read()
    literal = _extract_last_list_literal(raw)
    # Replace {,} with [,] and evaluate as a Python literal
    py_literal = literal.replace('{', '[').replace('}', ']')
    try:
        obj = ast.literal_eval(py_literal)
    except Exception as e:
        # Rethrow the exception with the beginning of the file attached to aid debugging
        snippet = py_literal[:2000]
        raise ValueError(f"{os.path.basename(path)} のリストを Python として解釈できません: {e}\n--- 先頭 2000 文字 ---\n{snippet}")
    return obj

# ---- Generate indexset (7 variables, total degree 3) in the same order as Mathematica's Table/Flatten ----
def make_indexset_deg3_dim7() -> List[List[int]]:
    """
    Mathematica:
        Flatten[Table[{i0,i1,i2,i3,i4,i5,i6}, {i0,0,3},...,{i6,0,3}], 6]
        // Select[weight[m]==3]
    To match the same order, enumerate nested loops with i0 the slowest and i6 the fastest, then
    keep only those whose total sum is 3.
    """
    out = []
    for i0 in range(0, 4):
      for i1 in range(0, 4):
        for i2 in range(0, 4):
          for i3 in range(0, 4):
            for i4 in range(0, 4):
              for i5 in range(0, 4):
                for i6 in range(0, 4):
                  if i0 + i1 + i2 + i3 + i4 + i5 + i6 == 3:
                      out.append([i0, i1, i2, i3, i4, i5, i6])
    return out

def dot_ge_zero(r: Sequence[int], m: Sequence[int]) -> bool:
    """Check whether r · m >= 0."""
    return sum(ri * mi for ri, mi in zip(r, m)) >= 0

def main():
    # Read files relative to the script's directory
    try:
        base_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        base_dir = os.getcwd()

    nv_path = os.path.join(base_dir, "normalVectors.m")
    ml_path = os.path.join(base_dir, "maximalList.m")

    normal_vectors = load_mathematica_list(nv_path)
    maximal_list   = load_mathematica_list(ml_path)

    # Basic shape check
    if not normal_vectors or not isinstance(normal_vectors[0], (list, tuple)):
        raise ValueError("normalVectors.m はベクトル（長さ7 の整数列）からなるリストである必要があります。")
    dim = len(normal_vectors[0])
    if dim != 7:
        raise ValueError(f"この移植コードは dim=7 を前提にしています（検出された dim={dim}）。")

    print(f"#loaded normalVectors = {len(normal_vectors)} 個")
    print(f"#loaded maximalList   = {len(maximal_list)} 個")

    indexset = make_indexset_deg3_dim7()
    print(f"#indexset (deg=3, dim=7) = {len(indexset)} 個  （期待値 84）")

    # Follow Mathematica's logic faithfully:
    # normal = {};
    # For[j=1, j<=Length[maximalList], j++,
    #   For[i=1, i<=Length[normalVectors], i++,
    #     r = normalVectors[[i]];
    #     positive[m_] := r.m >= 0;
    #     x = Select[indexset, positive];
    #     If[x == maximalList[[j]], normal = Join[{Join[r, {{j}}]}, normal]];
    #   ]
    # ]
    normal = []
    for j, target in enumerate(maximal_list, start=1):  # 1-origin
        # Ensure that target is a list of lists (just in case)
        if not isinstance(target, (list, tuple)) or not target or not isinstance(target[0], (list, tuple)):
            raise ValueError("maximalList.m の各要素は 7要素の指数ベクトルのリストである必要があります。")
        for r in normal_vectors:
            x = [m for m in indexset if dot_ge_zero(r, m)]
            if x == target:
                # Corresponds to Mathematica: Join[r, {{j}}]
                normal.append(list(r) + [[j]])

    print(f"#found normal = {len(normal)} 個")
    # Print the result
    # (also dump JSON for readability in Sage/Jupyter)
    print("normal =")
    try:
        # Compact display
        print(json.dumps(normal, ensure_ascii=False))
    except Exception:
        # Fallback
        print(normal)

    out_path = os.path.join(base_dir, "normal.json")
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(normal, f, ensure_ascii=False, indent=2)
    print(f"#saved -> {out_path}")

if __name__ == "__main__":
    main()