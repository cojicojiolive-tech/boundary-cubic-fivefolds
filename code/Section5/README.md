Section 5 — Singular Locus Analysis (Scripts)

This folder contains the scripts used to analyze the singular loci of the 21 closed-orbit normal forms in Section 5 of the paper. 
The computations use Macaulay2 (for Jacobian ideals, saturation, and set-theoretic primary decomposition) 
and Singular (for Arnold-type classification of isolated hypersurface singularities).

Folder Structure:
script/
  k1/
    components/
      components.m2
    isolated-points/
      generator.m2
      arnold_classify.sing (created by generator.m2)
  k2/
    components/
      components.m2
  ...
  k6/
    components/
      components.m2
    isolated-points/
      generator.m2
      arnold_classify.sing (created by generator.m2)
  ...
  k21/
    components/
      components.m2

There are 21 case folders, k1 … k21, corresponding to the normal forms φ_k^nf from Section 4.
Each case has a "components/" subfolder with components.m2.
For k = 1 and k = 6, there is also an "isolated-points/" subfolder, since isolated singularities occur in these cases.

Prerequisites:
- Macaulay2 (https://macaulay2.com)
- Singular (https://www.singular.uni-kl.de)

Make sure both executables are available on your PATH.

Instructions:

A. Compute the set-theoretic singular locus (all k)
--------------------------------------------------
For any case kX (e.g., k1, k7, …):
  cd script/kX/components
  M2 --script components.m2
This computes the saturated Jacobian ideal J(φ_k^nf) and decomposes the set-theoretic singular locus Sing(X_k).
The output reproduces the positive-dimensional singular locus given in Table 3 of the paper.

B. Analyze isolated singularities (only k = 1 and k = 6)
--------------------------------------------------------
1. Generate the Singular script:
   cd script/k1/isolated-points   # or script/k6/isolated-points
   M2 --script generator.m2
   This produces arnold_classify.sing.

2. Run the Arnold-type classification:
   Singular < arnold_classify.sing
This classifies the isolated point via Arnold’s normal form.
In both k = 1 and k = 6, the isolated point is of type QH(3)_19, with μ = τ = 19 and corank = 3.

C. Other cases (k ≠ 1, 6)
--------------------------
Only run components.m2. No isolated points occur, so there is no isolated-points folder.

Expected Output:
- components.m2 (Macaulay2): prints the decomposition of the singular locus, matching the entries in Table 3 of the paper.
- arnold_classify.sing (Singular, only for k = 1, 6): confirms the isolated singularity type QH(3)_19.

Reference:
This code reproduces the results in:
Y. Shibata, "The Boundary of the Moduli Space of Stable Cubic Fivefolds," Section 5 (singular loci) and Table 3.
