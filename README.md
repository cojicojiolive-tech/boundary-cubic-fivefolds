TOP-LEVEL README (Text)

Provenance
----------
This collection of scripts was used in the paper "Boundary of the Moduli Space of Stable Cubic Fivefolds".

Author
------
Yasutaka Shibata

Code Overview
-------------
Under the folder `code`, there are four subfolders that mirror the section numbering of the paper:
  - Section2
  - Section5
  - Section6
  - Section7

Each subfolder corresponds to the paper section with the same number and includes its own dedicated README with detailed instructions. This top-level README provides only a high-level overview.

Section Summaries
-----------------
• Section2
  - Purpose: Scripts that generate the list of 22 maximal 1-PS candidates used in the analysis.
  - Software: SageMath.

• Section5
  - Purpose: Computations of singularities (e.g., Jacobian ideals, saturation, eliminations).
  - Software: Macaulay2 and Singular.

• Section6
  - Purpose: Determination of adjacency relations.
  - Software: SageMath.

• Section7
  - Purpose: Machine verification that the 21 maximal lists have no inclusion relations under the SL(7) action.
  - Software: Magma. For users without access to Magma, precomputed logs are included.

Software Environment
--------------------
The following versions were used on a Mac Studio (Apple Silicon, arm64-Darwin).

• SageMath
  - Version: 10.7
  - Release Date: 2025-08-09

• Singular
  - Version: 4.4.1 (44102, 64 bit)
  - Platform: arm64-Darwin (macOS, Apple Silicon)
  - Build date: Jul 8 2025
  - Linked libraries:
    - GMP 6.3.0
    - NTL 11.5.1
    - FLINT 3.3.1

• Macaulay2
  - Version: 1.25.06
  - Platform: arm64-Darwin (macOS, Apple Silicon)

• Magma
  - Version: V2.28-25
  - Platform: macOS (arm64, Apple Silicon)
  - Host: shibatanoMac-Studio
  - Build date: Sat Sep 20 2025

Notes
-----
- All versions listed above are those installed and used on the Mac Studio (Apple Silicon, arm64-Darwin).
- Using Singular and Macaulay2, the computations in Sections 5 and 7 (e.g., saturated Jacobian ideals, Groebner-basis non-inclusion certificates) were reproduced.
- Using Magma, additional verification and algebraic manipulations were performed.
