Section 7 Scripts — Non-inclusion certificates for SL(7) orbits (Magma)

This folder contains the Section 7 scripts. All scripts are written for the Magma computer algebra system.

Contents
--------
- The folder scripts/O1/ contains 20 subfolders named O1-O2, …, O1-O21.
  Running the script in Section7/scripts/O1/O1-O2/ proves by machine that f1 ∉ SL(7)·f2.
- Similarly, for any folder named Ok/, each subfolder Ok-Oℓ/ contains a script 
  that proves f_k ∉ SL(7)·f_ℓ using Gröbner bases.

Output
------
Every run writes its log to a file named candidate_log in the same directory.

Prerequisites
-------------
- Magma installed and available on your PATH.

Directory layout (example)
--------------------------
Section7/
  scripts/
    O1/
      O1-O2/
        run.magma
      O1-O3/
        run.magma
      ...
      O1-O21/
        run.magma

The same pattern applies for other k (e.g., O2/, O3/, …): each Ok-Oℓ/ 
subfolder contains a run.magma that checks non-inclusion of the SL(7)-orbit 
of f_ℓ from that of f_k.

How to run
----------
From a shell:

  cd Section7/scripts/O1/O1-O2
  magma -b run.magma

- The script performs the Gröbner–basis computation (with the Rabinowitsch trick) 
  and certifies that f1 is not contained in the SL(7)-orbit of f2.
- A human-readable log is saved to candidate_log.

Run any other pair analogously; for example:

  cd ../O1-O21
  magma -b run.magma   # proves f1 ∉ SL(7)·f21

Background
----------
Each run.magma encodes the change-of-variables constraints for an SL(7) action 
and uses Gröbner bases to show that these constraints force an impossibility 
(typically det(P) = 0 in the Rabinowitsch extension), thereby proving 
non-inclusion of orbits. For details, see Section 7 of the paper.

Citation
--------
If you use these scripts in your work, please cite the accompanying paper and 
refer to Section 7 for the method and definitions.
