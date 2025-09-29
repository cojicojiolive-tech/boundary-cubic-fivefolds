# Adjacency Scripts for Cases k=1,...,21

## Overview
Under the `scripts` directory you will find parallel subfolders `k1` through `k21`.
Each folder contains a single script, `adjacency.sage`, for the corresponding case k.

Running `adjacency.sage` computes and prints the normal vector(s) of the hyperplane(s) that occur whenever
I(r)_{>0} ∩ I(r_{k})_{=0}
spans a hyperplane.

- I(r)_{>0}: degree-3 exponent vectors with strictly positive weight under a 1-PS r.
- I(r_{k})_{=0}: degree-3 exponent vectors with weight zero under the fixed 1-PS r_k of case k.

The structure and behavior are the same for every k.

## Directory Layout
scripts/
├─ k1/
│  └─ adjacency.sage
├─ k2/
│  └─ adjacency.sage
⋮
└─ k21/
   └─ adjacency.sage

## Requirements
- SageMath (version 9.x or later recommended)

## How to Run
From a case directory (e.g., k1):
cd scripts/k1
sage adjacency.sage

## Output
The script prints the normal vector(s) of the relevant hyperplane(s) as integer 7-tuples.
The same usage and output format apply to every case k.
