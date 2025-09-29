This folder contains two scripts for Section 2. 
Both scripts are intended for use with SageMath.



1.  Run GIT.sage first.
GIT.sage generates two files: normalVectors.m and maximalList.m.

normalVectors.m lists all candidate normal vectors.

maximalList.m collects those half-spaces determined by the normals that are maximal with respect to inclusion (there are 23 items).



2.  Next, run find_normals.sage.
This script finds a normal vector for each of the 23 maximal items.

For 22 of the 23, it returns exactly one normal vector.

For the 10th item, it returns six normal vectors; however, these correspond to an unstable cubic fivefold and should be discarded.