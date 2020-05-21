# PythonFHE
Python DGHV scheme parallelized with Dask library (see my FHE-DGHV C++ repo for more theory details)

Essentially a version of this code: https://github.com/coron/fhe, *without* C optimizations and *with* gmpy2/smpy libraries, parallelization, and the theoretical batching implementation.

Parallelized/optimized dask version is reccommended (dask2 branch), but a relatively fast sequential (and perhaps easier to skim) version is available on the no_parallel branch.
