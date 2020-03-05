# PythonFHE
Python DGHV scheme parallelized with Dask library (see my FHE-DGHV C++ repo for more details)

Essentially a version of this code: https://github.com/coron/fhe, *without* C optimizations or outside library dependencies and *with* parallelization and the theoretical batching implementation.

Parallelized/optimized dask version is reccommended, but a sequential (and perhaps easier to skim) version is available on the no_parallel branch.
