# FFLAS-FFPACK: Finite Field Linear Algebra Subroutines/Package

[![Build Status](https://ci.inria.fr/linbox/buildStatus/icon?job=FFLAS-FFPACK)](https://ci.inria.fr/linbox/job/FFLAS-FFPACK/)

## PURPOSE

The FFLAS-FFPACK library provides a set of basic routines for linear algebra over a finite field or the ring of integers with dense and sparse matrices.

It is inspired by the BLAS interface (Basic Linear Algebra Subprograms) and the LAPACK library for numerical linear algebra, and shares part of their design. Yet it differs in many aspects due to the specifities of computing over exact domains such as a finite fields and the rationals:
- it is generic with respect to the finite field, so as to accomodate a large variety of field sizes and implementations;
- consequently, all routines use the C++ template genericity and the library is primarily meant to be used as a source code library, to be included and compiled in the user's software.
- However, we also provide a compiled version instantiating most common routines over the most common finite fields.

## LICENSE

FFLAS-FFPACK is distributed unded the terms of the GNU LGPL v2.1 or later (see LICENSE).

## REQUIREMENTS:
- a C++ compiler supporting C++11 standard. This means g++ v4.7 or greater, clang++ v3.4 or greater, icpc v16 or greater (earlier versions of clang and icpc might also work but have not been tested)
- A BLAS library conforming to either the C or Fortran BLAS standard: OpenBLAS (recommended), or ATLAS. Make sure to use a single threaded version of the BLAS library.
- [Givaro](https://github.com/linbox-team/givaro) version at least 4.0.1, providing the implementations of the coefficient fields/rings. 

## INSTALLATION

In brief:
```./configure <options> && make && make install```

The most commonly used option include:
- `--with-blas-libs=<libs>` : to specify the arguments for the linker to find the BLAS
- `--enable-optimization` : to run configure-time optimizations

For example on a x86_64 architecture:
- Using OpenBLAS in Fedora: 
 - install the package `openblas-devel.x86_64`,
 - run `./configure --enable-optimization --with-blas-libs="-lopenblas"`
- Using OpenBLAS in Debian, Ubuntu, Mint, and all debian based distribution:
 - avoid using the distribution's package, as it is threaded by default. You need to
 - compile openblas yourself on these systems,
 - run `./configure --enable-optimization --with-blas-libs="-lopenblas"`
- Using ATLAS in Debian, Ubuntu, Mint: 
 - install the package `libatlas-dev`,
 - run `./configure --enable-optimization --with-blas-libs="-latlas -lcblas"`
- Using ATLAS in Fedora:
 - install the package `atlas-devel.x86_64`,
 - run `./configure --enable-optimization --with-blas-libs="-L/usr/lib64/atlas -lsatlas"`.

see INSTALL for further details.

## AVAILABILITY

 from [linbox-team/fflas-ffpack](https://github.com/linbox-team/fflas-ffpack)

## AUTHORS

The FFLAS-FFPACK group (see AUTHORS for a list of contributors).

## Contact and discussion

For any bug report, please file an issue on github's [issue tracker](https://github.com/linbox-team/fflas-ffpack/issues).

Please any other request, suggestion and comment to 
the discussion group [ffpack-devel](http://groups.google.com/group/ffpack-devel)
