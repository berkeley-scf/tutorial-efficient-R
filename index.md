---
title: Writing Efficient R Code
layout: default
author: Christopher Paciorek
---

## 1. This Tutorial

This tutorial covers strategies for writing efficient R code by taking advantage of the underlying structure of how R works. In addition it covers tools and strategies for timing and profiling R code.

While some of the strategies covered here are specific to R, many are built on principles that can guide your coding in other languages.

You should be able to work through this tutorial in any working R installation, including through RStudio. However, in many cases the R you are using may not be linked to a fast linear algebra package. 

This tutorial assumes you have a working knowledge of R. 

Materials for this tutorial, including the R markdown file and associated code files that were used to create this document are available on GitHub at (https://github.com/berkeley-scf/tutorial-efficient-R) in the `gh-pages` branch.
```
This tutorial by Christopher Paciorek is licensed under a Creative Commons Attribution 3.0 Unported License.


## 2. Background


In part because R is an interpreted language and in part because R
is very dynamic (objects can be modified essentially arbitrarily after
being created), R can be slow. Hadley Wickham's Advanced R book has
a [section on performance](https://adv-r.hadley.nz/perf-improve.html) that discusses this in detail. However, there
are a variety of ways that one can write efficient R code.

In general, try to make use of R's built-in functions (including matrix
operations and linear algebra), as these tend to be implemented
internally (i.e., via compiled code in C or Fortran).  Sometimes you
can figure out a trick to take your problem and transform it to make
use of the built-in functions.

Before you spend a lot of time trying to make your code go faster,
it's best to first write transparent, easy-to-read code to help avoid
bugs.  Then if it doesn't run fast enough, time the different parts of
the code (profiling) to assess where the bottlenecks are. Concentrate
your efforts on those parts of the code. Try out different
specifications, checking that the results are the same as your
original code.  And as you gain more experience, you'll get some
intuition for what approaches might improve speed, but even with
experience I find myself often surprised by what matters and what
doesn't.


## 3. Fast linear algebra

One way to speed up a variety of operations in R (sometimes by as much as an order of magnitude) is to make sure your installation of R uses an optimized BLAS (Basic Linear Algebra Subroutines). The BLAS underlies all linear algebra, including costly calculations such as matrix-matrix multiplication and matrix decompositions such as the SVD and Cholesky decomposition. Some optimized BLAS packages are:
 - Intel's *MKL*
 - *OpenBLAS*
 - AMD's *ACML*
 - *vecLib* for Macs

To use an optimized BLAS, talk to your systems adminstrator, see [Section A.3 of the R Installation and Administration Manual](https://cran.r-project.org/manuals.html), or see [these instructions to use *vecLib* BLAS on your own Mac](http://statistics.berkeley.edu/computing/blas).

Any calls to BLAS or to the LAPACK libraries that use BLAS to do higher-level linear algebra calculations will be nearly as fast as if you used C/C++ or Matlab, because R is using the compiled code from the BLAS and LAPACK libraries. 

In addition, the BLAS libraries above are parallelized -- they can use more than one core because they are 'threaded', and often will do so by default. More details in the [SCF tutorial on parallel programming](https://berkeley-scf.github.io/tutorial-parallelization). 

