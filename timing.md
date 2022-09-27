---
title: Timing and profiling R code
layout: default
---

## 1. Benchmarking

*system.time* is very handy for comparing the speed of different
implementations. Here’s a basic comparison of the time to calculate the
row means of a matrix using a for loop compared to the built-in
*rowMeans* function.

``` r
n <- 10000
m <- 1000
x <- matrix(rnorm(n*m), nrow = n)
system.time({
                mns <- rep(NA, n)
                for(i in 1:n) mns[i] <- mean(x[i , ])
         })
```

    ##    user  system elapsed 
    ##   0.215   0.028   0.243

``` r
system.time(rowMeans(x))
```

    ##    user  system elapsed 
    ##   0.021   0.000   0.021

In general, *user* time gives the CPU time spent by R and *system* time
gives the CPU time spent by the kernel (the operating system) on behalf
of R. Operations that fall under system time include opening files,
doing input or output, starting other processes, etc.

To time code that runs very quickly, you should use the *microbenchmark*
package. Of course one would generally only care about accurately timing
quick calculations if a larger operation does the quick calculation very
many times. Here’s a comparison of different ways of accessing an
element of a dataframe.

``` r
library(microbenchmark)
df <- data.frame(vals = 1:3, labs = c('a','b','c'))
microbenchmark(
  df[2,1],
  df$vals[2],
  df[2, 'vals']
)
```

    ## Unit: nanoseconds
    ##           expr  min     lq    mean median     uq   max neval cld
    ##       df[2, 1] 8187 8710.5 9349.14 9049.0 9339.0 33255   100   b
    ##     df$vals[2]  602  825.5 1031.48  889.5 1064.5  9222   100  a 
    ##  df[2, "vals"] 8143 8827.5 9498.09 9089.0 9478.5 41009   100   b

**Challenge**: Think about what order of operations in the first and
third approaches above might lead to the timing result seen above.

*microbenchmark* is also a good way to compare slower calculations,
e.g., doing a crossproduct via *crossproduct()* compared to “manually”
via transpose and matrix multiply:

``` r
library(microbenchmark)
n <- 1000
x <- matrix(rnorm(n^2), n)
microbenchmark(
    t(x) %*% x,
    crossprod(x),
    times = 10)
```

    ## Unit: milliseconds
    ##          expr      min       lq     mean   median       uq      max neval cld
    ##    t(x) %*% x 30.69937 31.28786 35.04459 32.86499 38.01497 48.31874    10   b
    ##  crossprod(x) 12.43864 14.12903 14.06831 14.24277 14.42822 14.46092    10  a

**Challenge**: What is inefficient about the manual approach above?

An alternative that also automates timings and comparisons is the
*benchmark* function from the *rbenchmark* package, but there’s not
really a reason to use it when *microbenchmark* is available.

``` r
library(rbenchmark)
# speed of one calculation
benchmark(t(x) %*% x,
          crossprod(x),
          replications = 10,
          columns=c('test', 'elapsed', 'replications'))
```

    ##           test elapsed replications
    ## 2 crossprod(x)   0.130           10
    ## 1   t(x) %*% x   0.292           10

In general, it’s a good idea to repeat (replicate) your timing, as there
is some stochasticity in how fast your computer will run a piece of code
at any given moment.

You might also checkout the *tictoc* package.

## 2. Profiling

The *Rprof* function will show you how much time is spent in different
functions, which can help you pinpoint bottlenecks in your code. The
output from *Rprof* can be hard to decipher, so you will probably want
to use the *proftools* package functions, which make use of *Rprof*
under the hood.

Here’s a function that does the linear algebra to find the least squares
solution in a linear regression, assuming `x` is the matrix of
predictors, including a column for the intercept.

``` r
lr_slow <- function(y, x) {
  xtx <- t(x) %*% x
  xty <- t(x) %*% y
  inv <- solve(xtx)  ## explicit matrix inverse is slow and usually a bad idea numerically
  return(inv %*% xty)
}
```

Let’s run the function with profiling turned on.

``` r
## generate random observations and random matrix of predictors
y <- rnorm(5000)
x <- matrix(rnorm(5000*1000), nrow = 5000)

library(proftools)

pd1 <- profileExpr(lr_slow(y, x))
hotPaths(pd1)
```

    ##  path               total.pct self.pct
    ##  lr_slow            100.00     0.00   
    ##  . %*% (<text>:2)    37.50    37.50   
    ##  . %*% (<text>:3)     1.56     1.56   
    ##  . solve (<text>:4)  29.69     0.00   
    ##  . . solve.default   29.69    28.12   
    ##  . . . diag           1.56     1.56   
    ##  . t (<text>:2)       1.56     0.00   
    ##  . . t.default        1.56     1.56   
    ##  . t (<text>:3)      29.69     0.00   
    ##  . . t.default       29.69    29.69

``` r
hotPaths(pd1, value = 'time')
```

    ##  path               total.time self.time
    ##  lr_slow            1.28       0.00     
    ##  . %*% (<text>:2)   0.48       0.48     
    ##  . %*% (<text>:3)   0.02       0.02     
    ##  . solve (<text>:4) 0.38       0.00     
    ##  . . solve.default  0.38       0.36     
    ##  . . . diag         0.02       0.02     
    ##  . t (<text>:2)     0.02       0.00     
    ##  . . t.default      0.02       0.02     
    ##  . t (<text>:3)     0.38       0.00     
    ##  . . t.default      0.38       0.38

The first call to *hotPaths* shows the percentage of time spent in each
call, while the second shows the actual time.

Note the nestedness of the results. For example, essentially all the
time spent in *solve* is actually spent in *solve.default*. In this case
*solve* is just an S3 generic method that immediately calls the specific
method *solve.default*.

We can see that a lot of time is spent in doing the two crossproducts.
So let’s try using *crossprod* to make those steps faster.

``` r
lr_medium <- function(y, x) {
  xtx <- crossprod(x)
  xty <- crossprod(x, y)
  inv <- solve(xtx)  ## explicit matrix inverse is slow and usually a bad idea numerically
  return(inv %*% xty)
}                   

pd2 <- profileExpr(lr_medium(y, x))
hotPaths(pd2)
```

    ##  path                   total.pct self.pct
    ##  lr_medium              100.00     0.00   
    ##  . crossprod (<text>:2)  43.75    43.75   
    ##  . crossprod (<text>:3)   6.25     6.25   
    ##  . solve (<text>:4)      50.00     0.00   
    ##  . . solve.default       50.00    46.88   
    ##  . . . diag               3.12     3.12

``` r
hotPaths(pd2, value = 'time')
```

    ##  path                   total.time self.time
    ##  lr_medium              0.64       0.00     
    ##  . crossprod (<text>:2) 0.28       0.28     
    ##  . crossprod (<text>:3) 0.04       0.04     
    ##  . solve (<text>:4)     0.32       0.00     
    ##  . . solve.default      0.32       0.30     
    ##  . . . diag             0.02       0.02

First note that this version takes less than half the time of the
previous one. Second note that a fair amount of time is spent computing
the explicit matrix inverse using *solve*. (There’s not much we can do
to speed up the *crossprod* operations, apart from making sure we are
using a fast BLAS and potentially a parallelized BLAS.) It’s well known
that one should avoid computing the explicit inverse if one can avoid
it. Here’s a faster version that avoids it.

``` r
lr_fast <- function(y, x) {
  xtx <- crossprod(x)
  xty <- crossprod(x, y)
  U <- chol(xtx)
  tmp <- backsolve(U, xty, transpose = TRUE)
  return(backsolve(U, tmp))
}

pd3 <- profileExpr(lr_fast(y, x))
hotPaths(pd3)
```

    ##  path                   total.pct self.pct
    ##  lr_fast                100.00     0.00   
    ##  . backsolve (<text>:5)   5.56     5.56   
    ##  . chol (<text>:4)       11.11     0.00   
    ##  . . chol.default         5.56     5.56   
    ##  . . get                  5.56     0.00   
    ##  . . . lazyLoadDBfetch    5.56     5.56   
    ##  . crossprod (<text>:2)  77.78    77.78   
    ##  . crossprod (<text>:3)   5.56     5.56

``` r
hotPaths(pd3, value = 'time')
```

    ##  path                   total.time self.time
    ##  lr_fast                0.36       0.00     
    ##  . backsolve (<text>:5) 0.02       0.02     
    ##  . chol (<text>:4)      0.04       0.00     
    ##  . . chol.default       0.02       0.02     
    ##  . . get                0.02       0.00     
    ##  . . . lazyLoadDBfetch  0.02       0.02     
    ##  . crossprod (<text>:2) 0.28       0.28     
    ##  . crossprod (<text>:3) 0.02       0.02

We can see we get another speedup from that final version of the code.
(But beware my earlier caution that if comparing times between
implementations, we should have replication.)

You might also check out *profvis* for an alternative to displaying
profiling information generated by *Rprof*.

Note that *Rprof* works by sampling - every little while (the *interval*
argument) during a calculation it finds out what function R is in and
saves that information to the file given as the argument to *Rprof*. So
if you try to profile code that finishes really quickly, there’s not
enough opportunity for the sampling to represent the calculation
accurately and you may get spurious results.

> *Warning*: *Rprof* conflicts with threaded linear algebra, so you may
> need to set OMP_NUM_THREADS to 1 to disable threaded linear algebra if
> you profile code that involves linear algebra.
