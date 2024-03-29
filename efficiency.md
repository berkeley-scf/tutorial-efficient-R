---
title: Writing efficient R code
layout: default
---

## 1. Pre-allocate memory

It is very inefficient to iteratively add elements to a vector, matrix,
data frame, array or list (e.g., using *c*, *cbind*, *rbind*, etc. to
add elements one at a time) in R. (Note that Python handles this sort of
thing much better.) Instead, create the full object in advance (this is
equivalent to variable initialization in compiled languages) and then
fill in the appropriate elements. The reason is that when R appends to
an existing object, it creates a new copy and as the object gets big,
most of the computation involves the repeated memory allocation to
create the new objects. Here’s an illustrative example, but of course we
would not fill a vector like this using loops because we would in
practice use vectorized calculations.

``` r
library(rbenchmark)

n <- 10000
z <- rnorm(n)

fun_append <- function(vals) {
   x <- exp(vals[1])
   n <- length(vals) 
   for(i in 2:n) x <- c(x, exp(vals[i]))
   return(x)
}
fun_prealloc <- function(vals) {
   n <- length(vals)
   x <- rep(as.numeric(NA), n)
   for(i in 1:n) x[i] <- exp(vals[i])
   return(x)
}
fun_vec <- function(vals) {
  x <- exp(vals)
  return(x)
}
benchmark(fun_append(z), fun_prealloc(z), fun_vec(z),
replications = 20, columns=c('test', 'elapsed', 'replications'))
```

    ##              test elapsed replications
    ## 1   fun_append(z)   2.016           20
    ## 2 fun_prealloc(z)   0.015           20
    ## 3      fun_vec(z)   0.003           20

It’s not necessary to use *as.numeric* above though it saves a bit of
time. **Challenge**: figure out why I have `as.numeric(NA)` and not just
`NA`. Hint: what is the type of `NA`?

In some cases, we can speed up the initialization by initializing a
vector of length one and then changing its length and/or dimension,
although in many practical circumstances this would be overkill.

For example, for matrices, start with a vector of length one, change the
length, and then change the dimensions

``` r
nr <- nc <- 2000
benchmark(
   x <- matrix(as.numeric(NA), nr, nc),
   {x <- as.numeric(NA); length(x) <- nr * nc; dim(x) <- c(nr, nc)},
replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##                                                                               test
    ## 2 {\n    x <- as.numeric(NA)\n    length(x) <- nr * nc\n    dim(x) <- c(nr, nc)\n}
    ## 1                                              x <- matrix(as.numeric(NA), nr, nc)
    ##   elapsed replications
    ## 2   0.038           10
    ## 1   0.142           10

For lists, we can do this

``` r
myList <- vector("list", length = n)
```

## 2. Vectorized calculations

One key way to write efficient R code is to take advantage of R’s
vectorized operations.

``` r
n <- 1e6
x <- rnorm(n)
benchmark(
    x2 <- x^2,
    { x2 <- as.numeric(NA)
      length(x2) <- n
      for(i in 1:n) { x2[i] <- x[i]^2 } },
    replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##                                                                                                        test
    ## 2 {\n    x2 <- as.numeric(NA)\n    length(x2) <- n\n    for (i in 1:n) {\n        x2[i] <- x[i]^2\n    }\n}
    ## 1                                                                                                 x2 <- x^2
    ##   elapsed replications
    ## 2   0.599           10
    ## 1   0.019           10

So what is different in how R handles the calculations above that
explains the huge disparity in efficiency? The vectorized calculation is
being done natively in C in a for loop. The explicit R for loop involves
executing the for loop in R with repeated calls to C code at each
iteration. This involves a lot of overhead because of the repeated
processing of the R code inside the loop. For example, in each iteration
of the loop, R is checking the types of the variables because it’s
possible that the types might change, such as in this loop:

    x <- 3
    for( i in 1:n ) {
         if(i == 7) {
              x <- 'foo'
         }
         y <- x^2
    }

You can usually get a sense for how quickly an R call will pass things
along to C or Fortran by looking at the body of the relevant function(s)
being called and looking for *.Primitive*, *.Internal*, *.C*, *.Call*,
or *.Fortran*. Let’s take a look at the code for `+`, *mean.default*,
and *chol.default*.

``` r
`+`
```

    ## function (e1, e2)  .Primitive("+")

``` r
mean.default
```

    ## function (x, trim = 0, na.rm = FALSE, ...) 
    ## {
    ##     if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
    ##         warning("argument is not numeric or logical: returning NA")
    ##         return(NA_real_)
    ##     }
    ##     if (isTRUE(na.rm)) 
    ##         x <- x[!is.na(x)]
    ##     if (!is.numeric(trim) || length(trim) != 1L) 
    ##         stop("'trim' must be numeric of length one")
    ##     n <- length(x)
    ##     if (trim > 0 && n) {
    ##         if (is.complex(x)) 
    ##             stop("trimmed means are not defined for complex data")
    ##         if (anyNA(x)) 
    ##             return(NA_real_)
    ##         if (trim >= 0.5) 
    ##             return(stats::median(x, na.rm = FALSE))
    ##         lo <- floor(n * trim) + 1
    ##         hi <- n + 1 - lo
    ##         x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
    ##     }
    ##     .Internal(mean(x))
    ## }
    ## <bytecode: 0x55fceb7ed558>
    ## <environment: namespace:base>

``` r
chol.default
```

    ## function (x, pivot = FALSE, LINPACK = FALSE, tol = -1, ...) 
    ## {
    ##     if (!missing(LINPACK)) 
    ##         stop("the LINPACK argument has been defunct since R 3.1.0")
    ##     if (is.complex(x)) 
    ##         stop("complex matrices not permitted at present")
    ##     .Internal(La_chol(as.matrix(x), pivot, tol))
    ## }
    ## <bytecode: 0x55fceb6c12e8>
    ## <environment: namespace:base>

Many R functions allow you to pass in vectors, and operate on those
vectors in vectorized fashion. So before writing a for loop, look at the
help information on the relevant function(s) to see if they operate in a
vectorized fashion. Functions might take vectors for one or more of
their arguments. Here we see that `nchar` is vectorized and that various
arguments to `substring` can be vectors.

``` r
address <- c("Four score and seven years ago our fathers brought forth",
             " on this continent, a new nation, conceived in Liberty, ",
             "and dedicated to the proposition that all men are created equal.")
nchar(address)
```

    ## [1] 56 56 64

``` r
# use a vector in the 2nd and 3rd arguments, but not the first
startIndices = seq(1, by = 3, length = nchar(address[1])/3)
startIndices
```

    ##  [1]  1  4  7 10 13 16 19 22 25 28 31 34 37 40 43 46 49 52 55

``` r
substring(address[1], startIndices, startIndices + 1)
```

    ##  [1] "Fo" "r " "co" "e " "nd" "se" "en" "ye" "rs" "ag" " o" "r " "at" "er" " b"
    ## [16] "ou" "ht" "fo" "th"

**Challenge**: Consider the chi-squared statistic involved in a test of
independence in a contingency table:

*χ*<sup>2</sup> = ∑<sub>*i*</sub>∑<sub>*j*</sub>(*y*<sub>*i*, *j*</sub>−*e*<sub>*i*, *j*</sub>)<sup>2</sup>/*e*<sub>*i*, *j*</sub>,    *e*<sub>*i*, *j*</sub> = (*y*<sub>*i*⋅</sub>*y*<sub>⋅*j*</sub>)/*y*<sub>⋅⋅</sub>

where *y*<sub>*i*⋅</sub> = ∑<sub>*j*</sub>*y*<sub>*i*, *j*</sub> and
*y*<sub>⋅*j*</sub> = ∑<sub>*i*</sub>*y*<sub>*i*, *j*</sub> and
*y*<sub>⋅⋅</sub> = ∑<sub>*i*</sub>∑<sub>*j*</sub>*y*<sub>*i*, *j*</sub>.
Write this in a vectorized way without any loops. Note that ‘vectorized’
calculations also work with matrices and arrays.

Vectorized operations can sometimes be faster than built-in functions
(note here the *ifelse* is notoriously slow), and clever vectorized
calculations even better, though sometimes the code is uglier. Here’s an
example of setting all negative values in a vector to zero.

``` r
x <- rnorm(1000000)
benchmark(
   truncx <- ifelse(x > 0, x, 0),
   {truncx <- x; truncx[x < 0] <- 0},
   truncx <- x * (x > 0),
   replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##                                            test elapsed replications
    ## 2 {\n    truncx <- x\n    truncx[x < 0] <- 0\n}   0.072           10
    ## 1                 truncx <- ifelse(x > 0, x, 0)   0.240           10
    ## 3                         truncx <- x * (x > 0)   0.042           10

Additional tips:

-   If you do need to loop over dimensions of a matrix or array, if
    possible loop over the smallest dimension and use the vectorized
    calculation on the larger dimension(s). For example if you have a
    10000 by 10 matrix, try to set up your problem so you can loop over
    the 10 columns rather than the 10000 rows.
-   In general, looping over columns is likely to be faster than looping
    over rows given R’s column-major ordering (matrices are stored in
    memory as a long array in which values in a column are adjacent to
    each other) (see more in Section 4.6 on the cache).
-   You can use direct arithmetic operations to
    add/subtract/multiply/divide a vector by each column of a matrix,
    e.g. `A*b` does element-wise multiplication of each column of *A* by
    a vector *b*. If you need to operate by row, you can do it by
    transposing the matrix.

Caution: relying on R’s recycling rule in the context of vectorized
operations, such as is done when direct-multiplying a matrix by a vector
to scale the rows relative to each other, can be dangerous as the code
is not transparent and poses greater dangers of bugs. In some cases you
may want to first write the code transparently and then compare the more
efficient code to make sure the results are the same. It’s also a good
idea to comment your code in such cases.

## 3. Using *apply* and specialized functions

Historically, another core efficiency strategy in R has been to use the
*apply* functionality (e.g., `apply`, `sapply`, `lapply`, `mapply`,
etc.).

### Some faster alternatives to `apply`

Note that even better than *apply* for calculating sums or means of
columns or rows (it also can be used for arrays) is
`rowSums`,`colSums`,`rowMeans`, and `colMeans`.

``` r
n <- 3000; x <- matrix(rnorm(n * n), nr = n)
benchmark(
   out <- apply(x, 1, mean),
   out <- rowMeans(x),
   replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##                       test elapsed replications
    ## 1 out <- apply(x, 1, mean)   1.930           10
    ## 2       out <- rowMeans(x)   0.202           10

We can ‘sweep’ out a summary statistic, such as subtracting off a mean
from each column, using *sweep*

``` r
system.time(out <- sweep(x, 2, STATS = colMeans(x), FUN = "-"))
```

    ##    user  system elapsed 
    ##   0.644   0.106   0.752

Here’s a trick for doing the sweep based on vectorized calculations,
remembering that if we subtract a vector from a matrix, it subtracts
each element of the vector from all the elements in the corresponding
ROW. Hence the need to transpose twice.

``` r
system.time(out2 <- t(t(x) - colMeans(x)))
```

    ##    user  system elapsed 
    ##   0.188   0.043   0.232

``` r
identical(out, out2)
```

    ## [1] TRUE

### Are *apply*, *lapply*, *sapply*, etc. faster than loops?

Using *apply* with matrices and versions of *apply* with lists or
vectors (e.g., `lapply`, `sapply`) may or may not be faster than looping
but generally produces cleaner code.

Whether looping and use of apply variants is slow will depend in part on
whether a substantial part of the work is in the overhead involved in
the looping or in the time required by the function evaluation on each
of the elements. If you’re worried about speed, it’s a good idea to
benchmark the *apply* variant against looping.

Here’s an example where *apply* is not faster than a loop. Similar
examples can be constructed where *lapply* or *sapply* are not faster
than writing a loop.

``` r
n <- 500000; nr <- 10000; nCalcs <- n/nr
mat <- matrix(rnorm(n), nrow = nr)
times <- 1:nr
system.time(
  out1 <- apply(mat, 2, function(vec) {
                         mod = lm(vec ~ times)
                         return(mod$coef[2])
                     })) 
```

    ##    user  system elapsed 
    ##   0.100   0.012   0.113

``` r
system.time({
                out2 <- rep(NA, nCalcs)
                for(i in 1:nCalcs){
                    out2[i] = lm(mat[ , i] ~ times)$coef[2]
                }
            }) 
```

    ##    user  system elapsed 
    ##   0.076   0.000   0.075

And here’s an example, where (unlike the previous example) the core
computation is very fast, so we might expect the overhead of looping to
be important. I believe that in old versions of R the *sapply* in this
example was faster than looping in R, but that doesn’t seem to be the
case currently. I think this may be related to various somewhat recent
improvements in R’s handling of loops, possibly including the use of the
byte compiler.

``` r
z <- rnorm(1e6)
fun_loop <- function(vals) {
    x <- as.numeric(NA)
    n <- length(vals)
    length(x) <- n
    for(i in 1:n) x[i] <- exp(vals[i])
    return(x)
}

fun_sapply <- function(vals) {
    x <- sapply(vals, exp)
    return(x)
}

fun_vec <- function(vals) {
    x <- exp(vals)
    return(x)
}

benchmark(fun_loop(z), fun_sapply(z), fun_vec(z),
replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##            test elapsed replications
    ## 1   fun_loop(z)   0.744           10
    ## 2 fun_sapply(z)   4.118           10
    ## 3    fun_vec(z)   0.065           10

You’ll notice if you look at the R code for *lapply* (*sapply* just
calls *lapply*) that it calls directly out to C code, so the for loop is
executed in compiled code. However, the code being executed at each
iteration is still R code, so there is still all the overhead of the R
interpreter.

``` r
print(lapply)
```

    ## function (X, FUN, ...) 
    ## {
    ##     FUN <- match.fun(FUN)
    ##     if (!is.vector(X) || is.object(X)) 
    ##         X <- as.list(X)
    ##     .Internal(lapply(X, FUN))
    ## }
    ## <bytecode: 0x55fce96a15b8>
    ## <environment: namespace:base>

## 4. Matrix algebra efficiency

Often calculations that are not explicitly linear algebra calculations
can be done as matrix algebra. If our R installation has a fast (and
possibly parallelized) BLAS, this allows our calculation to take
advantage of it.

For example, we can sum the rows of a matrix by multiplying by a vector
of ones. Given the extra computation involved in actually multiplying
each number by one, it’s surprising that this is faster than using R’s
heavily optimized *rowSums* function.

``` r
mat <- matrix(rnorm(500*500), 500)
ones <- rep(1, ncol(mat))
benchmark(apply(mat, 1, sum),
    mat %*% ones,
    rowSums(mat),
    replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##                 test elapsed replications
    ## 1 apply(mat, 1, sum)   0.031           10
    ## 2       mat %*% ones   0.003           10
    ## 3       rowSums(mat)   0.010           10

On the other hand, big matrix operations can be slow.

> **Challenge**: Suppose you want a new matrix that computes the
> differences between successive columns of a matrix of arbitrary size.
> How would you do this as matrix algebra operations? It’s possible to
> write it as multiplying the matrix by another matrix that contains 0s,
> 1s, and -1s in appropriate places. Here it turns out that the *for*
> loop is much faster than matrix multiplication. However, there is a
> way to do it faster as matrix direct subtraction.

### Order of operations and efficiency

When doing matrix algebra, the order in which you do operations can be
critical for efficiency. How should I order the following calculation?

``` r
n <- 5000
A <- matrix(rnorm(5000 * 5000), 5000)
B <- matrix(rnorm(5000 * 5000), 5000)
x <- rnorm(5000)
system.time(
  res1 <- A %*% B %*% x
)
```

    ##    user  system elapsed 
    ##  12.771   2.625   2.347

``` r
system.time(
  res2 <- A %*% (B %*% x)
)
```

    ##    user  system elapsed 
    ##   0.279   0.088   0.060

Why is the second order much faster?

### Avoiding unnecessary operations

We can use the matrix direct product (i.e., `A*B`) to do some
manipulations much more quickly than using matrix multiplication.
**Challenge**: How can I use the direct product to find the trace of a
matrix, `trace(XY)`?

Finally, when working with diagonal matrices, you can generally get much
faster results by being smart. The following operations: `X+D`, `DX`,
`XD` are mathematically the sum of two matrices and products of two
matrices. But we can do the computation without using two full matrices.
**Challenge**: How?

``` r
n <- 1000
X <- matrix(rnorm(n^2), n) 
diagvals <- rnorm(n)
D = diag(diagvals)
# the following lines are very inefficient
summedMat <- X + D
prodMat1 <- D %*% X
prodMat2 <- X %*% D
# How can we do each of those operations much more quickly?
```

More generally, sparse matrices and structured matrices (such as block
diagonal matrices) can generally be worked with MUCH more efficiently
than treating them as arbitrary matrices. The R packages *spam* (for
arbitrary sparse matrices), *bdsmatrix* (for block-diagonal matrices),
and *Matrix* (for a variety of sparse matrix types) can help, as can
specialized code available in other languages, such as C and Fortran
packages.

## 5. Fast mapping/lookup tables

Sometimes you need to map between two vectors. E.g.,
*y*<sub>*i*</sub> ∼ 𝒩(*μ*<sub>*j*\[*i*\]</sub>,*σ*<sup>2</sup>) is a
basic ANOVA type structure, where multiple observations are associated
with a common mean, *μ*<sub>*j*</sub>, via the `j[i]` mapping.

How can we quickly look up the mean associated with each observation? A
good strategy is to create a vector, `grp`, that gives a numeric mapping
of the observations to their cluster, playing the role of `j[i]` above.
Then you can access the *μ* value relevant for each observation as:
`mus[grp]`. This requires that `grp` correctly map to the right elements
of `mus`.

The *match* function can help in creating numeric indices that can then
be used for lookups. Here’s how you would create an index vector, *grp*,
if it doesn’t already exist.

``` r
df <- data.frame(
    id = 1:5,
    clusterLabel = c('C', 'B', 'B', 'A', 'C')) 
info <- data.frame(
    grade = c('A', 'B', 'C'),
    numGrade = c(95, 85, 75),
    fail = c(FALSE, FALSE, TRUE) )
grp <- match(df$clusterLabel, info$grade) 
df$numGrade <- info$numGrade[grp]
df
```

    ##   id clusterLabel numGrade
    ## 1  1            C       75
    ## 2  2            B       85
    ## 3  3            B       85
    ## 4  4            A       95
    ## 5  5            C       75

### Lookup by name versus index

In the example above we looked up the `mu` values based on `grp`, which
supplies the needed indexes as numeric indexes.

R also allows you to look up elements of vector by name, as illustrated
here by rearranging the code above a bit:

``` r
info2 <- info$numGrade
names(info2) <- info$grade
info2
```

    ##  A  B  C 
    ## 95 85 75

``` r
info2[df$clusterLabel]
```

    ##  C  B  B  A  C 
    ## 75 85 85 95 75

You can do similar things in terms of looking up by name with dimension
names of matrices/arrays, row and column names of dataframes, and named
lists.

However, looking things up by name can be slow relative to looking up by
index. Here’s a toy example where we have a vector or list with 1000
elements and the character names of the elements are just the character
versions of the indices of the elements.

``` r
library(microbenchmark)

n <- 1000
x <- 1:n
xL <- as.list(x)
nms <- paste0("var", as.character(x))
names(x) <- nms
names(xL) <- nms
x[1:3]
```

    ## var1 var2 var3 
    ##    1    2    3

``` r
xL[1:3]
```

    ## $var1
    ## [1] 1
    ## 
    ## $var2
    ## [1] 2
    ## 
    ## $var3
    ## [1] 3

``` r
microbenchmark(
    x[500],  # index lookup in vector
    x["var500"], # name lookup in vector
    xL[[500]], # index lookup in list
    xL[["var500"]]) # name lookup in list
```

    ## Unit: nanoseconds
    ##            expr  min     lq    mean median     uq   max neval cld
    ##          x[500]  267  285.0  309.25  301.0  310.0   589   100 a  
    ##     x["var500"] 2153 2195.5 2356.37 2221.0 2239.5 15519   100  b 
    ##       xL[[500]]  122  146.0  160.12  155.5  166.5   364   100 a  
    ##  xL[["var500"]] 3100 3129.0 3208.02 3145.0 3167.0  8999   100   c

Lookup by name is slow because R needs to scan through the objects one
by one until it finds the one with the name it is looking for. In
contrast, to look up by index, R can just go directly to the position of
interest.

Side note: there is a lot of variability across the 100 replications
shown above. This might have to do with cache effects (see next
section).

In contrast, we can look up by name in an environment very quickly,
because environments in R use hashing, which allows for fast lookup that
does not require scanning through all of the names in the environment.
In fact, this is how R itself looks for values when you specify
variables in R code, because the global environment, function frames,
and package namespaces are all environments.

``` r
xEnv <- as.environment(xL)  # convert from a named list
xEnv$var500  
```

    ## [1] 500

``` r
microbenchmark(
  x[500],
  xL[[500]],
  xEnv[["var500"]],
  xEnv$var500
)
```

    ## Unit: nanoseconds
    ##              expr min    lq   mean median    uq  max neval cld
    ##            x[500] 267 278.0 446.92  287.0 468.5 3732   100   c
    ##         xL[[500]] 123 132.5 204.50  140.0 191.5 3028   100 ab 
    ##  xEnv[["var500"]] 119 129.0 180.70  137.5 193.5  694   100 a  
    ##       xEnv$var500 171 187.0 319.57  197.5 278.0 4687   100  b

## 6. Cache-aware programming

In addition to main memory (what we usually mean when we talk about
RAM), computers also have memory caches, which are small amounts of fast
memory that can be accessed very quickly by the processor. For example
your computer might have L1, L2, and L3 caches, with L1 the smallest and
fastest and L3 the largest and slowest. The idea is to try to have the
data that is most used by the processor in the cache.

If the next piece of data needed for computation is available in the
cache, this is a *cache hit* and the data can be accessed very quickly.
However, if the data is not available in the cache, this is a *cache
miss* and the speed of access will be a lot slower. *Cache-aware
programming* involves writing your code to minimize cache misses.
Generally when data is read from memory it will be read in chunks, so
values that are contiguous will be read together.

How does this inform one’s programming? For example, if you have a
matrix of values stored in column-major order, computing on a column
will be a lot faster than computing on a row, because the column can be
read into the cache from main memory and then accessed in the cache. In
contrast, if the matrix is large and therefore won’t fit in the cache,
when you access the values of a row, you’ll have to go to main memory
repeatedly to get the values for the row because the values are not
stored contiguously.

There’s a nice example of the importance of the cache at [the bottom of
this blog
post](https://wrathematics.github.io/2016/10/28/comparing-symmetric-eigenvalue-performance/).

If you know the size of the cache, you can try to design your code so
that in a given part of your code you access data structures that will
fit in the cache. This sort of thing is generally more relevant if
you’re coding in a language like C. But it can matter sometimes in R
too. Here’s an example:

``` r
nr <- 800000
nc <- 100
## large matrix that won't fit in cache
A <- matrix(rnorm(nr * nc), nrow = nr)
tA <- t(A)
benchmark(
    apply(A, 2, mean),   ## operate by column
    apply(tA, 1, mean),  ## exact same calculation, but by row
    replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##                 test elapsed replications
    ## 1  apply(A, 2, mean)   9.594           10
    ## 2 apply(tA, 1, mean)  21.462           10

Now let’s compare things when we make the matrix small enough that it
fits in the cache. In this case it fits into the largest (L3) cache but
not the smaller (L2 and L1) caches, and we see that the difference in
speed disappears.

``` r
nr <- 800
nc <- 100
## small matrix that should fit in cache
A <- matrix(rnorm(nr * nc), nrow = nr)
## Yep, the size is less than the L3 cache:
object.size(A)
```

    ## 640216 bytes

``` r
memuse::Sys.cachesize()
```

    ## L1I:   32.000 KiB 
    ## L1D:   32.000 KiB 
    ## L2:   256.000 KiB 
    ## L3:     8.000 MiB

``` r
tA <- t(A)
benchmark(apply(A, 2, mean),  ## by column
          apply(tA, 1, mean),  ## by row
  replications = 10, columns=c('test', 'elapsed', 'replications'))
```

    ##                 test elapsed replications
    ## 1  apply(A, 2, mean)   0.012           10
    ## 2 apply(tA, 1, mean)   0.009           10
