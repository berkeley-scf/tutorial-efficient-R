---
layout: default
title: Writing efficient R code
toc-title: Table of contents
---

## 1. Pre-allocate memory

It is very inefficient to iteratively add elements to a vector, matrix,
data frame, array or list (e.g., using *c*, *cbind*, *rbind*, etc. to
add elements one at a time). Instead, create the full object in advance
(this is equivalent to variable initialization in compiled languages)
and then fill in the appropriate elements. The reason is that when R
appends to an existing object, it creates a new copy and as the object
gets big, most of the computation involves the repeated memory
allocation to create the new objects. Here's an illustrative example,
but of course we would not fill a vector like this using loops because
we would in practice use vectorized calculations.

::: {.cell hash="efficiency_cache/markdown/preallocate_a20407768fc348fc5e867785ddc8ce20"}
``` {.r .cell-code}
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

::: {.cell-output .cell-output-stdout}
                 test elapsed replications
    1   fun_append(z)   2.026           20
    2 fun_prealloc(z)   0.014           20
    3      fun_vec(z)   0.003           20
:::
:::

It's not necessary to use *as.numeric* above though it saves a bit of
time. **Challenge**: figure out why I have `as.numeric(NA)` and not just
`NA`. Hint: what is the type of `NA`?

In some cases, we can speed up the initialization by initializing a
vector of length one and then changing its length and/or dimension,
although in many practical circumstances this would be overkill.

For example, for matrices, start with a vector of length one, change the
length, and then change the dimensions

::: cell
``` {.r .cell-code}
nr <- nc <- 2000
benchmark(
   x <- matrix(as.numeric(NA), nr, nc),
   {x <- as.numeric(NA); length(x) <- nr * nc; dim(x) <- c(nr, nc)},
replications = 10, columns=c('test', 'elapsed', 'replications'))
```

::: {.cell-output .cell-output-stdout}
                                                                                  test
    2 {\n    x <- as.numeric(NA)\n    length(x) <- nr * nc\n    dim(x) <- c(nr, nc)\n}
    1                                              x <- matrix(as.numeric(NA), nr, nc)
      elapsed replications
    2   0.039           10
    1   0.159           10
:::
:::

For lists, we can do this

::: cell
``` {.r .cell-code}
myList <- vector("list", length = n)
```
:::

## 2. Vectorized calculations

One key way to write efficient R code is to take advantage of R's
vectorized operations.

::: {.cell hash="efficiency_cache/markdown/vectorize_7c658677f8cfd2d1f9f1342ade9855a4"}
``` {.r .cell-code}
n <- 1e6
x <- rnorm(n)
benchmark(
    x2 <- x^2,
    { x2 <- as.numeric(NA)
      length(x2) <- n
      for(i in 1:n) { x2[i] <- x[i]^2 } },
    replications = 10, columns=c('test', 'elapsed', 'replications'))
```

::: {.cell-output .cell-output-stdout}
                                                                                                           test
    2 {\n    x2 <- as.numeric(NA)\n    length(x2) <- n\n    for (i in 1:n) {\n        x2[i] <- x[i]^2\n    }\n}
    1                                                                                                 x2 <- x^2
      elapsed replications
    2   0.752           10
    1   0.028           10
:::
:::

So what is different in how R handles the calculations above that
explains the huge disparity in efficiency? The vectorized calculation is
being done natively in C in a for loop. The explicit R for loop involves
executing the for loop in R with repeated calls to C code at each
iteration. This involves a lot of overhead because of the repeated
processing of the R code inside the loop. For example, in each iteration
of the loop, R is checking the types of the variables because it's
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
or *.Fortran*. Let's take a look at the code for `+`, *mean.default*,
and *chol.default*.

::: cell
``` {.r .cell-code}
`+`
```

::: {.cell-output .cell-output-stdout}
    function (e1, e2)  .Primitive("+")
:::

``` {.r .cell-code}
mean.default
```

::: {.cell-output .cell-output-stdout}
    function (x, trim = 0, na.rm = FALSE, ...) 
    {
        if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
            warning("argument is not numeric or logical: returning NA")
            return(NA_real_)
        }
        if (isTRUE(na.rm)) 
            x <- x[!is.na(x)]
        if (!is.numeric(trim) || length(trim) != 1L) 
            stop("'trim' must be numeric of length one")
        n <- length(x)
        if (trim > 0 && n) {
            if (is.complex(x)) 
                stop("trimmed means are not defined for complex data")
            if (anyNA(x)) 
                return(NA_real_)
            if (trim >= 0.5) 
                return(stats::median(x, na.rm = FALSE))
            lo <- floor(n * trim) + 1
            hi <- n + 1 - lo
            x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
        }
        .Internal(mean(x))
    }
    <bytecode: 0x56513781d148>
    <environment: namespace:base>
:::

``` {.r .cell-code}
chol.default
```

::: {.cell-output .cell-output-stdout}
    function (x, pivot = FALSE, LINPACK = FALSE, tol = -1, ...) 
    {
        if (!missing(LINPACK)) 
            stop("the LINPACK argument has been defunct since R 3.1.0")
        if (is.complex(x)) 
            stop("complex matrices not permitted at present")
        .Internal(La_chol(as.matrix(x), pivot, tol))
    }
    <bytecode: 0x5651397a2f28>
    <environment: namespace:base>
:::
:::

Many R functions allow you to pass in vectors, and operate on those
vectors in vectorized fashion. So before writing a for loop, look at the
help information on the relevant function(s) to see if they operate in a
vectorized fashion. Functions might take vectors for one or more of
their arguments. Here we see that `nchar` is vectorized and that various
arguments to `substring` can be vectors.

::: cell
``` {.r .cell-code}
address <- c("Four score and seven years ago our fathers brought forth",
             " on this continent, a new nation, conceived in Liberty, ",
             "and dedicated to the proposition that all men are created equal.")
nchar(address)
```

::: {.cell-output .cell-output-stdout}
    [1] 56 56 64
:::

``` {.r .cell-code}
# use a vector in the 2nd and 3rd arguments, but not the first
startIndices = seq(1, by = 3, length = nchar(address[1])/3)
startIndices
```

::: {.cell-output .cell-output-stdout}
     [1]  1  4  7 10 13 16 19 22 25 28 31 34 37 40 43 46 49 52 55
:::

``` {.r .cell-code}
substring(address[1], startIndices, startIndices + 1)
```

::: {.cell-output .cell-output-stdout}
     [1] "Fo" "r " "co" "e " "nd" "se" "en" "ye" "rs" "ag" " o" "r " "at" "er" " b"
    [16] "ou" "ht" "fo" "th"
:::
:::

**Challenge**: Consider the chi-squared statistic involved in a test of
independence in a contingency table:

\[
`\chi`{=tex}\^{2}=`\sum`{=tex}*{i}`\sum`{=tex}*{j}`\frac{(y_{ij}-e_{ij})^{2}}{e_{ij}}`{=tex},,,,,
e\_{ij}=`\frac{y_{i\cdot}y_{\cdot j}}{y_{\cdot\cdot}}`{=tex} \]

where $y_{i\cdot}=\sum_{j}y_{ij}$ and $y_{\cdot j} = \sum_{i} y_{ij}$.
Write this in a vectorized way without any loops. Note that 'vectorized'
calculations also work with matrices and arrays.

Vectorized operations can sometimes be faster than built-in functions
(note here the *ifelse* is notoriously slow), and clever vectorized
calculations even better, though sometimes the code is uglier. Here's an
example of setting all negative values in a vector to zero.

::: {.cell hash="efficiency_cache/markdown/vec-tricks_2147592b9a565c3acf200e3fad74d719"}
``` {.r .cell-code}
x <- rnorm(1000000)
benchmark(
   truncx <- ifelse(x > 0, x, 0),
   {truncx <- x; truncx[x < 0] <- 0},
   truncx <- x * (x > 0),
   replications = 10, columns=c('test', 'elapsed', 'replications'))
```

::: {.cell-output .cell-output-stdout}
                                               test elapsed replications
    2 {\n    truncx <- x\n    truncx[x < 0] <- 0\n}   0.071           10
    1                 truncx <- ifelse(x > 0, x, 0)   0.278           10
    3                         truncx <- x * (x > 0)   0.044           10
:::
:::

Additional tips: - If you do need to loop over dimensions of a matrix or
array, if possible loop over the smallest dimension and use the
vectorized calculation on the larger dimension(s). For example if you
have a 10000 by 10 matrix, try to set up your problem so you can loop
over the 10 columns rather than the 10000 rows. - In general, looping
over columns is likely to be faster than looping over rows given R's
column-major ordering (matrices are stored in memory as a long array in
which values in a column are adjacent to each other) (see more in
Section 4.6 on the cache). - You can use direct arithmetic operations to
add/subtract/multiply/divide a vector by each column of a matrix,
e.g. `A*b` does element-wise multiplication of each column of *A* by a
vector *b*. If you need to operate by row, you can do it by transposing
the matrix.

Caution: relying on R's recycling rule in the context of vectorized
operations, such as is done when direct-multiplying a matrix by a vector
to scale the rows relative to each other, can be dangerous as the code
is not transparent and poses greater dangers of bugs. In some cases you
may want to first write the code transparently and then compare the more
efficient code to make sure the results are the same. It's also a good
idea to comment your code in such cases.

## 3. Using *apply* and specialized functions

Historically, another core efficiency strategy in R has been to use the
*apply* functionality (e.g., `apply`, `sapply`, `lapply`, `mapply`,
etc.).

### Some faster alternatives to `apply`

Note that even better than *apply* for calculating sums or means of
columns or rows (it also can be used for arrays) is
{row,col}{Sums,Means}.

::: {.cell hash="efficiency_cache/markdown/apply_f8a3e3e8e3a32629daf41a1832991a51"}
``` {.r .cell-code}
n <- 3000; x <- matrix(rnorm(n * n), nr = n)
benchmark(
   out <- apply(x, 1, mean),
   out <- rowMeans(x),
   replications = 10, columns=c('test', 'elapsed', 'replications'))
```

::: {.cell-output .cell-output-stdout}
                          test elapsed replications
    1 out <- apply(x, 1, mean)   2.094           10
    2       out <- rowMeans(x)   0.200           10
:::
:::

We can 'sweep' out a summary statistic, such as subtracting off a mean
from each column, using *sweep*

::: cell
``` {.r .cell-code}
system.time(out <- sweep(x, 2, STATS = colMeans(x), FUN = "-"))
```

::: {.cell-output .cell-output-stdout}
       user  system elapsed 
      0.143   0.024   0.167 
:::
:::

Here's a trick for doing the sweep based on vectorized calculations,
remembering that if we subtract a vector from a matrix, it subtracts
each element of the vector from all the elements in the corresponding
ROW. Hence the need to transpose twice.

::: cell
``` {.r .cell-code}
system.time(out2 <- t(t(x) - colMeans(x)))
```

::: {.cell-output .cell-output-stdout}
       user  system elapsed 
      0.269   0.044   0.314 
:::

``` {.r .cell-code}
identical(out, out2)
```

::: {.cell-output .cell-output-stdout}
    [1] TRUE
:::
:::

### Are *apply*, *lapply*, *sapply*, etc. faster than loops?

Using *apply* with matrices and versions of *apply* with lists or
vectors (e.g., `lapply`, `sapply`) may or may not be faster than looping
but generally produces cleaner code.

Whether looping and use of apply variants is slow will depend in part on
whether a substantial part of the work is in the overhead involved in
the looping or in the time required by the function evaluation on each
of the elements. If you're worried about speed, it's a good idea to
benchmark the *apply* variant against looping.

Here's an example where *apply* is not faster than a loop. Similar
examples can be constructed where *lapply* or *sapply* are not faster
than writing a loop.

::: cell
``` {.r .cell-code}
n <- 500000; nr <- 10000; nCalcs <- n/nr
mat <- matrix(rnorm(n), nrow = nr)
times <- 1:nr
system.time(
  out1 <- apply(mat, 2, function(vec) {
                         mod = lm(vec ~ times)
                         return(mod$coef[2])
                     })) 
```

::: {.cell-output .cell-output-stdout}
       user  system elapsed 
      0.082   0.000   0.082 
:::

``` {.r .cell-code}
system.time({
                out2 <- rep(NA, nCalcs)
                for(i in 1:nCalcs){
                    out2[i] = lm(mat[ , i] ~ times)$coef[2]
                }
            }) 
```

::: {.cell-output .cell-output-stdout}
       user  system elapsed 
      0.077   0.000   0.077 
:::
:::

And here's an example, where (unlike the previous example) the core
computation is very fast, so we might expect the overhead of looping to
be important. I believe that in old versions of R the *sapply* in this
example was faster than looping in R, but that doesn't seem to be the
case currently. I think this may be related to various somewhat recent
improvements in R's handling of loops, possibly including the use of the
byte compiler.

::: {.cell hash="efficiency_cache/markdown/apply-vs-for-part2_ba67656b594585c560785e4e1101339f"}
``` {.r .cell-code}
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

::: {.cell-output .cell-output-stdout}
               test elapsed replications
    1   fun_loop(z)   0.719           10
    2 fun_sapply(z)   4.063           10
    3    fun_vec(z)   0.062           10
:::
:::

You'll notice if you look at the R code for *lapply* (*sapply* just
calls *lapply*) that it calls directly out to C code, so the for loop is
executed in compiled code. However, the code being executed at each
iteration is still R code, so there is still all the overhead of the R
interpreter.

::: cell
``` {.r .cell-code}
print(lapply)
```

::: {.cell-output .cell-output-stdout}
    function (X, FUN, ...) 
    {
        FUN <- match.fun(FUN)
        if (!is.vector(X) || is.object(X)) 
            X <- as.list(X)
        .Internal(lapply(X, FUN))
    }
    <bytecode: 0x565133605438>
    <environment: namespace:base>
:::
:::

## 4. Matrix algebra efficiency

Often calculations that are not explicitly linear algebra calculations
can be done as matrix algebra. If our R installation has a fast (and
possibly parallelized) BLAS, this allows our calculation to take
advantage of it.

For example, we can sum the rows of a matrix by multiplying by a vector
of ones. Given the extra computation involved in actually multiplying
each number by one, it's surprising that this is faster than using R's
heavily optimized *rowSums* function.

::: cell
``` {.r .cell-code}
mat <- matrix(rnorm(500*500), 500)
benchmark(apply(mat, 1, sum),
    mat %*% rep(1, ncol(mat)),
    rowSums(mat),
    replications = 10, columns=c('test', 'elapsed', 'replications'))
```

::: {.cell-output .cell-output-stdout}
                           test elapsed replications
    1        apply(mat, 1, sum)   0.027           10
    2 mat %*% rep(1, ncol(mat))   0.003           10
    3              rowSums(mat)   0.010           10
:::
:::

On the other hand, big matrix operations can be slow. **Challenge**:
Suppose you want a new matrix that computes the differences between
successive columns of a matrix of arbitrary size. How would you do this
as matrix algebra operations? It's possible to write it as multiplying
the matrix by another matrix that contains 0s, 1s, and -1s in
appropriate places. Here it turns out that the *for* loop is much faster
than matrix multiplication. However, there is a way to do it faster as
matrix direct subtraction.

### Order of operations and efficiency

When doing matrix algebra, the order in which you do operations can be
critical for efficiency. How should I order the following calculation?

::: {.cell hash="efficiency_cache/markdown/linalg-order_703218449b96601964de98bd53dd968c"}
``` {.r .cell-code}
n <- 5000
A <- matrix(rnorm(5000 * 5000), 5000)
B <- matrix(rnorm(5000 * 5000), 5000)
x <- rnorm(5000)
system.time(
  res1 <- A %*% B %*% x
)
```

::: {.cell-output .cell-output-stdout}
       user  system elapsed 
     11.809   1.210   1.754 
:::

``` {.r .cell-code}
system.time(
  res2 <- A %*% (B %*% x)
)
```

::: {.cell-output .cell-output-stdout}
       user  system elapsed 
      0.314   0.053   0.061 
:::
:::

Why is the second order much faster?

### Avoiding unnecessary operations

We can use the matrix direct product (i.e., `A*B`) to do some
manipulations much more quickly than using matrix multiplication.
**Challenge**: How can I use the direct product to find the trace of a
matrix, $XY$?

Finally, when working with diagonal matrices, you can generally get much
faster results by being smart. The following operations: $X+D$, $DX$,
$XD$ are mathematically the sum of two matrices and products of two
matrices. But we can do the computation without using two full matrices.
**Challenge**: How?

::: {.cell hash="efficiency_cache/markdown/diag_65d46660cb0e7bebabd81248d29c7726"}
``` {.r .cell-code}
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
:::

More generally, sparse matrices and structured matrices (such as block
diagonal matrices) can generally be worked with MUCH more efficiently
than treating them as arbitrary matrices. The R packages *spam* (for
arbitrary sparse matrices), *bdsmatrix* (for block-diagonal matrices),
and *Matrix* (for a variety of sparse matrix types) can help, as can
specialized code available in other languages, such as C and Fortran
packages.

## 5. Fast mapping/lookup tables

Sometimes you need to map between two vectors. E.g.,
$y_{i}\sim\mathcal{N}(\mu_{j[i]},\sigma^{2})$ is a basic ANOVA type
structure, where multiple observations are associated with a common
mean, $\mu_j$, via the `j[i]` mapping.

How can we quickly look up the mean associated with each observation? A
good strategy is to create a vector, *grp*, that gives a numeric mapping
of the observations to their cluster, playing the role of `j[i]` above.
Then you can access the $\mu$ value relevant for each observation as:
`mus[grp]`. This requires that *grp* correctly map to the right elements
of *mus*.

The *match* function can help in creating numeric indices that can then
be used for lookups. Here's how you would create an index vector, *grp*,
if it doesn't already exist.

::: cell
``` {.r .cell-code}
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

::: {.cell-output .cell-output-stdout}
      id clusterLabel numGrade
    1  1            C       75
    2  2            B       85
    3  3            B       85
    4  4            A       95
    5  5            C       75
:::
:::

### Lookup by name versus index

In the example above we looked up the `mu` values based on `grp`, which
supplies the needed indexes as numeric indexes.

R also allows you to look up elements of vector by name, as illustrated
here by rearranging the code above a bit:

::: cell
``` {.r .cell-code}
info2 <- info$numGrade
names(info2) <- info$grade
info2
```

::: {.cell-output .cell-output-stdout}
     A  B  C 
    95 85 75 
:::

``` {.r .cell-code}
info2[df$clusterLabel]
```

::: {.cell-output .cell-output-stdout}
     C  B  B  A  C 
    75 85 85 95 75 
:::
:::

You can do similar things in terms of looking up by name with dimension
names of matrices/arrays, row and column names of dataframes, and named
lists.

However, looking things up by name can be slow relative to looking up by
index. Here's a toy example where we have a vector or list with 1000
elements and the character names of the elements are just the character
versions of the indices of the elements.

::: cell
``` {.r .cell-code}
library(microbenchmark)

n <- 1000
x <- 1:n
xL <- as.list(x)
nms <- paste0("var", as.character(x))
names(x) <- nms
names(xL) <- nms
x[1:3]
```

::: {.cell-output .cell-output-stdout}
    var1 var2 var3 
       1    2    3 
:::

``` {.r .cell-code}
xL[1:3]
```

::: {.cell-output .cell-output-stdout}
    $var1
    [1] 1

    $var2
    [1] 2

    $var3
    [1] 3
:::

``` {.r .cell-code}
microbenchmark(
    x[500],  # index lookup in vector
    x["var500"], # name lookup in vector
    xL[[500]], # index lookup in list
    xL[["var500"]]) # name lookup in list
```

::: {.cell-output .cell-output-stdout}
    Unit: nanoseconds
               expr  min     lq    mean median     uq  max neval  cld
             x[500]  271  297.0  344.25  309.5  324.5 3329   100  b  
        x["var500"] 2163 2191.0 2273.04 2207.5 2236.5 7394   100   c 
          xL[[500]]  124  152.0  188.67  162.0  175.0 2354   100 a   
     xL[["var500"]] 3061 3188.5 3241.54 3244.5 3289.0 3819   100    d
:::
:::

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

::: {.cell hash="efficiency_cache/markdown/env-lookup_52016cacc00205b972416156e5e6e063"}
``` {.r .cell-code}
xEnv <- as.environment(xL)  # convert from a named list
xEnv$var500  
```

::: {.cell-output .cell-output-stdout}
    [1] 500
:::

``` {.r .cell-code}
microbenchmark(
  x[500],
  xL[[500]],
  xEnv[["var500"]],
  xEnv$var500
)
```

::: {.cell-output .cell-output-stdout}
    Unit: nanoseconds
                 expr min    lq   mean median    uq  max neval cld
               x[500] 259 274.0 360.31  282.0 345.0 3967   100   b
            xL[[500]] 124 137.0 181.97  143.5 193.0 1934   100  a 
     xEnv[["var500"]] 111 120.5 166.00  129.0 146.5 2344   100  a 
          xEnv$var500 174 189.0 246.45  194.0 270.5  979   100  a 
:::
:::

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

How does this inform one's programming? For example, if you have a
matrix of values stored in column-major order, computing on a column
will be a lot faster than computing on a row, because the column can be
read into the cache from main memory and then accessed in the cache. In
contrast, if the matrix is large and therefore won't fit in the cache,
when you access the values of a row, you'll have to go to main memory
repeatedly to get the values for the row because the values are not
stored contiguously.

There's a nice example of the importance of the cache at [the bottom of
this blog
post](https://wrathematics.github.io/2016/10/28/comparing-symmetric-eigenvalue-performance/).

If you know the size of the cache, you can try to design your code so
that in a given part of your code you access data structures that will
fit in the cache. This sort of thing is generally more relevant if
you're coding in a language like C. But it can matter sometimes in R
too. Here's an example:

::: {.cell hash="efficiency_cache/markdown/cache-aware_61f56d0fb4e03b9d82bdcf6d2016571c"}
``` {.r .cell-code}
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

::: {.cell-output .cell-output-stdout}
                    test elapsed replications
    1  apply(A, 2, mean)  10.742           10
    2 apply(tA, 1, mean)  22.539           10
:::
:::

Now let's compare things when we make the matrix small enough that it
fits in the cache. In this case it fits into the largest (L3) cache but
not the smaller (L2 and L1) caches, and we see that the difference in
speed disappears.

::: {.cell hash="efficiency_cache/markdown/cache-aware2_8595e5a3662a65bf4443c46a625e790e"}
``` {.r .cell-code}
nr <- 800
nc <- 100
## small matrix that should fit in cache
A <- matrix(rnorm(nr * nc), nrow = nr)
## Yep, the size is less than the L3 cache:
object.size(A)
```

::: {.cell-output .cell-output-stdout}
    640216 bytes
:::

``` {.r .cell-code}
memuse::Sys.cachesize()
```

::: {.cell-output .cell-output-stdout}
    L1I:   32.000 KiB 
    L1D:   32.000 KiB 
    L2:   256.000 KiB 
    L3:     8.000 MiB 
:::

``` {.r .cell-code}
tA <- t(A)
benchmark(apply(A, 2, mean),  ## by column
          apply(tA, 1, mean),  ## by row
  replications = 10, columns=c('test', 'elapsed', 'replications'))
```

::: {.cell-output .cell-output-stdout}
                    test elapsed replications
    1  apply(A, 2, mean)   0.011           10
    2 apply(tA, 1, mean)   0.011           10
:::
:::
