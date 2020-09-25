## @knitr system-time
n <- 10000
m <- 1000
x <- matrix(rnorm(n*m), nrow = n)
system.time({
                mns <- rep(NA, n)
                for(i in 1:n) mns[i] <- mean(x[i , ])
         })
system.time(rowMeans(x))

## @knitr microbenchmark
library(microbenchmark)
df <- data.frame(vals = 1:3, labs = c('a','b','c'))
microbenchmark(
  df[2,1],
  df$vals[2],
  df[2, 'vals']
)

## @knitr microbenchmark2
library(microbenchmark)
n <- 1000
x <- matrix(rnorm(n^2), n)
microbenchmark(
    t(x) %*% x,
    crossprod(x),
    times = 10)

## @knitr benchmark
library(rbenchmark)
# speed of one calculation
benchmark(t(x) %*% x,
          crossprod(x),
          replications = 10,
          columns=c('test', 'elapsed', 'replications'))



## @knitr Rprof-fun
lr_slow <- function(y, x) {
  xtx <- t(x) %*% x
  xty <- t(x) %*% y
  inv <- solve(xtx)   ## explicit matrix inverse is slow and generally a bad idea numerically
  return(inv %*% xty)
}

## @knitr Rprof-run1

## generate random observations and random matrix of predictors
y <- rnorm(5000)
x <- matrix(rnorm(5000*1000), nrow = 5000)

library(proftools)

pd1 <- profileExpr(lr_slow(y, x))
hotPaths(pd1)
hotPaths(pd1, value = 'time')


## @knitr Rprof-run2
lr_medium <- function(y, x) {
  xtx <- crossprod(x)
  xty <- crossprod(x, y)
  inv <- solve(xtx)   ## explicit matrix inverse is slow and generally a bad idea numerically
  return(inv %*% xty)
}                   

pd2 <- profileExpr(lr_medium(y, x))
hotPaths(pd2)
hotPaths(pd2, value = 'time')


## @knitr Rprof-run3
lr_fast <- function(y, x) {
  xtx <- crossprod(x)
  xty <- crossprod(x, y)
  U <- chol(xtx)
  tmp <- backsolve(U, xty, transpose = TRUE)
  return(backsolve(U, tmp))
}

pd3 <- profileExpr(lr_fast(y, x))
hotPaths(pd3)
hotPaths(pd3, value = 'time')


## @knitr Rprof-old

## old approach:
library(fields)

Rprof("makeRegr.prof", interval = 0.005, line.profiling = TRUE)
out <- lr_slow(y, x)
Rprof(NULL)
summaryRprof("makeRegr.prof")




## @knitr preallocate

n <- 10000
z <- rnorm(n)

fun1 <- function(vals) {
   x <- exp(vals[1])
   n <- length(vals) 
   for(i in 2:n) x <- c(x, exp(vals[i]))
   return(x)
}
fun2 <- function(vals) {
   n <- length(vals)
   x <- rep(as.numeric(NA), n)
   for(i in 1:n) x[i] <- exp(vals[i])
   return(x)
}
fun3 <- function(vals) {
  x <- exp(vals)
  return(x)
}
benchmark(fun1(z), fun2(z), fun3(z),
replications = 20, columns=c('test', 'elapsed', 'replications'))

## @knitr init-matrix
nr <- nc <- 2000
benchmark(
   x <- matrix(as.numeric(NA), nr, nc),
   {x <- as.numeric(NA); length(x) <- nr * nc; dim(x) <- c(nr, nc)},
replications = 10, columns=c('test', 'elapsed', 'replications'))

## @knitr init-list
myList <- vector("list", length = n)

## @knitr vectorize
n <- 1e6
x <- rnorm(n)
benchmark(
    x2 <- x^2,
    { x2 <- as.numeric(NA)
      length(x2) <- n
      for(i in 1:n) { x2[i] <- x[i]^2 } },
    replications = 10, columns=c('test', 'elapsed', 'replications'))

## @knitr primitive
`+`
mean.default
chol.default

## @knitr vectorized


address <- c("Four score and seven years ago our fathers brought forth",
             " on this continent, a new nation, conceived in Liberty, ",
             "and dedicated to the proposition that all men are created equal.")
nchar(address)
# use a vector in the 2nd and 3rd arguments, but not the first
startIndices = seq(1, by = 3, length = nchar(address[1])/3)
startIndices
substring(address[1], startIndices, startIndices + 1)

## @knitr vec-tricks
x <- rnorm(1000000)
benchmark(
   truncx <- ifelse(x > 0, x, 0),
   {truncx <- x; truncx[x < 0] <- 0},
   truncx <- x * (x > 0),
   replications = 10, columns=c('test', 'elapsed', 'replications'))

## @knitr apply
n <- 3000; x <- matrix(rnorm(n * n), nr = n)
benchmark(
   out <- apply(x, 1, mean),
   out <- rowMeans(x),
   replications = 10, columns=c('test', 'elapsed', 'replications'))

## @knitr sweep
system.time(out <- sweep(x, 2, STATS = colMeans(x), FUN = "-"))

## @knitr vectorized-sweep
system.time(out2 <- t(t(x) - colMeans(x)))
identical(out, out2)

## @knitr apply-vs-for
n <- 500000; nr <- 10000; nCalcs <- n/nr
mat <- matrix(rnorm(n), nrow = nr)
times <- 1:nr
system.time(
  out1 <- apply(mat, 2, function(vec) {
                         mod = lm(vec ~ times)
                         return(mod$coef[2])
                     })) 
system.time({
                out2 <- rep(NA, nCalcs)
                for(i in 1:nCalcs){
                    out2[i] = lm(mat[ , i] ~ times)$coef[2]
                }
            }) 


## @knitr apply-vs-for-part2

z <- rnorm(10000)
fun2 <- function(vals) {
    x <- as.numeric(NA)
    n <- length(vals)
    length(x) <- n
    for(i in 1:n) x[i] <- exp(vals[i])
    return(x)
}

fun4 <- function(vals) {
    x <- sapply(vals, exp)
    return(x)
}

fun3 <- function(vals) {
    x <- exp(vals)
    return(x)
}

benchmark(fun2(z), fun4(z), fun3(z),
replications = 10, columns=c('test', 'elapsed', 'replications'))


## @knitr matrix-calc
mat <- matrix(rnorm(500*500), 500)
benchmark(apply(mat, 1, sum),
	mat %*% rep(1, ncol(mat)),
	rowSums(mat),
	replications = 10, columns=c('test', 'elapsed', 'replications'))

## @knitr linalg-order
n <- 5000
A <- matrix(rnorm(5000 * 5000), 5000)
B <- matrix(rnorm(5000 * 5000), 5000)
x <- rnorm(5000)
system.time(
  res1 <- A %*% B %*% x
)
system.time(
  res2 <- A %*% (B %*% x)
)

## @knitr diag
n <- 1000
X <- matrix(rnorm(n^2), n) 
diagvals <- rnorm(n)
D = diag(diagvals)
# the following lines are very inefficient
summedMat <- X + D
prodMat1 <- D %*% X
prodMat2 <- X %*% D
# How can we do each of those operations much more quickly?

## @knitr match-lookup
df <- data.frame(
    id = 1:5,
    clusterLabel = c('C', 'B', 'B', 'A', 'C')) 
info <- data.frame(
    grade = c('A', 'B', 'C'),
    numGrade = c(95, 85, 75),
    fail = c(FALSE, FALSE, TRUE) )
grp <- match(df$clusterLabel, info$grade) 
df$numGrade = info$numGrade[grp]
df

## @knitr name-lookup
vals <- rnorm(10)
names(vals) <- letters[1:10]
select <- c("h", "h", "a", "c")
vals[select]

## @knitr index-lookup
n <- 1000
x <- 1:n
xL <- as.list(x)
nms <- as.character(x)
names(x) <- nms
names(xL) <- nms
microbenchmark(
    x[500],  # index lookup in vector
    x["500"], # name lookup in vector
    xL[[500]], # index lookup in list
    xL[["500"]]) # name lookup in list


## @knitr env-lookup
xEnv <- as.environment(xL)  # convert from a named list
xEnv$"500"  
# I need quotes above because numeric; otherwise xEnv$nameOfObject is fine
xEnv[["500"]]
microbenchmark(
  x[500],
  xL[[500]],
  xEnv[["500"]],
  xEnv$"500")


## @knitr cache-aware

nr = 800000
nc = 100
## large matrix that won't fit in cache
A = matrix(rnorm(nr * nc), nrow = nr)
system.time(apply(A, 2, mean))  ## operate by column
A = t(A)
system.time(apply(A, 1, mean))  ## same calculation, but by row

nr = 800
nc = 100
## small matrix that should fit in cache
A = matrix(rnorm(nr * nc), nrow = nr)
tA = t(A)
benchmark(apply(A, 2, mean),  ## by column
          apply(tA, 1, mean),  ## by row
  replications = 10, columns=c('test', 'elapsed', 'replications'))


