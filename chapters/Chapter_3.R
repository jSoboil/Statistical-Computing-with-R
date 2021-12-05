# Methods for Generating Random Variables ---------------------------------

## The Inverse Transform Method --------------------------------------------
# The inverse transform method of generating RVs is based on the following 
# well-known result:

## Theorem 3.1 (Probability Integral Transformation) If X is a continuous random 
# variable with cdf F_X(x), then U = F_X(X) ~ Uniform(0, 1).

# The inverse transform method of generating RVs applied the probability 
# integral transformation, where the inverse transformation is defined as

## F_x^-1(u) = inf{x: F_X(x) = u}, 0 < u < 1

# Hence, if U ~ Uniform(0, 1), then for all xER
## P(F_X^-1(U) <= x) = P(inf{t: F_X(t) = U} <= x)
##                   = P(U <= F_X(x))
##                   = F_U(F_X(x)) = F_X(x),

# and therefore F_X^-1(u) has the same distribution as X. Thus, to generate a 
# random observation X, first generate a Uniform(0, 1) variate y and deliver the 
# inverse value F_X^-1(u). The method is easy to apply *provided that the inverse
# density function is easy to compute

### Example 3.2 -------------------------------------------------------------
# This example use the inverse method to simulate a random sample from the 
# distribution with density f(x) = 3x^2, 0 < x < 1.

# Here F_X(x) = x^3 for 0 < x < 1, and F_X^-1(u) = u^{1/3}. First we generate all n
# required random uniform numbers as vector u. Then u^{1/3} is a vector of length
# n containing the sample x_1, ... , x_n.
n <- 1000
u <- runif(n)
x <- u ^ (1 / 3)
hist(x, probability = TRUE, main = expression(f(x) == 3 * x ^ 2)) # density hist of sample
y <- seq(0, 1, 0.01)
lines(y, 3 * y ^ 2) # density curve f(x)
# ... the hist and density plot suggest that the empirical and theoretical 
# distributions approximately agree.

### Example 3.3 -------------------------------------------------------------
# (Exponential distribution) This example applied the inverse transform method to
# generate a random sample from the exponential distribution with mean 1/lambda.

# If X ~ Exp(lambda), then for x > 0 the cdf of X is F_X(x) = 1 - e^(-lambda * x).
# The inverse transformation us F_X^-1(u) = -(1/lambda) * log(1 - u). Note that
# U and 1 - U have the same distribution and it is simpler to set 
# x = -(1/lambda) * log(u). To generate a random sample of size n with parameter
# lambda:
lambda <- 0.5 # set lambda (mean) constant, i.e. E[X] = 1 / lambda
x <- (-log(runif(n)) / lambda)
hist(x, probability = TRUE, main = expression(
 f(x) == lambda * e ^ {-lambda * x}))
y <- seq(min(x), max(x), 0.1)
lines(y, lambda * exp(-lambda * y)) # density curve f(x) 

## Inverse Transform Method,  Discrete Case --------------------------------
# The inverse transform method can also be applied to the discrete case. If X is 
# a discrete random variable and ... < x_[i-1] < x_[i] < x_[i + 1] < ... are the
# points of disconuity of F_[X](u), then the inverse transformation is 
# F_[X]^-1(u) = x_[i], where F_[X](x_[i - 1]) < u <= F_[X](x_[i]). So, for each
# random variate required:

## 1. Generate a random u from Uniform(0, 1)
## 2. Deliver x_[i] where F(x_[i - 1]) < u <= F(x_[i])

### Example 3.4 -------------------------------------------------------------
# (Two point distribution) This example applied the inverse transform method to
# generate a random sample of Bernoulli(p = 0.4) variates.

# In this example, F_X(0) = f_X(0) = 1 - p and F_X(1) = 1. Thus, F_X^-1(u) = 1
# if u > 0.6 and F_X^-1(u) = 0 if u <= 0.6. The generator should therefore 
# deliver the numerical value of the logical expression u > 0.6:
set.seed(300)
n <- 1000
p <- 0.4
u <- runif(n)
x <- as.integer(u > 0.6)
mean(x)

### Example 3.6 -------------------------------------------------------------
rlogarithmic <- function(n, theta) {
 # returns a random logarithmic(theta) sample size n:
 u <- runif(n)
 # set the initial length of the cdf vector:
 N <- ceiling(-16 / log10(theta))
 k <- 1:N
 a <- -1 / log(1 - theta)
 fk <- exp(log(a) + k * log(theta) - log(k))
 Fk <- cumsum(fk)
 x <- integer(n)
 for (i in 1:n) {
  x[i] <- as.integer(sum(u[i] > Fk)) #F^{-1}(u) = 1
  while (x[i] == N) {
   # if x == N we need to extend the cdf
   logf <- log(a) + (N + 1) * log(theta) - log(N + 1)
   fk <- c(fk, exp(logf))
   N <- N + 1
   x[i] <- as.integer(sum(u[i] > Fk))
  }
 }
 x + 1
}

# Using the above function, we can then generate random samples from the 
# Logarithmic(0.5) distribution:
n <- 1000
theta <- 0.5
x <- rlogarithmic(n = n, theta = theta)
k <- sort(unique(x))
p <- -1 / log(1 - theta) * theta ^ k / k
se <- sqrt(p * (1 - p) / n)

## The Acceptance-Rejection Method -----------------------------------------
# Suppose that X and Y are random variables with with density or pmf f and g,
# respectively, and there exists a constant c such that

# \frac{f(t)}{g(t)} <= c

# for all t such that f(t) > 0. Then the acceptance-rejection method (or 
# rejection method) can be applied to generate the random variable X.

##  The Acceptance-Rejection Method ----------------------------------------
# 1. Find a random variable Y with density g satisfying f(t) / g(t) <= c, for
# all t such that f(t) > 0. Provide a method to generate random Y.

# 2. For each random variate required:
## a) Generate a random y from the distribution with density g
## b) Generate a random u from the Uniform(0, 1) distribution
## c) If u < f(y) / cg(y), accept y and deliver x = y; otherwise reject y.

# Note that in step 2c,

# P(accept | Y) = P(U < \frac{f(Y)}{cg(Y)} | Y) = \frac{f(Y)}{cg(Y)}

# The last equality is simply evaluating the cdf of U. The total probability of
# acceptance for any iteration is therefore

# \Sigma_{y} P(accept | y)P(Y = y) = |Sigma_{y} \frac{f(y)}{cg(y) g(y) = 1/c}.

# and the number of iterations until acceptance has the geometric distribution 
# with mean c. Hence, on average each sample value of X requires c iterations.
# For efficiency, Y should be easy to simulate and c small.

# To see that the accepted sample has the same distribution as X, apply Bayes'
# theorem.

### Example 3.7 -------------------------------------------------------------
# This example illustrates the method for the Beta distribution. On average, 
# how many random numbers must be simulated to generate 1000 variables from the
# Beta(/alpha = 2, \beta = 2) distribution by this method? It depends on the
# upper bound c of \frac{f(x)}{g(x)}, which depends on the choice of the function 
# g(x).

# The Beta(2, 2) density is f(x) = 6x(1 - x), where 0 < x < 1. Let g(x) be the
# Uniform(0, 1) density. Then \frac{f(x)}{g(x)} <= 6 for all 0 < x < 1, and so
# c = 6. A random x dfrom g(x) is accepted if

# \frac{f(x)}{cg(x)} = \frac{6x(1 - x)}{6(1)} = x(1 - x) > u

# On average, cn = 6000 iterations (12000 random numbers) will be required for a
# sample size 1000. In the following simulation below, the counter j for 
# iterations is not necessary, but included to record how many iterations were
# actually needed to generate the 1000 beta variates:
n <- 1000
k <- 0
j <- 0 # counter for accepted samples
y <- numeric(n)

while (k < n) {
 u <- runif(1)
 j <- j + 1
 x <- runif(1) # random variate from g
 
 if (x * (1 - x) > u) {
  # accept sample x:
  k <- k + 1
  y[k] <- x
 }
}
j

# Now we can compare the empirical and theoretical percentiles:
p <- seq(0.1, 0.9, 0.1)
Q_hat <- quantile(y, p)
Q <- qbeta(p, 2, 2)
se <- sqrt(p * (1 - p)) / (n * dbeta(Q, 2, 2) ^ 2)

# Below, the sample percentiles (first line) approximately match the Beta(2, 2)
# percentiles computed by qbeta (second line):
round(rbind(Q_hat, Q, se), 3)

## Transformation Methods --------------------------------------------------






















