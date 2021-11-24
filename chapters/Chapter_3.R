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

# Here F_X(x) = 3 for 0 < x < 1, and F_X^1(u) = u^{1/3}. First we generate all n
# required random uniform numbers as vector u. Then u^{1/3} is a vector of length
# n containing the sample x_1, ... , x_n.
n <- 1000
u <- runif(n)
x <- u ^ (1 / 3)
hist(x, probability = TRUE, main = expression(f(x) == 3 * x ^ 2)) # density hist of sample
y <- seq(0, 1, 0.1)
lines(y, 3 * y ^ 2) # density curve f(x) 
# ... the hist and density plot suggest that the empirical and theoretical 
# distributions approximately agree.

### Example 3.3 -------------------------------------------------------------
# (Exponential distribution) This example applied the inverse transform method to
# generate a random sample from the exponential distribution with mean 1/delta.

# If X ~ Exp(delta), then for x > 0 the cdf of X is F_X(x) = 1 - e^(-delta * x).
# The inverse transformation us F_X^-1(u) = -(1/delta) * log(1 - u). Note that
# U and 1 - U have the same distribution and it is simpler to set 
# x = -(1/delta) * log(u). To generate a random sample of size n with parameter
# lambda:
lambda <- 0.5 # set lambda (mean) constant, i.e. E[X] = 1 / lambda
x <- (-log(runif(n)) / lambda)
hist(x, probability = TRUE, main = expression(
 f(x) == lambda * e ^ {-lambda * x}))
y <- seq(min(x), max(x), 0.1)
lines(y, lambda * exp(-lambda * y)) # density curve f(x) 

## Inverse Transform Method,  Discrete Case --------------------------------