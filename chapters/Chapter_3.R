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
set.seed(41513)
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
set.seed(123)
n <- 1000
a <- 3
b <- 2
u <- rgamma(n, shape = a, rate = 1)
v <- rgamma(n, shape = b, rate = 1)
x <- u / (u + v)
q <- qbeta(ppoints(n), a, b)
qqplot(q, x, cex = 0.25, xlab = "Beta(3, 2)", ylab = "Sample")
abline(0, 1)

## Sums & Mixtures ---------------------------------------------------------
# Sums and mixtures of random variables are special types of transformations.

### Convolutions ------------------------------------------------------------
# Let X_1, ..., X_n be i.i.d with X_j \sim X, and let S = X_1 + ... + X_n. The
# distribution function of the sum S is called the n-fold convolution of X and 
# denoted F^{star(n)}_X. It is straightforward to simulate a convolution by directly
# generating X_1, ..., X_n and computing the sum.

# Several distributions are related by convolution. If \nu > 0 is an integer, 
# the \chi^{2} distribution with \nu degrees of freedom is the convolution of
# \nu i.i.d squared standard normal variables. The negative binomial distribution
# NegBin(r, p) is the convolution of r i.i.d Geom(p) r.v's. The convolution of r
# independent Exp(\lambda) r.v's has the Gamma(r, \lambda) distribution.

# In R it is of course easier to generate r.v's from these distributions using
# the rchisq, rgeom, and rbinom functions. Nevertheless, these examples 
# illustrate a general method that can be applied whenever distributions are
# related by convolutions.

#### Example: \chi^{2}
# This example generates a \chi^{2}(\nu) random variable as the convolution of 
# \nu squared normals. If Z_1, ..., Z_{\nu} are i.i.d N(0, 1) r.v's, then V = 
# Z^{2}_1 + ... + Z^{2}_2 has the \chi^{2}(\nu) distribution. To generate a
# random sample of size n from \chi^{2}(\nu), we:

# 1. Fill an n \times \nu matrix with n\nu random N(0, 1) variates
# 2. Square each entry in the matrix(1)
# 3. Compute the row sums of the squared normals. Each row sum is one random 
# observation from the \chi^{2}(\nu) distribution
# 4. Deliver the vector of row sums.

# We develop an example with n = 1000 and \nu = 2 degrees of freedom below.
set.seed(4213)
n <- 1000
nu <- 2
X <- matrix(rnorm(n = n*nu), nrow = n, ncol = nu)^2 # matrix of sq. normals
# Then sum the sq. normals across each row of matrix:
# method 1
y_1 <- rowSums(X)
# method 2:
y_2 <- apply(X, MARGIN = 1, FUN = sum)
identical(y_1, y_2) # test if same... doesn't matter which one is used

# A \chi^{2}(\nu) random variable has mean \nu and variance 2\nu. The samples we
# have generated match this closely with negligible error...
mean(y_1)
mean(y_1^2)

### Mixtures ----------------------------------------------------------------
# A random variable X has a discrete mixture if the distribution of X is a 
# weighted sum F_X(x) = \sum \theta_i F_X(x) for some sequence of random 
# variables X_1, X_2, .... and \theta_i > 0 such that \sum_i \theta_i = 1. The
# constants \theta_i are called the 'mixing' weights or 'mixing' probabilities.
# Although the notation is similar for sums and mixtures, the distributions 
# represented are different.

# A r.v X is a continuous mixture if the distribution of X is 
# F_X(x) = \int_{\inf}^{-\inf} F_{X\mid Y = y}(x) f_Y(y) for a family 
# X \mid Y = y indexed by the real numbers y and weighting function f_Y such 
# that \int_{\inf}^{-\inf} F_{X\mid Y = y}(x) f_Y(y) dy = 1, i.e. it is a valid 
# cumulative density.

# Let's compare the methods for simulation of a convolution and a mixture of
# normal variables. Suppose X_1 \sim N(0, 1) and X_2 \sim N(3, 1) are 
# independent. The notation S = X_1 + X_2 denotes the *convolution* of X_1 and
# X_2. The distribution of S is normal with mean \mu_1 + \mu_2 = 3 and variance
# \sigma_1^{2} + \sigma_2^{2} = 2. Hence, to simulate the convolution, we must:

# 1. Generate x_1 \sim N(0, 1)
# 2. Generate x_2 \sim N(3, 1)
# 3. Deliver s = x_1 + x_2

# We can also define a 50% normal *mixture* X, denoted F_X(x) = 0.5 * F_{X_1}(x) 
# + 0.5 * F_{X_2}(x). Unlike the convolution, the distribution of the mixture X
# is distinctly *non-normal*; it is bimodal. Hence, to simulate the mixture, we
# must:

# 1. Generate an integer k \in \{1, 2\}, where P(1) = P(2) = 0.5
# 2. If k = 1 deliver random x from N(0, 1) else if k = 2 deliver random from
# N(3, 1)

#### Example: Convolutions and mixtures
# Let X_1 \sim Gamma(2, 2) and X_2 \sim Gamma(2, 4) be independent. Let's see 
# the visual difference between the convolution S = X_1 + X_2 and mixture 
# F_X(x) = 0.5 * F_{X_1}(x) + 0.5 * F_{X_2}(x):
n <- 1000
x_1 <- rgamma(n, 2, 2)
x_2 <- rgamma(n, 2, 4)
s <- x_1 + x_2 # convolution
u <- runif(n)
k <- as.integer(u > 0.5) # vector of 0's and 1's (binary logic)
x <- k * x_1 + (1 - k) * x_2 # the mixture
# Plot:
par(mfcol = c(1, 2))
hist(s, probability = TRUE, xlim = c(0, 5), ylim = c(0, 1))
hist(x, probability = TRUE, xlim = c(0, 5), ylim = c(0, 1))
par(mfcol = c(1, 1)) # restore previous setting

#### Example: Mixture of several gamma distributions
# Here there are several components to the mixture and the mixing weights are 
# not uniform. The mixture is

# F_X = \sum_{i = 1}^5 \theta_j F_{X_j}

# where X_j \sim Gamma(r = 3, \lambda_{j} = \frac{1}{j}) are independent and the
# mixing probabilities are \theta_j = \frac{j}{15}, j = 1, ..., 5.

# To simulate one random variate from the mixture F_X:

# 1. Generate an integer k \sim \{1, 2, 3, 4, 5\}, where P(k) = \theta_k, 
# k = 1, ..., 5
# 2. Deliver a random Gamma(r, \lambda_{k}) variate

# Hence, to generate sample size n, steps 1 and 2 are repeated n times. Since
# for loops can be costly in R, we can translate the above into a vectorised
# approach. So then, we

# 1. Generate a random sample k_1, ..., k_n of integers in a vector k, where
# P(k) = \theta_k, k = 1, ..., 5. Then k[i] indicates which of the the five 
# gamma distributions will be sampled to get the ith element of the sample.
# 2. Set rate equal to the length n vector \lambda = (\lambda_k)
# 3. Generate a gamme sample size n, with shape r and rate vector rate

# Thus, and efficient way to implement this in R is shown below:
n <- 5000
k <- sample(1:5, size = n, replace = TRUE, prob = (1:5) / 15)
rate <- 1 / k
x <- rgamma(n = n, shape = 3, rate = rate)
# Plot the density of the mixture with densities of the compoenents...
plot(x = density(x),
     xlim = c(0, 40), ylim = c(0, 0.3), 
     lwd = 3, xlab = "x", main = "Mixture (thick) of Independent Gammas (thin)")
for (i in 1:5) {
 lines(density(rgamma(n = n, shape = 3, rate = 1/i)))
}
# Thick line is the final mixture; thin lines are the density of each X_j.

# The rest of the examples are repetitions of the same theme. Onto the next 
# section!

# Multivariate Distributions ----------------------------------------------
# A d-dimensional multivariate normal requires covariance matrix, since 
# X ~ N_d(\mu, \Sigma),where \Sigma has entries \sigma_{ij} = Cov(X_i, X_j).
# One of the ways - there are several - to obtain a covariance matrix is 
# through Choleski factorization. The Choleski factorization of a real symmetric
# positive-definite matrix is X = Q^TQ, where Q is an upper triangular matrix.
# In R, the basic syntax to implement Choleski factorization is chol(X).
X <- matrix(data = c(1, 0.9, 0.9, 1), nrow = 2)
# Check is if matyric is symmetric and positive definite:
X - t(X) ## This should be 0 if A is symmetric
eigen(X)$values ## These should all be non-negative.
# Note that if any are 0 then use chol(A,pivot=TRUE)

## Calculate Choleski factorization:
R <- chol(X)
# Thus, if correct, XX = R^T%*%R:
t(R)%*%R

