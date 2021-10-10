# mtool
An R package for Structural Penalization in Cox Model with Time-dependent Covariates

How to load the package:
Under the mtool directory, type `R CMD build mtool` to load the package.

## Examples:
The example data generating mechanism is seen in GenerateData.R

How to perform cross-validation and the test of code is seen in TEST.R

An example of obtained regularization path is seen in S1_Regularization Path.pdf

## Functions:

1) `lik (betas, covariates, Id, Event, Fup, Start, Stop)`

*Return the minus log likelihood of Cox model divided by the sample size*

2) `der_lik (betas, covariates, Id, Event, Fup, Start, Stop)`

*Return the minus derivative log likelihood of Cox model divided by the sample size*

3) `SurvGraphSelect(covariates, Id, Event, Fup, Start, Stop, grp, grpV, etaG, regul, betas, t, alpha, epsilon, lam1, lam2 = 0.0, lam3 = 0.0, num_threads = -1L, intercept = FALSE, resetflow = FALSE, verbose = FALSE, pos = FALSE, clever = TRUE, eval = TRUE, size_group = 1L, transpose = FALSE)`

*Return the coefficients estimate for a specific penalization hyper-parameter lambda*

## Arguments:

betas: a vector of initial values of coefficients

covariates: a matrix, which has p columns of covariates

Id: ID for each of the subject

Event: a vector of binary values, 0 indicates in the period of time, no event occured; 1 indicates in the period of time, event occurs.

Fup: a vector of integers, the same value for each of subject, indicates the follow-up time

Start: a vector indicates the start time of the period of time

Stop: a vector indicates the stop time of a period of time

g: the number of groups

grp: a g by g sparse matrix, if the i-th group is nested in the j-th goup, then the i-th row, j-th column of the matrix is 1, otherwise 0

grpV: a p by g sparse matrix, if the i-th covariate is in the j-th group, then the i-th row, j-th column of the matrix is 1, otherwise 0. However, if the k-th group is nested in the l-th group, then the covariates overlapped in the two groups should not be specified in the l-th group

etaG: the weights of each group

t: step size

alpha: shrinkage rate

epsilon: convergence control

lam1: penalization hyper-parameter

numThreads: (optional, number of threads for exploiting multi-core / multi-cpus. By default, it takes the value -1, which automatically selects all the available CPUs/cores)

# Authors
Guanbo wang (guanbo.wang@mail.mcgill.ca), Yi Lian, Yi Yang

McGill University

