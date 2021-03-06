% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sarim.R
\name{sarim}
\alias{sarim}
\title{Generate the structure/penalty matrix for sarim formulae}
\usage{
sarim(formula, data = list(), intercept = "FALSE", nIter = 1000L,
  burnin = 100L, family = "gaussian", link = "identity", sigma = 0.1,
  sigma_a = 1e-04, sigma_b = 1e-04, Ntrials = 1L, m = 250L,
  thr = 1e-04)
}
\arguments{
\item{formula}{usual formula, y ~ x1 + sx(x2), be aware an intercept can only be
included by setting intercept = "TRUE".}

\item{data}{a data.frame.}

\item{intercept}{decide if an intercept should be included or not.}

\item{nIter}{number of MCMC-Iteration.}

\item{burnin}{number of Burn-In iterations (currently NOT available).}

\item{family}{family can currently only be "gaussian", "binomial" or "poisson".}

\item{link}{link-function, currently only be "identity", "logit" or "log" implemented,
if family = "binomial", please choose "logit", if family = "poisson", please choose "log".}

\item{sigma}{start value for the dispersion parameter, it is assumed that sigma ~ IG(sigma_a, sigma_b),
i.e. it is assumed that sigma is Inverse Gamma distributed with parameters sigma_a and sigma_b.}

\item{sigma_a}{prior-value, see sigma.}

\item{sigma_b}{prior-value, see sigma.}

\item{Ntrials}{Number of trails, especially for binomial distribution.}

\item{m}{number of maximal Lanczos-iteration for start iteration, after first iteration, it use the 
number of the first + 10, in the third iteration it use the number of the second + 10; 
for saving memory.}

\item{thr}{threshold, when the Lanczos-algorithm or conjugate gradient-algorithm should break.}
}
\description{
Function used in definition of smooth terms within sarim model formulae.
}
