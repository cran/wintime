% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perm.R
\name{perm}
\alias{perm}
\title{Resample using permutations}
\usage{
perm(
  type,
  time_restriction,
  model,
  n,
  m,
  Time,
  Delta,
  trt,
  cov,
  z_ewtr,
  z_comp,
  resample_num,
  seed,
  nimp
)
}
\arguments{
\item{type}{A string value indicating the wintime package method that will run with resampling.}

\item{time_restriction}{The time cutoff value (days).}

\item{model}{A string value indicating the model used on observed data ('markov' or 'km').}

\item{n}{The total number of trial participants.}

\item{m}{The number of events in the hierarchy.}

\item{Time}{A m x n matrix of event times (days). Rows should represent events and columns should represent participants. Event rows should be
in increasing order of clinical severity.}

\item{Delta}{A m x n matrix of event indicators. Rows should represent events and columns should represent participants. Event rows should
be in increasing order of clinical severity.}

\item{trt}{A numeric vector of treatment arm indicators (1 for treatment, 0 for control).}

\item{cov}{A n x p matrix of covariate values, where p is the number of covariates. Rows should represent participants and columns should
represent covariate values.}

\item{z_ewtr}{The Z-statistic of EWTR.}

\item{z_comp}{The Z-statistic of the composite event approach.}

\item{resample_num}{The number of desired permutations.}

\item{seed}{The seed used for random number generation.}

\item{nimp}{The number of random imputations for Redistribution-to-the-right.}
}
\value{
A list of a vector of length resample_num containing the treatment effect estimates (for type='max' these are z-statistics) for each permutation,
a m x resample_num matrix of the components of the treatment effect..
}
\description{
This function reruns the desired wintime package method on a given number of permutations. This resampling method is recommended
for the Expected win time (EWT) and Restricted mean survival in favor of treatment (RMT) methods.
}
