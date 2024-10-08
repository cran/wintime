% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/comp.R
\name{COMP}
\alias{COMP}
\title{Run composite analysis}
\usage{
COMP(n, Time, Delta, cov, trt)
}
\arguments{
\item{n}{The total number of trial participants.}

\item{Time}{A m x n matrix of event times (days), where m is the number of events in the hierarchy. Rows should represent events and columns
should represent participants. Event rows should be in increasing order of clinical severity.}

\item{Delta}{A m x n matrix of event indicators, where m is the number of events in the hierarchy. Rows should represent events and
columns should represent participants. Event rows should be in increasing order of clinical severity.}

\item{cov}{A n x p matrix of covariate values, where p is the number of covariates. Rows should represent participants and columns
should represent covariate values.}

\item{trt}{A vector of length n containing treatment arm indicators (1 for treatment, 0 for control).}
}
\value{
The z-statistic of the treatment effect from the Cox Model fit.
}
\description{
This function fits a Cox Model to time-to-event data and calculates the z statistic. In the wintime package, this function is used for the
EWTR-composite max test (MAX) method.
}
