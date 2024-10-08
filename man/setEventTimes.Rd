% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{setEventTimes}
\alias{setEventTimes}
\title{Created a sorted vector of event times}
\usage{
setEventTimes(m, delta0, delta1, time0, time1, follow)
}
\arguments{
\item{m}{The number of events in the hierarchy.}

\item{delta0}{A vector of event indicators for the control person.}

\item{delta1}{A vector of event indicators for the treatment person.}

\item{time0}{A vector of event times (days) for the control person.}

\item{time1}{A vector of event times (days) for the treatment person.}

\item{follow}{The maximum follow up time (days).}
}
\value{
A sorted vector of event times (days) for a given pair.
}
\description{
This function creates a sorted vector of event times for a pair. This function is used in all pairwise functions.
}
