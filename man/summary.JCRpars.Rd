% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functionlist.R
\name{summary.JCRpars}
\alias{summary.JCRpars}
\title{Calculate joint hypothesis testing for the model parameters}
\usage{
\method{summary}{JCRpars}(object, alpha = 0.05, ...)
}
\arguments{
\item{object}{an object of LReg}

\item{alpha}{level of significance (default is 0.05)}

\item{...}{redundant argument}
}
\value{
F statistic, quantile of F distribution at the level of 1-alpha, and significance status
}
\description{
Calculate joint hypothesis testing for the model parameters
}
