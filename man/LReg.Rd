% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LReg.R
\name{LReg}
\alias{LReg}
\title{A main function that calculates useful values used to fit the linear model}
\usage{
LReg(myformula, mydata)
}
\arguments{
\item{myformula}{a formula of the model that includes response and covariates}

\item{mydata}{data where we fit the linear model}
}
\value{
a list of useful values
}
\description{
A main function that calculates useful values used to fit the linear model
}
\examples{
LReg(Species ~ Area +Elevation+ Nearest+ Scruz+ Adjacent, mydata1)
}
