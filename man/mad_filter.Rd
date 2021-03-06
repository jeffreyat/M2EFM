% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m2eQTL.R
\name{mad_filter}
\alias{mad_filter}
\title{Takes an ExpressionData object and filters it by the given MAD threshold.}
\usage{
mad_filter(exp, thres = 0.05)
}
\arguments{
\item{exp}{an ExpressionData object}

\item{thres}{MAD threshold above which to keep genes}
}
\value{
a data.frame with gene locations
}
\description{
Takes an ExpressionData object and filters it by the given MAD threshold.
}

