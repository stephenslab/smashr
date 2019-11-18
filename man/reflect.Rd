% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ashsmooth.R
\name{reflect}
\alias{reflect}
\title{Reflect and extend a vector.}
\usage{
reflect(x)
}
\arguments{
\item{x}{An n-vector.}
}
\value{
A list with two list elements: \code{"x"} containing the
  reflected signal; and \code{"idx"} containing the indices of the
  original signal.
}
\description{
Reflects a vector if it has length a power of 2;
  otherwise extends the vector to have length a power of 2 and then
  reflects it.
}