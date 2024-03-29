% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CCCNetwork.R
\name{CCCNetwork}
\alias{CCCNetwork}
\title{Cell-cell Communiction Network}
\usage{
CCCNetwork(
  restab,
  ligand = NULL,
  receptor = NULL,
  pathway = NULL,
  cutoff = 0.2,
  weight = 10,
  ...
)
}
\arguments{
\item{restab}{A \code{data.frame} returned by \code{FRETCOM()}.}

\item{ligand}{A specific ligand name in \code{restab}, should be paired with \code{receptor}.}

\item{receptor}{A specific receptor name in \code{restab}, should be paired with \code{ligand}.}

\item{pathway}{A specific pathway name in \code{restab} (if any). Once \code{pathway} is
defined, \code{ligand} and \code{receptor} will not be considered.}

\item{cutoff}{A numeric value within the range of \verb{(0, 1)} indicates the cutoff
for the adjusted p-values. If the adjusted p-value is less than
the specified cutoff, we treat it as significant, and it will be
displayed in the heatmap. The default is \code{0.2}.}

\item{weight}{A numeric value indicates the weight that determines the edge width.
Edge width is calculated as the product of the estimate and the weight.
The default is \code{10}.}

\item{...}{Other arguments passed on to \code{igraph::plot()}.}
}
\value{
A network produced by \code{igraph}.
}
\description{
Cell-cell Communiction Network
}
