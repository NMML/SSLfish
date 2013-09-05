\name{set.A.array}
\alias{set.A.array}
\title{Set sea lion projection matrix}
\usage{
  set.A.array(Svars, yr)
}
\arguments{
  \item{Svars}{Named list created with the
  \code{\link{set.ssl.pars}} function.}

  \item{yr}{Current year.}
}
\description{
  Takes survival and natality information from the sea lion
  simulation variables object and create a Leslie
  projection matrix for each site.
}
\details{
  Takes the survival and natality information for each site
  in year \code{yr} and places the information into a
  Leslie matrix for population projection to year
  \code{yr}+1.
}
\author{
  Devin Johnson
}

