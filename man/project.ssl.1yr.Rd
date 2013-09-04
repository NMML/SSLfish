\name{project.ssl.1yr}
\alias{project.ssl.1yr}
\title{Project sea lion populations}
\usage{
  project.ssl.1yr(Svars, yr)
}
\arguments{
  \item{Svars}{Named list of sea lion variables craeted
  from \code{\link{set.ssl.pars}}.}

  \item{yr}{Current year.}
}
\description{
  This function updates the sea lion population by one
  year.
}
\details{
  This function takes survival and natality information
  from the \code{Svars} object to update the sea lion
  population from year 'yr' to 'yr+1'. The simulated survey
  sample of the population is accomplished at the same
  time.
}
\author{
  Devin Johnson
}

