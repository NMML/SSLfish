\name{run.1.sim}
\alias{run.1.sim}
\title{Execute one simulation run}
\usage{
  run.1.sim(par, sim.opts, n.sites, n.yrs, burnin)
}
\arguments{
  \item{par}{Vector of simulation variables subject to
  calibration}

  \item{sim.opts}{???}

  \item{n.sites}{Number of sites within the simulation}

  \item{n.yrs}{Number of total years in the simulation}

  \item{burnin}{???}
}
\description{
  This function takes the arguments of \code{par} and runs
  the population simulation 1 time for the purposes of
  calibration.
}
\details{
  This function runs a single simulation run so that the
  free variables in the par vector can be calibrated to
  provide realistic population growth and interaction
  between sea lions, prey and fisheries components.
}
\author{
  Paul Conn
}

