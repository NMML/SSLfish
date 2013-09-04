\name{update.fleet.pars}
\alias{update.fleet.pars}
\title{Update fishing fleet simulation variables}
\usage{
  update.fleet.pars(n.sites, iyr, Fvars, FLvars, sim.opts,
    burnin)
}
\arguments{
  \item{n.sites}{Number of sites in the simulation}

  \item{iyr}{Current year in the simulation}

  \item{Fvars}{Fish population variables list created by
  call to \code{\link{set.fish.pars}}}

  \item{FLvars}{Fishing fleets variables list created by
  call to \code{\link{set.fleet.pars}}}

  \item{sim.opts}{???}

  \item{burnin}{Number of years in the simulation before
  'fishing' begins.}
}
\description{
  Update fishing fleet simulation variables
}
\author{
  Paul Conn
}

