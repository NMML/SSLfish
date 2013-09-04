\name{set.ssl.pars}
\alias{set.ssl.pars}
\title{Set sea lion simulation variables}
\usage{
  set.ssl.pars(n.sites, n.yrs, SSLvital = NULL)
}
\arguments{
  \item{n.sites}{Number of sea lion populations}

  \item{n.yrs}{Number of years to run simulation.}

  \item{SSLvital}{Optional list of scaling parameters to
  adjust baseline survival and fecundity from the HFYS
  1970s values.}
}
\description{
  Set sea lion simulation variables
}
\details{
  This function sets the initial simulation variables for
  the sea lion populations. As a baseline the 1970s
  survival and natality parameters of Holmes et al. (2007)
  are used to initialize the female populations. These
  values can be adjusted via the \code{SSLvital} argument
  which allows the user to provide multipliers for these
  parameters on the logit scale. The male survival
  parameters are taken from Calkins and Pitcher (???). The
  age structured mass is calculated from the fitted
  Richards growth models presented in Calkins and Pitcher
  (???) as well.
}
\author{
  Devin Johnson and Paul Conn
}

