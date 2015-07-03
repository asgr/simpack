\name{simpack-package}
\alias{simpack-package}
\alias{simpack}
\docType{package}
\title{
 Functions for calculating orbits and interacting with the Gadget N-body code
~~ package title ~~
}
\description{
General purpose orbital functions. Functions for operating with, and generating, Gadget binary snapshots. This includes read/write functions and IC generators.
}
\details{
\tabular{ll}{
Package: \tab simpack\cr
Type: \tab Package\cr
Version: \tab 0.13\cr
Date: \tab 2013-07-04\cr
License: \tab GPL-2 \cr
}

}
\author{
Aaron Robotham, Chris Power

Maintainer: Aaron Robotham <aaron.robotham@uwa.edu.au>

}

\examples{
\dontrun{
temp=snapread('snapshot_XXX')
temp$part[,'x']=temp$part[,'x']+10
snapwrite(temp$part,temp$head,'snapshot_XXX_mod')
}
}
