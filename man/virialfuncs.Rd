\name{virialfuncs}
\alias{virialfuncs}
\alias{MvirToSigma}
\alias{SigmaToMvir}
\alias{MvirToRvir}
\alias{RvirToMvir}
\alias{SigmaToRvir}
\alias{RvirToSigma}
\alias{Mpart}
\alias{Npart}

\title{
Virial conversion functions
}
\description{
All 6 Virial paramter conversion functions. Each can map precisely to the other as a one paramter function.
}
\usage{
MvirToSigma(Mvir = 1e+12, Munit = 1, Lunit = 1000, Vunit = 1, DeltaVir = 200, h = 1)
SigmaToMvir(Sigma = 230, Munit = 1, Lunit = 1000, Vunit = 1, DeltaVir = 200, h = 1)
MvirToRvir(Mvir=1e12, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, h=1)
RvirToMvir(Rvir=162.635, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, h=1)
SigmaToRvir(Sigma=230, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, h=1)
RvirToSigma(Rvir=162.635, Munit=1, Lunit=1e3, Vunit=1, DeltaVir=200, h=1)
Mpart(VirialParam=1e12,Npart=1000,type='mass',Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1,int=FALSE)
Npart(VirialParam=1e12,Mpart=1e9,type='mass',Munit=1,Lunit=1e3,Vunit=1,DeltaVir=200,h=1,int=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Mvir}{
Mass within virial radius in units on 'Munit'
}
  \item{Sigma}{
Velocity dispersion within virial radius in units on 'Vunit'
}
  \item{Rvir}{
Size of virial radius in units on 'Lunit'
}
  \item{VirialParam}{
The desired virial parameter, either in Munits of mass (type='mass'), Vunits of velocity dispersion (type='sig') or Lunits of virial radius (type='rad')- default is mass
}
  \item{Mpart}{
The mass of each particle in units Munit
}
  \item{Npart}{
The number of particles within the virial radius
}
  \item{type}{
Tells the code which tpye of VirialParam has been entered, either mass (mass), velosity dispersion (sig) or Virial radius (rad)- the user should be careful that the units are as expected from Munit / Vunit or Lunit respectively
}
  \item{Munit}{
Base mass unit in multiples of Msun
}
  \item{Lunit}{
Base length unit in multiples of parsecs
}
  \item{Vunit}{
Base velocity unit in multiples of km/s
}
  \item{DeltaVir}{
Relative over-density virial radius
}
  \item{h}{
Hubble constant in units of H0/100
}
  \item{int}{
Should arguments by returned as integers, for Npart in particlar this is probably important
}
}

\details{
These functions allow for various analytic conversions between the 3 major properties related to virial radius: the mass, velocity dispresion and size. The default properties calculate properties for 1e12 Msun halos and assume masses in Msun, velocities in km/s and distances in Kpc.
}

\value{
MvirToSigma outputs velocity (in units of Vunit) dispersion given mass
SigmaToMvir outputs mass (in units of Munit) given velocity dispersion
MvirToRvir outputs radius (in units of Lunit) given mass
RvirToMvir outputs mass (in units of Munit) given radius
SigmaToRvir outputs radius (in units of Lunit) given velocity dispersion
RvirToSigma outputs velocity (in units of Vunit) dispersion given radius
Mpart outputs the required mass of each particle
Npart outputs the number of particles within the virial radius
}

\author{
Aaron Robotham, Chris Power
}

\seealso{
\code{\link{snapwrite}},\code{\link{genhalo}},\code{\link{addhead}},\code{\link{genparam}},\code{\link{cosdist}}
}

\examples{
temp=MvirToSigma(1e13)
print(temp) # Default
temp2=MvirToSigma(1e13, Vunit=1e-3)
print(temp2) # Velocity in m/s
temp3=SigmaToMvir(temp2,Vunit=1e-3)
print(temp3)
temp4=SigmaToRvir(temp2,Vunit=1e-3)
print(temp4)
#How many particles to we need for a 200 km/s halo where each particle is 1e8 solar masses?
temp5=Npart(200,1e8,type='sig')
print(temp5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{gadget}

