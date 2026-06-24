\name{Setup.Grid}
\alias{Setup.Grid}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Grid Setup Function}
\description{
  Produces a list of grid discretisation values in a standardised format.
}
\usage{
Setup.Grid(mmin, m1, mmax, mstep, mout, t1, tmax, tstep, toutmin, toutmax, toutstep, xmin, xmax, xstep, xout, ymin, ymax, ystep, yout)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mmin}{minimum mass value for entire size spectra (in log(g))}
  \item{m1}{maximum mass value for plankton size spectra (in log(g))}
  \item{mmax}{maximum mass value for entire size spectra (in log(g))}
  \item{mstep}{mass step size used in finite differncing calculation. This must be an exact divisor of the total mass range (mmax-mmin)}
  \item{mout}{mass step size that is outputed to results text file. This must be a multiple of the mstep value (but not necessarily a divisor of the total mass range)}
  \item{t1}{time value at which spatial movement starts (in days). (All calculations start at a time value of 0)}
  \item{tmax}{maximum time value to calculate to (in days)}
  \item{tstep}{time step size used in finite differencing calculation. This must be an exact divisor of the total time range (tmax)}
  \item{toutmin}{minimum time value that is outputed to results text file}
  \item{toutmax}{maximum time value that is outputed to results text file}
  \item{toutstep}{time step size that is outputed to results text file}
  \item{xmin}{minimum x-space value (in m)}
  \item{xmax}{maximum x-space value (in m)}
  \item{xstep}{x-space step size used in finite differencing calculation. This must be an exact divisor of the total x-space range (xmax-xmin)}
  \item{xout}{x-space step size that is outputed to results text file}
  \item{ymin}{minimum y-space value (in m)}
  \item{ymax}{maximum y-space value (in m)}
  \item{ystep}{y-space step size used in finite differencing calculation. This must be an exact divisor of the total y-space range (ymax-ymin)}
  \item{yout}{y-space step size that is outputed to results text file}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns an object of class \code{grid.params} containing the input arguments as slots either as specified or containing default values.
}
\references{ ~put references to the literature/web site here ~ }
\author{Matt Castle}
\note{ ~~further notes~~

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{}
