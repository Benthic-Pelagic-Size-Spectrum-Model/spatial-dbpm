\name{Setup.Initial}
\alias{Setup.Initial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Initial Values Setup Function}
\description{
  Produces a text file containing the initial values for a given species
}
\usage{
Setup.Initial(grid.in,mfunc,xfunc,yfunc)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{grid.in}{minimum mass value for entire size spectra (in log(g))}
  \item{mfunc}{maximum mass value for plankton size spectra (in log(g))}
  \item{xfunc}{maximum mass value for entire size spectra (in log(g))}
  \item{yfunc}{mass step size used in finite differncing calculation. This must be an exact divisor of the total mass range (mmax-mmin)}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns the filename of the text file produced.
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
