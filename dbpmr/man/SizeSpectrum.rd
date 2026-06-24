\name{SizeSpectrum}
\alias{SizeSpectrum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Dynamic SizeSpectrum Calculation Function}
\description{
Outputs 3 text files per species into the working directory. One containing the actual number densities, one containg the groth values and one containing the mortality values.
}
\usage{
SizeSpectrum(run.in, grid.in, species.in, sp.initial)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{run.in}{An object of class \code{run.params}}
  \item{grid.in}{An object of class \code{grid.params}}
  \item{species.in}{A list of objects of class \code{species.params}, one for each species to be included}
  \item{sp.initial}{An optional list of filenames containing the initial distributions of the species to be included}
}
\details{
This function makes a call to C code to perform the dynamic size spectrum calculation, using the values specified in the input arguments.
}
\value{
  Returns a list of txt files containing the output produced.
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
