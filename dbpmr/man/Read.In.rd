\name{Read.In}
\alias{Read.In}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Data Reading Function}
\description{
This reads in the output from the text files produced by the function \code{SizeSpectrum}
}
\usage{
Read.In(filename)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{A full filename and path to a text file produced by the function \code{SizeSpectrum}}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns an object of class \code{species.results} which has 8 slots containing all of the information used to produce the output file.
  These are
    \item{results}{a data frame containing the actual species results in a standardised format}
    \item{run}{a data frame containing the values and strings used in the model run}
    \itme{grid}{an object of class \code{grid.params} containing all the values used in the grid discretisation}
    \item{species}{an object of class \code{species.params} containing all the parameters used to generate the species results}
    \item{trange}{a numeric vector containing all of the t values outputed}
    \item{mrange}{a numeric vector containing all of the m values outputed}
    \item{xrange}{a numeric vector containing all of the x values outputed}
    \item{yrange}{a numeric vector containing all of the y values outputed}
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
