\name{Average.Time}
\alias{Average.Time}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Time Averaging Function}
\description{
A function to average the output from several timesteps given a time range
}
\usage{
Average.Time(species,time.lim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{species}{an object of class \code{species.results}}
  \item{time.lim}{a vector of length 2 specifying the start and end times of the range to be averaged over}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns an object of class \code{timestep.data}which has six slots;
  \item{data}{a data frame containg the a timestep value of the averaged results}
  \item{spatial.dim}{an integer specifying the spatial dimension of the results}
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
