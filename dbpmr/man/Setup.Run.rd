\name{Setup.Run}
\alias{Setup.Run}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run Setup Function}
\description{
  Produces a list of the main run parameters for the model under investigation
}
\usage{
Setup.Run(no.species, spatial.dim, coupled.flag, initial.flag, filename)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{no.species}{an integer value for the number of pelagic species in the model}
  \item{spatial.dim}{a value determining the spatial dimension of the model. Either 0, 1 or 2}
  \item{coupled.flag}{a logical value specifying whether the species are to interact}
  \item{initial.flag}{a logical value specifying whether the initial distributions are to be inputed by the user}
  \item{filename}{a root filename for the output text files}
}
\details{

}
\value{
  Returns an object of class \code{run.params} containing the input arguments as slots.
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
