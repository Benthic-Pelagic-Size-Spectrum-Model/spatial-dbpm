\name{Setup.Species}
\alias{Setup.Species}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Species Setup Function}
\description{
  Produces a list of species specific constants in a standardised format.
}
\usage{
Setup.Species(alpha, beta, lambda, A, K, mu_0, mu_s, epsilon, u_0, sig, q_0, prey, pred, comp, gamma_prey, gamma_pred, gamma_comp, speciesname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{The volume search allometric exponent}
  \item{beta}{The intrinsic mortality allometric exponent}
  \item{lambda}{The slope of the default plankton size spectrum}
  \item{A}{The volume search coefficient}
  \item{K}{The energy conversion efficiency}
  \item{mu_0}{The intrinsic mortality coefficient}
  \item{mu_s}{The senescence mortality coefficient}
  \item{epsilon}{The senescence mortality 'sharpness' coefficient}
  \item{u_0}{The intercept of the default plankton size spectrum}
  \item{sig}{The feeding mass ratio kernel standard deviation}
  \item{q_0}{The feeding mass ratio kernel mean}
  \item{prey}{The prey seeking spatial movement coefficient}
  \item{pred}{The predator avoiding spatial movement coefficient}
  \item{comp}{The competition avoiding spatial movement coefficient}
  \item{gamma_prey}{The prey seeking spatial movement allometric exponent}
  \item{gamma_pred}{The predator avoiding spatial movement allometric exponent}
  \item{gamma_comp}{The competition avoiding spatial movement allometric exponent}
  \item{speciesname}{The name of the species. This must be unique if more than one species is under investigation as it used to construct the output results text files.}
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  Returns an object of class \code{species.params} list containing the arguments either as specified or containing default values.
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
