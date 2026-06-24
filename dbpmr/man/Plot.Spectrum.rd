\name{Plot.Spectrum}
\alias{Plot.Spectrum}
\alias{Points.Spectrum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{A Multi Dimensional Size Spectrum Plotting Function}
\description{
This is a generic function that can plot various 3D surface plots, or 2D transects or Size Spectrum plots.
}
\usage{
Plot.Spectrum(timestep, mass.lim, x.lim, y.lim,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{timestep}{an object of class \code{time.step}}
  \item{mass.lim}{an optional numeric vector of length 1 or 2 specifying the end points of the mass range to plot}
  \item{x.lim}{an optional numeric vector of length 1 or 2 specifying the end points of the x space range to plot}
  \item{y.lim}{an optional numeric vector of length 1 or 2 specifying the end points of the y space range to plot}
  \item{...}{standard graphical arguments to be passed to lower level plot and persp functions}
}
\details{
If no ranges are given for the mass, x and y values then the maximum possible ranges will be used. If the spatial dimension is 2 then this will produce an error (4D plotting is not supported by R yet)
Specifying single values for any of these ranges will reduce the dimensions of the plot. I.e given a 1D spatial time.step object the default plot will be a 3D persp plot with mass and x-space
on the bottom and (log) abundance vertically. If a single x-space value is given then the function will plot the size spectrum for that spatial point.
}
\value{
The function does not return anything and is used for its side effect of producing the plots.
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
