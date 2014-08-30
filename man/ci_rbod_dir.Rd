\name{ci_rbod_dir}
\alias{ci_rbod_dir}

\title{Directional Robust Benefit of the Doubt approach (D-RBoD)}

\description{Directional Robust Benefit of the Doubt approach (D-RBoD) is the directional robust version of the BoD method.}

\usage{ci_rbod_dir(x,indic_col,M,B,dir)}

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{M}{The number of elements in each of the bootstrapped samples.}
  \item{B}{The number of bootstap replicates.}
  \item{dir}{Main direction. For example you can set the average rates of substitution.}  
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_rbod_dir_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="rbod_dir".}
}

\references{

Daraio C., Simar L., "\emph{Introducing environmental variables in nonparametric frontier models: a probabilistic approach}", Journal of productivity analysis, 2005, 24(1), 93 121.

Simar L., Vanhems A., "\emph{Probabilistic characterization of directional distances and their robust versions}", Journal of Econometrics, 2012, 166(2), 342 354.

Vidoli F., Fusco E., Mazziotta C., "\emph{Non-compensability in composite indicators: a robust directional frontier method}", Social Indicators Research, Springer Netherlands.

}

\author{Fusco E., Vidoli F.}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_rbod}}}

\examples{
data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,c(2:3),polarity = c("POS","POS"), method=2)
CI = ci_rbod_dir(data_norm$ci_norm, c(1:2), M = 25, B = 50, c(1,0.1)) 
}
