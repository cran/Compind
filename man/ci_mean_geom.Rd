\name{ci_mean_geom}
\alias{ci_mean_geom}
\title{Weighting method based on geometric aggregation}

\description{Geometric aggregation lets to bypass the full compensability hypothesis using geometric mean.}

\usage{ci_mean_geom(x, indic_col, na.rm=TRUE)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{na.rm}{Remove NA values before processing; default is TRUE.}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_mean_geom_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="mean_geom".}
}

\references{OECD, \emph{Handbook on constructing composite indicators: methodology and user guide}, 2008.}

\author{
  Vidoli F.
}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_factor}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_mean_geom(Indic)

data(EU_NUTS1)
CI = ci_mean_geom(EU_NUTS1,c(2:3))
}