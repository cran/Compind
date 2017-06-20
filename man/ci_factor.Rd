\name{ci_factor}
\alias{ci_factor}
\title{Weighting method based on Factor Analysis}

\description{Factor analysis groups together collinear simple indicators to estimate a composite indicator that captures as much as possible of the information common to individual indicators.}

\usage{
ci_factor(x,indic_col,method="ONE",dim)
}

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{method}{If method = "ONE" (default) the composite indicator estimated values are equal to first component scores; if method = "ALL" the composite indicator estimated values are equal to component score multiplied by its proportion variance;  if method = "CH" it can be choose the number of the component to take into account.}
  \item{dim}{Number of chosen component (if method = "CH", default is 3).}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_factor_est}{Composite indicator estimated values.}
  \item{loadings_fact}{Variance explained by principal factors (in percentage terms).}
  \item{ci_method}{Method used; for this function ci_method="factor".}
}

\references{OECD (2008) "\emph{Handbook on constructing composite indicators: methodology and user guide}".}

\author{
  Vidoli F.
}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_mpi}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_factor(Indic)

data(EU_NUTS1)
CI = ci_factor(EU_NUTS1,c(2:3), method="ALL")

data(EU_2020)
data_norm = normalise_ci(EU_2020,c(47:51),polarity = c("POS","POS","POS","POS","POS"), method=2)
CI3 = ci_factor(data_norm$ci_norm,c(1:5),method="CH", dim=3)
}
