\name{ci_factor_mixed}
\alias{ci_factor_mixed}
\title{Weighting method based on Factor analysis of mixed data (FAMD)}

\description{Factor analysis of mixed data (FAMD) can be seen as a principal component method dedicated to analyze a data set containing both quantitative and qualitative variables making possible to compute composite indicators taking into account continous, dummy, or factor variables}

\usage{
ci_factor_mixed(x,indic_col,method="ONE",dim)
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
  \item{ci_method}{Method used; for this function ci_method="factor_mixed".}
}


\author{Luis Carlos Castillo Tellez}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_factor}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
i3 <- seq(0, 1, len = 100)
i3 = as.factor(ifelse(i3>0.5,1,0))
Indic = data.frame(i1, i2, i3)

CI  = ci_factor_mixed(Indic,c(1:3))
CI2 = ci_factor_mixed(Indic,c(1:3), method="ALL")
CI3 = ci_factor_mixed(Indic,c(1:3), method="CH", dim=2)
}
