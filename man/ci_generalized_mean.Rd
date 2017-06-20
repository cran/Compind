\name{ci_generalized_mean}
\alias{ci_generalized_mean}
\title{Weighting method based on generalized mean}

\description{Generalized means are a family of functions for aggregating sets of numbers (it include as special cases the Pythagorean means, arithmetic, geometric, and harmonic means). The generalized mean is also known as power mean or Holder mean.}

\usage{ci_generalized_mean(x, indic_col, p, na.rm=TRUE)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{p}{Exponent \eqn{p} (real number).}
  \item{na.rm}{Remove NA values before processing; default is TRUE.}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_generalized_mean_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="generalized_mean".}
}

%\references{??}

\note{
The \emph{generalized mean} with the exponent \eqn{p} can be espressed as:
\deqn{M_p(I_1,\dots,I_n) = \left( \frac{1}{n} \sum_{i=1}^n I_i^p \right)^{\frac{1}{p}}}
Particular case are:
\eqn{p=-\infty}: minimum, \eqn{p=-1}: harmonic mean, \eqn{p=0}: geometric mean, 
\eqn{p=1}: arithmetic mean, \eqn{p=2}: root-mean-square and \eqn{p=\infty}: maximum.
}

\author{
  Vidoli F.
}

\seealso{ \code{\link{ci_geom_gen}}, \code{\link{ci_factor}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_generalized_mean(Indic, p=-1) # harmonic mean

data(EU_NUTS1)
CI = ci_generalized_mean(EU_NUTS1,c(2:3),p=2) # geometric mean
}