\name{ci_rbod}
\alias{ci_rbod}

\title{Robust Benefit of the Doubt approach (RBoD)}
\description{Robust Benefit of the Doubt approach (RBoD) is the robust version of the BoD method. It is based on the concept of the expected minimum input function of order-\emph{m} so "\emph{in place of looking for the lower boundary of the support of F, as was typically the case for the full-frontier (DEA or FDH), the order-m efficiency score can be viewed as the expectation of the maximal score, when compared to m units randomly drawn from the population of units presenting a  greater level of simple indicators}", Daraio and Simar (2005).}

\usage{ci_rbod(x,indic_col,M,B)}

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{M}{The number of elements in each of the bootstrapped samples.}
  \item{B}{The number of bootstrap replicates.}
  
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_rbod_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="rbod".}
}

\references{

Daraio, C., Simar, L. "\emph{Introducing environmental variables in nonparametric frontier models: a probabilistic approach}", Journal of productivity analysis, 2005, 24(1), 93 - 121.

Vidoli F., Mazziotta C., "\emph{Robust weighted composite indicators by means of frontier methods with an application to European infrastructure endowment}", Statistica Applicata, Italian Journal of Applied Statistics, 2013.

}

\author{Vidoli F.}

\seealso{ \code{\link{ci_bod}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_rbod(Indic,B=10)

data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,c(2:3),polarity = c("POS","POS"), method=2)
CI = ci_rbod(data_norm$ci_norm,c(1:2),M=10,B=20)
}
