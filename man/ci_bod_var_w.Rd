\name{ci_bod_var_w}
\alias{ci_bod_var_w}
\title{Variance weighted Benefit of the Doubt approach (BoD variance weighted)}
\description{Variance weighted Benefit of the Doubt approach (BoD variance weighted) is a particular form of BoD method with additional information in the optimization problem. In particular it has been added weight constraints (in form of an Assurance region type I (AR I)) endogenously determined in order to take into account the ratio of the vertical variability of each simple indicator relative to one another.}

\usage{ci_bod_var_w(x,indic_col,boot_rep = 5000)}

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{boot_rep}{The number of bootstrap replicates (default=5000) for the estimates of the nonparametric bootstrap (first order normal approximation) confidence intervals for the variances of the simple indicators.}
}

\details{For more informations about the estimation of the confidence interval for the variances, please see function \emph{boot.ci}, package \emph{boot}.}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_bod_var_w_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="bod_var_w".}
}

\references{Vidoli F., Mazziotta C. (2013) "\emph{Robust weighted composite indicators by means of frontier methods with an application to European infrastructure endowment}", Statistica Applicata, Italian Journal of Applied Statistics.}

\author{Vidoli F.}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_rbod}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_bod_var_w(Indic)

}
