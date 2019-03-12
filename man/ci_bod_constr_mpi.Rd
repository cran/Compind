\name{ci_bod_constr_mpi}
\alias{ci_bod_constr_mpi}
\title{Non Compensative Constrained Benefit of the Doubt approach (BoD)}

\description{The constrained Benefit of the Doubt function lets to introduce additional constraints to the weight variation in the optimization procedure so that all the weights obtained are greater than a lower value (low_w) and less than an upper value (up_w). In a second step the composite indicator is adjusted by a "penalty" coefficient related to the variability of each unit (see, method of the coefficient of variation penalty - ci_mpi)}

\usage{ci_bod_constr_mpi(x,indic_col,up_w,low_w,penalty="POS")}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{up_w}{Importance weights upper bound.}
  \item{low_w}{Importance weights lower bound.}
  \item{penalty}{Penalty direction; Use "POS" (default) in case of 'increasing' or 'positive' composite index (e.g., well-being index)), "NEG" in case of 'decreasing' or 'negative' composite index (e.g., poverty index).}
}

%\details{...}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_bod_constr_est_mpi}{Constrained composite indicator estimated values.}
  \item{ci_bod_constr_mpi_pen}{Penalized constrained composite indicator estimated values.}
    \item{ci_bod_constr_mpi_weights}{Raw constrained weights (not penalized) assigned to the simple indicators.}
  \item{ci_method}{Method used; for this function ci_method="bod_constrained_mpi".}
}

\references{ Van Puyenbroeck T. and Rogge N. (2017) "\emph{Geometric mean quantity index numbers with Benefit-of-the-Doubt weights}", European Journal of Operational Research, Volume 256, Issue 3, Pages 1004 - 1014.

De Muro P., Mazziotta M., Pareto A. (2011), "\emph{Composite Indices of Development and Poverty: An Application to MDGs}", Social Indicators Research, Volume 104, Number 1, pp. 1-18.
}

\author{Vidoli F.}

\seealso{ \code{\link{ci_bod_constr}},\code{\link{ci_mpi}} }

\examples{
data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,
                         c(2:3),
                         c("NEG","POS"),
                         method=1,
                         z.mean=100,
                         z.std=10)

CI = ci_bod_constr_mpi(data_norm$ci_norm,
                   c(1:2),
                   up_w=1,
                   low_w=0.1,
                   penalty="POS")
                   
  }
