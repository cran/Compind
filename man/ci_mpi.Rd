\name{ci_mpi}
\alias{ci_mpi}
\title{Mazziotta-Pareto Index (MPI) method}
\description{Mazziotta-Pareto Index (MPI) is a non-linear composite index method which transforms a set of individual indicators in standardized variables and summarizes them using an arithmetic mean adjusted by a "penalty" coefficient related to the variability of each unit (method of the coefficient of variation penalty).}

\usage{ci_mpi(x, indic_col, penalty="POS")}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{penalty}{Penalty direction; Use "POS" (default) in case of 'increasing' or 'positive' composite index (e.g., well-being index)), "NEG" in case of 'decreasing' or 'negative' composite index (e.g., poverty index).}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_mpi_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="mpi".}
}

\references{
De Muro P., Mazziotta M., Pareto A. (2011), "\emph{Composite Indices of Development and Poverty: An Application to MDGs}", Social Indicators Research, Volume 104, Number 1, pp. 1-18.}

\author{Vidoli F.}

\seealso{ \code{\link{ci_bod}}, \code{\link{normalise_ci}}}

\examples{
data(EU_NUTS1)

# Please, pay attention. MPI can be calculated only with two standardizations methods:
# Classic MPI - method=1, z.mean=100 and z.std=10
# Correct MPI - method=1, z.mean=70  and z.std=130
# For more info, please see references.

data_norm = normalise_ci(EU_NUTS1,c(2:3),c("NEG","POS"),method=1,z.mean=0, z.std=1)
CI = ci_mpi(data_norm$ci_norm, penalty="NEG")

data(EU_NUTS1)
CI = ci_mpi(EU_NUTS1,c(2:3),penalty="NEG")
}