\name{bandwidth_CI}
\alias{bandwidth_CI}

\title{Multivariate mixed bandwidth selection for exogenous variables (in presence of all desirable  indicators)}

\description{A function for the selection of optimal multivariate mixed bandwidths for the kernel
density estimation of continuous and discrete exogenous variables.}

\usage{bandwidth_CI(x, indic_col, Q=NULL, Q_ord=NULL)}

\arguments{
  \item{x}{A data frame containing simple indicators.}
  \item{indic_col}{Simple indicators column number; it has to be greater than 1}
  \item{Q}{A matrix containing continuous exogenous variables.}
  \item{Q_ord}{A matrix containing discrete exogenous variables.}
}


\value{
  \item{bandwidth}{A matrix containing the optimal bandwidths for the exogenous variables indicate in Q and Q_ord.}
  \item{ci_method}{"bandwidth_CI"}
}

\details{Author thanks Nicky Rogge for his help and for making available the original code of the bandwidth function.}

\author{Fusco E., Rogge N.}

%\seealso{ \code{\link{ci_bod}}, \code{\link{ci_mpi}}}

\examples{
data(EU_2020)
indic <- c("employ_2011", "gasemiss_2011","deprived_2011")  
dat   <- EU_2020[-c(10,18),indic]
Q_GDP <- EU_2020[-c(10,18),"percGDP_2011"]

# Conditional robust BoD Constrained VWR
band = bandwidth_CI(dat, Q = Q_GDP)
}