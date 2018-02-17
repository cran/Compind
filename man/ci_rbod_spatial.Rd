\name{ci_rbod_spatial}
\alias{ci_rbod_spatial}

\title{Spatial robust Benefit of the Doubt approach (Sp-RBoD)}
\description{The Spatial robust Benefit of the Doubt approach (Sp-RBoD) method allows to take into account the spatial contextual condition into the robust Benefit of the Doubt method.}

\usage{ci_rbod_spatial(x, indic_col, M=20, B=100, W) }

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{M}{The number of elements in each of the bootstrapped samples; default is 20.}
  \item{B}{The number of bootstrap replicates; default is 100.}
  \item{W}{The spatial weights matrix. A square non-negative matrix with no NAs representing spatial weights; may be a matrix of class "sparseMatrix" (spdep package)}
  
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_rbod_spatial_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="rbod_spatial".}
}

\references{
Fusco E., Vidoli F., Sahoo B.K. (2018) "\emph{Spatial heterogeneity in composite indicator: a methodological proposal}", Omega, Vol. 77, pp. 1-14
}

\author{Fusco E., Vidoli F.}

\seealso{ \code{\link{ci_rbod}}}

\examples{

data(EU_NUTS1)

coord = EU_NUTS1[,c("Long","Lat")]
k<-knearneigh(as.matrix(coord), k=5)
k_nb<-knn2nb(k)
W_mat <-nb2mat(k_nb,style="W",zero.policy=TRUE)

CI = ci_rbod_spatial(EU_NUTS1,c(2:3),M=10,B=20, W=W_mat)
}
