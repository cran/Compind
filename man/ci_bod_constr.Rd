\name{ci_bod_constr}
\alias{ci_bod_constr}
\title{Constrained Benefit of the Doubt approach (BoD)}

\description{The constrained Benefit of the Doubt function lets to introduce additional constraints to the weight variation in the optimization procedure so that all the weights obtained are greater than a lower value (low_w) and less than an upper value (up_w).}

\usage{ci_bod_constr(x,indic_col,up_w,low_w)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{up_w}{Importance weights upper bound.}
  \item{low_w}{Importance weights lower bound.}
}

%\details{...}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_bod_constr_est}{Constrained composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="bod_constrained".}
  \item{ci_bod_constr_weights}{Raw constrained weights assigned to the simple indicators.}
}

\references{ Van Puyenbroeck T. and Rogge N. (2017) "\emph{Geometric mean quantity index numbers with Benefit-of-the-Doubt weights}", European Journal of Operational Research, Volume 256, Issue 3, Pages 1004 - 1014.}

\author{Rogge N., Vidoli F.}

\seealso{ \code{\link{ci_bod_dir}},\code{\link{ci_bod}} }

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_bod_constr(Indic,up_w=1,low_w=0.05)

data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,c(2:3),polarity = c("POS","POS"), method=2)
CI = ci_bod_constr(data_norm$ci_norm,c(1:2),up_w=1,low_w=0.05)
}
