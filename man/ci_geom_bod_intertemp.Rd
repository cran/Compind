\name{ci_geom_bod_intertemp}
\alias{ci_geom_bod_intertemp}
\title{Intertemporal analysis for geometric mean quantity index numbers}

\description{Intertemporal analysis for geometric mean quantity index numbers with Benefit-of-the-Doubt weights - see function \code{\link{ci_bod_constr}}.}

\usage{ci_geom_bod_intertemp(x0,x1,indic_col,up_w,low_w,bench)}

\arguments{
  \item{x0}{A data.frame containing simple indicators - time 0}
  \item{x1}{A data.frame containing simple indicators - time 1}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{up_w}{Weights upper bound.}
  \item{low_w}{Weights lower bound.}
  \item{bench}{Row number of the benchmark unit}
}

%\details{...}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_geom_bod_intertemp_est}{A matrix containing the Overall Change (period t1 vs t0), the Change Effect (period t1 vs t0), the Benchmark Effect (period t1 vs t0) and Weight Effect (period t1 vs t0).}
  \item{ci_method}{Method used; for this function ci_method="Intertemporal_effects_Geometric_BoD".}
}

\references{ Van Puyenbroeck T. and Rogge N. (2017) "\emph{Geometric mean quantity index numbers with Benefit-of-the-Doubt weights}", European Journal of Operational Research, Volume 256, Issue 3, Pages 1004 - 1014}

\author{Rogge N., Vidoli F.}

\seealso{ \code{\link{ci_bod_constr}},\code{\link{ci_bod}} }

\examples{

i1_t1 <- seq(0.3, 0.5, len = 100)
i2_t1 <- seq(0.3, 1, len = 100)  
Indic_t1 = data.frame(i1_t1, i2_t1)

i1_t0 <- i1_t1 - rnorm (100, 0.2, 0.03)
i2_t0 <- i2_t1   - rnorm (100, 0.2, 0.03)
Indic_t0 = data.frame(i1_t0, i2_t0)

intertemp = ci_geom_bod_intertemp(Indic_t0,Indic_t1,c(1:2),up_w=0.95,low_w=0.05,1)
intertemp 

}
