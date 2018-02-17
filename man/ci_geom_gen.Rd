\name{ci_geom_gen}
\alias{ci_geom_gen}
\title{Generalized geometric mean quantity index numbers}

\description{This function use the geometric mean to aggregate the single indicators. Two weighting criteria has been implemented: EQUAL: equal weighting and BOD: Benefit-of-the-Doubt weights following the Puyenbroeck and Rogge (2017) approach.}

\usage{ci_geom_gen(x,indic_col,meth,up_w,low_w,bench)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{meth}{"EQUAL" = Equal weighting set, "BOD" =  Benefit-of-the-Doubt weighting set.}
  \item{up_w}{if meth="BOD"; upper bound of the weighting set.}
  \item{low_w}{if meth="BOD"; lower bound of the weighting set.}
  \item{bench}{Row number of the benchmark unit used to normalize the data.frame \code{x}.}
}


\value{An object of class "CI". This is a list containing the following elements:

If meth = "EQUAL":

\itemize{

  \item{ci_mean_geom_est}{: Composite indicator estimated values.}
  \item{ci_method}{: Method used; for this function ci_method="mean_geom".}
}

If meth = "BOD":  
\itemize{

  \item{ci_geom_bod_est}{: Constrained composite indicator estimated values.}
\item{ci_geom_bod_weights}{: Raw constrained weights assigned to the simple indicators.}
  \item{ci_method}{: Method used; for this function ci_method="geometric_bod".}
  }

}

\references{ Van Puyenbroeck T. and Rogge N. (2017) "\emph{Geometric mean quantity index numbers with Benefit-of-the-Doubt weights}", European Journal of Operational Research, Volume 256, Issue 3, Pages 1004 - 1014}

\author{Rogge N., Vidoli F.}

\seealso{ \code{\link{ci_bod_dir}},\code{\link{ci_bod}} }

\examples{
i1 <- seq(0.3, 1, len = 100) - rnorm (100, 0.1, 0.03)
i2 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.1, 0.03)
i3 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.1, 0.03)
Indic = data.frame(i1, i2,i3)

geom1 = ci_geom_gen(Indic,c(1:3),meth = "EQUAL")
geom1$ci_mean_geom_est
geom1$ci_method


geom2 = ci_geom_gen(Indic,c(1:3),meth = "BOD",0.7,0.3,100)
geom2$ci_geom_bod_est
geom2$ci_geom_bod_weights


}
