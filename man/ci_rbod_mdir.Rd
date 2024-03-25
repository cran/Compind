\name{ci_rbod_mdir}
\alias{ci_rbod_mdir}
\title{Robust multi-directional Benefit of the Doubt approach (MDRBoD)}
\description{Robust Multi-directional Benefit of the Doubt (MDRBoD) allows to introduce the non-compensability among simple indicators in a standard Robust BOD in an objective manner: the preference structure, i.e., the direction, is determined directly from the data and is specific for each unit and these estimated values are calculated as the reference sample varies in order to smooth out the effect of outliers or out-of-range data.}

\usage{ci_rbod_mdir(x,indic_col,M, B, interval)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{M}{The number of elements in each of the bootstrapped samples.}
  \item{B}{The number of bootstrap replicates.}
  \item{interval}{Desired probability for Student distribution [see function qt()]; default = 0.05.}

}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_rbod_mdir_est}{Composite indicator estimated values.}
  \item{conf}{lower_ci and upper_ci; Estimated confidence interval for the composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="rbod_mdir".}
  \item{ci_rbod_mdir_spec}{Simple indicators specific scores.}
  \item{ci_rbod_mdir_dir}{Directions for each simple indicator and unit.}
  
}

\references{

F. Vidoli, E. Fusco, G. Pignataro, C. Guccio (2024) "\emph{Multi-directional Robust Benefit of the Doubt model: An application to the measurement of the quality of acute care services in OECD countries}", Socio-Economic Planning Sciences. https://doi.org/10.1016/j.seps.2024.101877  

}

\author{Vidoli F.}

\seealso{ \code{\link{ci_rbod}}, \code{\link{ci_bod_mdir}}}

\examples{
data(BLI_2017)
CI <- ci_rbod_mdir(BLI_2017,c(2:12), M=10,B=20, interval=0.05)
}
