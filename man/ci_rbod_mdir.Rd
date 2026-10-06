\name{ci_rbod_mdir}
\alias{ci_rbod_mdir}
\title{Robust multi-directional Benefit of the Doubt approach (MDRBoD)}
\description{Robust Multi-directional Benefit of the Doubt (MDRBoD) allows to introduce the non-compensability among simple indicators in a standard Robust BOD in an objective manner: the preference structure, i.e., the direction, is determined directly from the data and is specific for each unit and these estimated values are calculated as the reference sample varies in order to smooth out the effect of outliers or out-of-range data.}

\usage{ci_rbod_mdir(x,indic_col,M, B, interval)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{M}{The number of elements in each of the bootstrapped samples (M lower than the number of units).}
  \item{B}{The number of bootstrap replicates.}
  \item{interval}{Significance level of the confidence intervals of the estimated scores; default = 0.05, i.e. 95\% confidence intervals.}
}

\details{At each replicate \eqn{b = 1, \dots, B} a subsample of \code{M} units is drawn with replacement and every unit is evaluated against the order-\eqn{m} frontier spanned by the subsample, as in Vidoli et al. (2024), eq. 3-6: the potential improvements of each simple indicator define the unit-specific direction and \eqn{\beta \in [0,1]} measures the proportion of the improvement needed to reach the frontier (the evaluated unit is included in its own reference set, so that a unit dominating the subsample lies on the frontier, \eqn{\beta = 0}). Directions and \eqn{\beta} are averaged over the \eqn{B} replicates and the composite indicator is computed with the averaged quantities (eq. 7-10); scores are therefore bounded by 1 and converge to \code{\link{ci_bod_mdir}} as \code{M} grows.

Two confidence intervals are returned. \code{conf} follows eq. 11 (Student-t, \eqn{\pm t \cdot s/\sqrt{B}}, where \eqn{s} is the standard deviation of the replicate scores), centred on the estimated score: it measures the Monte Carlo precision of the estimate and narrows as \eqn{B} increases. \code{conf_perc} contains the percentile interval of the replicate scores, which describes the variability of the score with respect to the reference subsample and does not shrink with \eqn{B}.}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_rbod_mdir_est}{Composite indicator estimated values (eq. 10).}
  \item{conf}{lower_ci and upper_ci; Student-t confidence interval centred on the estimated values (eq. 11).}
  \item{conf_perc}{lower_ci and upper_ci; percentile interval of the replicate scores.}
  \item{ci_rbod_mdir_spec}{Simple indicators specific scores (eq. 9).}
  \item{ci_rbod_mdir_dir}{Directions for each simple indicator and unit, averaged over the replicates (eq. 7).}
  \item{ci_rbod_mdir_beta}{Proportion of the improvement needed to reach the frontier, averaged over the replicates (eq. 8).}
  \item{ci_rbod_mdir_boot}{Matrix (units x replicates) of the replicate scores.}
  \item{ci_method}{Method used; for this function ci_method="rbod_mdir".}
}

\references{

F. Vidoli, E. Fusco, G. Pignataro, C. Guccio (2024) "\emph{Multi-directional Robust Benefit of the Doubt model: An application to the measurement of the quality of acute care services in OECD countries}", Socio-Economic Planning Sciences. https://doi.org/10.1016/j.seps.2024.101877  

}

\author{Vidoli F.}

\seealso{ \code{\link{ci_rbod}}, \code{\link{ci_bod_mdir}}}

\examples{
data(BLI_2017)
CI <- ci_rbod_mdir(BLI_2017, c(2:5), M = 10, B = 10, interval = 0.05)
CI$ci_rbod_mdir_est
CI$conf_perc
\donttest{
CI <- ci_rbod_mdir(BLI_2017, c(2:12), M = 6, B = 20, interval = 0.05)
}
}
