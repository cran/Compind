\name{plot_M_robust_ci}
\alias{plot_M_robust_ci}
\title{Visualize the percentage of outperforming units in robust methods}
\description{This function generates a plot showing the proportion of outperforming units (i.e., units with a score greater than or equal to 1) across varying subset sizes (M) used in robust Benefit of the Doubt (BoD) models. It helps to analyze how the choice of subset size impacts the number of outperforming units identified by the model. For \code{method = "rbod_mdir"}, whose scores are bounded by 1, the plot shows the proportion of units on the robust frontier (score equal to 1).}

\usage{plot_M_robust_ci(x,indic_col,method,mvector,B,dir, interval)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{method}{The method to be tested. Implemented methods are "rbod", "rbod_dir", "rbod_mdir".}
  \item{mvector}{The vector of different elements in each of the bootstrapped samples to be tested.}
  \item{B}{The number of bootstrap replicates.}
  \item{dir}{Main direction. For example you can set the average rates of substitution.}
  \item{interval}{Only for \code{method = "rbod_mdir"}: significance level of the confidence intervals [see \code{\link{ci_rbod_mdir}}]; default = 0.05.}
}

\author{Fusco E.}

\seealso{ \code{\link{ci_rbod}}, \code{\link{ci_rbod_dir}}, \code{\link{ci_rbod_mdir}}}

\examples{
data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,c(2:3),polarity = c("POS","POS"), method=2)
plot_M_robust_ci(data_norm$ci_norm, c(1:2), method = "rbod", 
                 mvector = seq(5,25,10), B = 10)
plot_M_robust_ci(data_norm$ci_norm, c(1:2), method = "rbod_dir", 
                 mvector = seq(5,25,10), B = 10, dir=c(1,0.1))


\donttest{
data(BLI_2017)
plot_M_robust_ci(x = BLI_2017, indic_col = c(2:12), method = "rbod_mdir", 
                 mvector = c(5, 10), B = 5, interval = 0.05)
}
}
