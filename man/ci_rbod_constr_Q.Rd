\name{ci_rbod_constr_Q}
\alias{ci_rbod_constr_Q}
\title{Conditional robust constrained Benefit of the Doubt approach (BoD)}
\description{The Conditional robust constrained Benefit of the Doubt function introduces additional constraints to the weight variation in the optimization procedure (Constrained Virtual Weights Restriction) allowing to restrict the importance attached to a single indicator expressed in percentage terms, ranging between a lower and an upper bound (VWR); this function, furthermore, allows to calculate the composite indicator allowing to impose a preference structure (ordVWR). This function, in addition to being robust against outlier data (see \code{ci_rbod_constr_bad} function) allows to take into account external contextual continuous (\code{Q}) or/and ordinal (\code{Q_ord}) variables.}

\usage{ci_rbod_constr_Q(x, indic_col, 
low_w=0, pref=NULL, M, B, Q=NULL, Q_ord=NULL, bandwidth)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{low_w}{Importance weights lower bound.}
  \item{pref}{The preference vector among indicators; For example if \code{Indic1} is the most important, \code{Indic2},\code{Indic3} are more important than \code{Indic4} and no preference judgment on \code{Indic5} (= not included in the vector), the \code{pref} vector can be written as: \code{c("Indic1", "Indic2","Indic3","Indic4")}}
  \item{M}{The number of elements in each of the bootstrapped samples.}
  \item{B}{The number of bootstrap replicates.}
  \item{Q}{A matrix containing continuous exogenous variables.}
  \item{Q_ord}{A matrix containing discrete exogenous variables.}
  \item{bandwidth}{Multivariate mixed bandwidth for exogenous variables; it can be calculated by \code{bandwidth_CI} function.}
  
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_rbod_constr_Q_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="rbod_constr_Q".}
  \item{ci_rbod_constr_Q_weights}{Raw weights assigned to each simple indicator.}
  \item{ci_rbod_constr_Q_target}{Indicator target values.}
}

\references{
Rogge N., de Jaeger S. and Lavigne C. (2017) "\emph{Waste Performance of NUTS 2-regions in the EU: A Conditional Directional Distance Benefit-of-the-Doubt Model}", Ecological Economics, vol.139, pp. 19-32.

Zanella A., Camanho A.S. and Dias T.G. (2015) "\emph{Undesirable outputs and weighting schemes in composite indicators based on data envelopment analysis}", European Journal of Operational Research, vol. 245(2), pp. 517-530.
}

\author{Fusco E., Rogge N., Vidoli F.}

\seealso{ \code{\link{ci_rbod_constr_bad}}, \code{\link{ci_bod_constr_bad}} }



\examples{
data(EU_2020)

indic <- c("employ_2011", "gasemiss_2011","deprived_2011")  
dat   <- EU_2020[-c(10,18),indic]
Q_GDP <- EU_2020[-c(10,18),"percGDP_2011"]

# Conditional robust BoD Constrained VWR
band = bandwidth_CI_bad(dat, ngood=1, nbad=2, Q = Q_GDP)

CI_BoD_C = ci_rbod_constr_Q(dat, 
                                low_w=0.05, 
                                pref=NULL, 
                                M=10, 
                                B=50, 
                                Q=Q_GDP, 
                                bandwidth = band$bandwidth)
CI_BoD_C$ci_rbod_constr_Q_est

}
