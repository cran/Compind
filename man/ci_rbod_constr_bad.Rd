\name{ci_rbod_constr_bad}
\alias{ci_rbod_constr_bad}
\title{Robust constrained Benefit of the Doubt approach (BoD) in presence of undesirable  indicators}
\description{The Robust constrained Benefit of the Doubt function introduces additional constraints to the weight variation in the optimization procedure (Constrained Virtual Weights Restriction) allowing to restrict the importance attached to a single indicator expressed in percentage terms, ranging between a lower and an upper bound (VWR); this function, furthermore, allows to calculate the composite indicator simultaneously in presence of undesirable (bad) and desirable (good) indicators allowing to impose a preference structure (ordVWR). This function is the robust version of the \code{ci_bod_constr_bad}: it is based on the concept of the expected minimum input function of order-\emph{m} (Daraio and Simar, 2005) allowing to compare the unit under analysis against \code{M} peers by extracting \code{B} samples with replacement.}

\usage{ci_rbod_constr_bad(x, indic_col, ngood=1, nbad=1, low_w=0, pref=NULL, M, B) }

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{ngood}{The number of desirable outputs; it has to be greater than 0.}
  \item{nbad}{The number of undesirable  outputs; it has to be greater than 0.}
  \item{low_w}{Importance weights lower bound.}
  \item{pref}{The preference vector among indicators; For example if \code{Indic1} is the most important, \code{Indic2},\code{Indic3} are more important than \code{Indic4} and no preference judgment on \code{Indic5} (= not included in the vector), the \code{pref} vector can be written as: \code{c("Indic1", "Indic2","Indic3","Indic4")}}
  \item{M}{The number of elements in each of the bootstrapped samples.}
  \item{B}{The number of bootstrap replicates.}
  
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_rbod_constr_bad_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="rbod_constr_bad".}
  \item{ci_rbod_constr_bad_weights}{Raw weights assigned to each simple indicator.}
  \item{ci_rbod_constr_bad_target}{Indicator target values.}
}

\references{
Rogge N., de Jaeger S. and Lavigne C. (2017) "\emph{Waste Performance of NUTS 2-regions in the EU: A Conditional Directional Distance Benefit-of-the-Doubt Model}", Ecological Economics, vol.139, pp. 19-32.

Zanella A., Camanho A.S. and Dias T.G. (2015) "\emph{Undesirable outputs and weighting schemes in composite indicators based on data envelopment analysis}", European Journal of Operational Research, vol. 245(2), pp. 517-530.
}

\author{Fusco E., Rogge N.}

\seealso{ \code{\link{ci_bod_constr}}, \code{\link{ci_bod_constr_bad}} }



\examples{
data(EU_2020)
indic <- c("employ_2011", "percGDP_2011", "gasemiss_2011","deprived_2011")  
dat <- EU_2020[-c(10,18),indic]

# Robust BoD Constrained VWR
CI_BoD_C = ci_rbod_constr_bad(dat, ngood=2, nbad=2, low_w=0.05, pref=NULL, M=10, B=50)
CI_BoD_C$ci_rbod_constr_bad_est

# Robust BoD Constrained ordVWR
importance <- c("gasemiss_2011","percGDP_2011","employ_2011")
CI_BoD_C = ci_rbod_constr_bad(dat, ngood=2, nbad=2, low_w=0.05, pref=importance, M=10, B=50)
CI_BoD_C$ci_rbod_constr_bad_est
}
