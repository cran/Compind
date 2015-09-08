\name{ci_mean_min}
\alias{ci_mean_min}
\title{Mean-Min Function}
\description{The Mean-Min Function (MMF) is an intermediate case between arithmetic \code{mean}, according to which no unbalance is penalized, and \code{min} function, according to which the penalization is maximum. It depends on two parameters that are respectively related to the intensity of penalization of unbalance (\eqn{\alpha}) and intensity of complementarity (\eqn{\beta}) among indicators.}

\usage{ci_mean_min(x, indic_col, alpha, beta)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{alpha}{The intensity of penalisation of unbalance among indicators, \eqn{0 \le \alpha \le 1}}
  \item{beta}{The intensity of complementarity among indicators, \eqn{\beta \ge 0}}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_mean_min_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="mean_min".}
}

\references{Casadio Tarabusi, E., & Guarini, G. (2013) "\emph{An unbalance adjustment method
for development indicators}", Social indicators research, 112(1), 19-45.}

\author{Vidoli F.}

\seealso{ \code{\link{ci_mpi}}, \code{\link{normalise_ci}}}

\examples{
data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,c(2:3),c("NEG","POS"),method=2)
CI = ci_mean_min(data_norm$ci_norm, alpha=0.5, beta=1)

}