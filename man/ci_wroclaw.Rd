\name{ci_wroclaw}
\alias{ci_wroclaw}
\title{Wroclaw Taxonomic Method}

\description{Wroclaw taxonomy method (also known as the dendric method), originally developed at the University of Wroclaw, is based on the distance from a theoretical unit characterized by the best performance for all indicators considered; the composite indicator is therefore based on the sum of euclidean distances from the ideal unit and normalized by a measure of variability of these distance (mean + 2*std).}

\usage{ci_wroclaw(x,indic_col)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{Simple indicators column number.}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_wroclaw_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="wroclaw".}
}

\details{Please pay attention that \emph{ci_wroclaw_est} is the distance from the "ideal" unit; so, units with higher values for the simple indicators get lower values of composite indicator.}

\references{

UNESCO, "\emph{Social indicators: problems of definition and of selection}", Paris 1974.

Mazziotta C., Mazziotta M., Pareto A., Vidoli F., "\emph{La sintesi di indicatori territoriali di dotazione
infrastrutturale: metodi di costruzione e procedure di ponderazione a confronto}", Rivista di
Economia e Statistica del territorio, n.1, 2010.
}

\author{ Vidoli F.}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_mpi}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_wroclaw(Indic)

data(EU_NUTS1)
CI = ci_wroclaw(EU_NUTS1,c(2:3))

data(EU_2020)
data_selez = EU_2020[,c(1,22,191)]
data_norm = normalise_ci(data_selez,c(2:3),c("POS","NEG"),method=3)
ci_wroclaw(data_norm$ci_norm,c(1:2))

}