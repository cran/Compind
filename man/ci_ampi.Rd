\name{ci_ampi}
\alias{ci_ampi}
\title{Adjusted Mazziotta-Pareto Index (AMPI) method}
\description{Adjusted Mazziotta-Pareto Index (AMPI) is a non-compensatory composite index that allows to take into account the time dimension, too. The calculation part is similat to the MPI framework, but the standardization part make the scores obtained over the years comparable.}

\usage{ci_ampi(x, indic_col, gp, time, polarity, penalty = "NEG")}

\arguments{
  \item{x}{A data.frame containing simple indicators in a Long Data Format.}
  \item{indic_col}{Simple indicators column number.}
  \item{gp}{Goalposts; to facilitate the interpretation of results, the \code{goalposts} can be chosen so that 100 represents a reference value (e.g., the average in a given year).}
  \item{time}{The time variable (mandatory); if the analysis is carried out over a single year, it is necessary to create a constant variable (i.e. \code{dataframe@year <- 2014}).}
  \item{polarity}{Polarity vector: "POS" = positive, "NEG" = negative. The polarity of a individual indicator is the sign of the relationship between the indicator and the phenomenon to be measured (e.g., in a well-being index, "GDP per capita" has 'positive' polarity and "Unemployment rate" has 'negative' polarity).}
  \item{penalty}{Penalty direction; Use "NEG" (default) in case of 'increasing' or 'positive' composite index (e.g., well-being index)), "POS" in case of 'decreasing' or 'negative' composite index (e.g., poverty index).}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_ampi_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="ampi".}
}

\details{Author thanks Leonardo Alaimo for their help and for making available the original code of the AMPI function.}

\references{
Mazziotta, M., Pareto, A. (2013) "\emph{A Non-compensatory Composite Index for Measuring Well-being over Time}", Cogito. Multidisciplinary Research Journal Vol. V, no. 4, pp. 93-104

Mazziotta, M., Pareto, A. (2016)."\emph{On a Generalized Non-compensatory Composite Index for Measuring Socio-economic Phenomena}", Cogito. Social Indicators Research, Vol. 127, no. 3, pp. 983-1003
}

\author{Fusco E., Alaimo L., Giovagnoli C., Patelli L.}

\seealso{ \code{\link{ci_bod}}, \code{\link{normalise_ci}}}

\examples{
data(EU_2020)

data_test = EU_2020[,c("employ_2010","employ_2011","finalenergy_2010","finalenergy_2011")] 

EU_2020_long<-reshape(data_test, 
                      varying=c("employ_2010","employ_2011","finalenergy_2010","finalenergy_2011"), 
                      direction="long", 
                      idvar="geo", 
                      sep="_")

CI <- ci_ampi(EU_2020_long, 
              indic_col=c(2:3),
              gp=c(50, 100), 
              time=EU_2020_long[,1], 
              polarity= c("POS", "POS"), 
              penalty="POS")
CI$ci_ampi_est

}