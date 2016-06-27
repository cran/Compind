\name{EU_NUTS1}
\alias{EU_NUTS1}
\title{EU NUTS1 Transportation data}
\description{Eurostat regional transport statistics (reg_tran) data, year 2012.

For more info, please see \url{http://ec.europa.eu/eurostat/data/browse-statistics-by-theme}.
}
\usage{data(EU_NUTS1)}
\format{EU_NUTS1 is a dataset with 34 observations and two indicators describing transportation infrastructure endowment of the main (in terms of population and GDP) European NUTS1 regions: France, Germany, Italy, Spain (United Kingdom has been omitted, due to lack of data concerning railways).

\describe{
    \item{roads}{Calculated as (2 * Motorways - Kilometres per 1000 km2 + Other roads - Kilometres per 1000 km2 )/3}	
    \item{trains}{Calculated as (2 *Railway lines double+Electrified railway lines)/3}	
}}

\references{
 Vidoli F., Mazziotta C., "\emph{Robust weighted composite indicators by means of frontier methods with an application to European infrastructure endowment}", Statistica Applicata, Italian Journal of Applied Statistics, 2013.}
\author{Vidoli F.}
\examples{data(EU_NUTS1)}