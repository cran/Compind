\name{ci_owa}
\alias{ci_owa}

\title{Ordered Weighted Average (OWA)}
\description{The Ordered Weighted Averaging (OWA) operator is a multi-criteria decision aggregation method that is structurally non-compensatory (Yager, 1988).}

\usage{ci_owa(x, id, indic_col, atleastjp)}

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
  \item{id}{Units' unique identifier.}
  \item{indic_col}{Simple indicators column number.}
  \item{atleastjp}{Fuzzy linguistic quantifier "At least j".}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{CI_OWA_n}{Composite indicator estimated values for OWA-.}
  \item{CI_OWA_p}{Composite indicator estimated values for OWA+.}
  \item{wp}{OWA weights' vector "More than j".}
  \item{wn}{OWA weights' vector "At least j".}
  \item{ci_method}{Method used; for this function ci_method="owa".}
}

\references{

Yager, R. R. (1988). On ordered weighted averaging aggregation operators in multicriteria decision making. IEEE Transactions on systems, Man, and Cybernetics, 18(1), 183-190.

}

\author{Fusco E., Liborio M.P.}

\seealso{ \code{\link{ci_ogwa}}}

\examples{
data(data_HPI)

data_HPI = data_HPI[complete.cases(data_HPI),]
data_HPI_2019 = data_HPI[data_HPI$year==2019,]

Indic_name = c("Life_Expectancy","Ladder_of_life","Ecological_Footprint")
Indic_norm = data.frame("ISO"=data_HPI_2019$ISO, 
                        normalise_ci(data_HPI_2019[, Indic_name], 
                        c(1:3), 
                        c("POS","POS","NEG"),
                        method=2)$ci_norm)
                        
Indic_norm = Indic_norm[Indic_norm$Life_Expectancy>0 & 
                         Indic_norm$Ladder_of_life>0 & 
                         Indic_norm$Ecological_Footprint >0 ,]

atleast = 2
CI_owa_n = ci_owa(Indic_norm, id="ISO", 
                   indic_col=c(2:4), 
                   atleastjp=atleast)$CI_OWA_n
CI_owa_p = ci_owa(Indic_norm, id="ISO", 
                   indic_col=c(2:4), 
                   atleastjp=atleast)$CI_OWA_p
}
