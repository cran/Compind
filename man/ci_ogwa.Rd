\name{ci_ogwa}
\alias{ci_ogwa}

\title{Ordered Geographically Weighted Average (OWA)}
\description{The Ordered Geographically Weighted Averaging (OWA) operator is an extension of the multi-criteria decision aggregation method called OWA (Yager, 1988) that accounts for spatial heterogeneity.}

\usage{ci_ogwa(x, id, indic_col, atleastjp, coords,
                    kernel = "bisquare", adaptive = F, bw, 
                    p = 2, theta = 0, longlat = F, dMat)}

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
  \item{id}{Units' unique identifier.}
  \item{indic_col}{Simple indicators column number.}
  \item{coords}{A two-column matrix of latitude and longitude coordinates.}
  \item{atleastjp}{Fuzzy linguistic quantifier "At least j".}
  \item{kernel}{function chosen as follows:
                gaussian: wgt = exp(-.5*(vdist/bw)^2);
                exponential: wgt = exp(-vdist/bw);
                bisquare: wgt = (1-(vdist/bw)^2)^2 if vdist < bw, wgt=0 otherwise;
                tricube: wgt = (1-(vdist/bw)^3)^3 if vdist < bw, wgt=0 otherwise;
                boxcar: wgt=1 if dist < bw, wgt=0 otherwise.}
  \item{adaptive}{if TRUE calculate an adaptive kernel where the bandwidth (bw) corresponds to                      the number of nearest neighbours (i.e. adaptive distance); default is FALSE,                      where a fixed kernel is found (bandwidth is a fixed distance).}
  \item{bw}{bandwidth used in the weighting function.}
  \item{p}{the power of the Minkowski distance, default is 2, i.e. the Euclidean distance.}	
  \item{theta}{an angle in radians to rotate the coordinate system, default is 0.}	
  \item{longlat}{if TRUE, great circle distances will be calculated.}	
  \item{dMat}{a pre-specified distance matrix, it can be calculated by the function \code{\link[GWmodel]{gw.dist}}.}	

}

\value{An object of class "CI". This is a list containing the following elements:
  \item{CI_OGWA_n}{Composite indicator estimated values for OGWA-.}
  \item{CI_OGWA_p}{Composite indicator estimated values for OGWA+.}
  \item{wp}{OGWA weights' vector "More than j".}
  \item{wn}{OGWA weights' vector "At least j".}
  \item{ci_method}{Method used; for this function ci_method="ogwa".}
}

\references{
Fusco, E., Liborio, M.P., Rabiei-Dastjerdi, H., Vidoli, F., Brunsdon, C. and Ekel, P.I. (2023), Harnessing Spatial Heterogeneity in Composite Indicators through the Ordered Geographically Weighted Averaging (OGWA) Operator. Geographical Analysis. https://doi.org/10.1111/gean.12384
}

\author{Fusco E., Liborio M.P.}

\seealso{ \code{\link{ci_owa}}}

\examples{
data(data_HPI)

data_HPI_2019 = data_HPI[data_HPI$year==2019,]
Indic_name = c("Life_Expectancy","Ladder_of_life","Ecological_Footprint")

Indic_norm = normalise_ci(data_HPI_2019, Indic_name, c("POS","POS","NEG"),method=2)$ci_norm
Indic_norm = Indic_norm[Indic_norm$Life_Expectancy>0 & 
              Indic_norm$Ladder_of_life>0 & 
              Indic_norm$Ecological_Footprint >0,]

Indic_CI = data.frame(Indic_norm, 
                       data_HPI_2019[rownames(Indic_norm),
                       c("lat","long","HPI","ISO","Country")])

atleast = 2

coord = Indic_CI[,c("lat","long")]

CI_ogwa_n = ci_ogwa(Indic_CI, id="ISO", 
                      indic_col=c(1:3), 
                      atleastjp=atleast,
                      coords=as.matrix(coord),
                      kernel = "gaussian",
                      adaptive=FALSE, 
                      longlat=FALSE)$CI_OGWA_n


#CI_ogwa_p = ci_ogwa(Indic_CI, id="ISO", 
#                      indic_col=c(1:3), 
#                     atleastjp=atleast,
#                      coords=as.matrix(coord),
#                      kernel = "gaussian",
#                      adaptive=FALSE, 
#                      longlat=FALSE)$CI_OGWA_p

}
