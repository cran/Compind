\name{normalise_ci}
\alias{normalise_ci}

\title{
Normalisation and polarity functions
}

\description{This function lets to normalise simple indicators according to the polarity of each one.}

\usage{
normalise_ci(x, indic_col, polarity, method=1, z.mean=0, z.std=1, ties.method ="average")
}

\arguments{
  \item{x}{A data frame containing simple indicators.}
  \item{indic_col}{Simple indicators column number.}
  \item{method}{Normalisation methods: 
						\itemize{
								\item \bold{1} (default) = standardization or z-scores using the following formulation: \deqn{z_{ij}=z.mean \pm \frac{x_{ij}-M_{x_j}}{S_{x_j}}\cdot z.std} where \eqn{\pm} depends on \emph{polarity} parameter and \emph{z.mean} and \emph{z.std} represent the shifting parameters.
                
								\item \bold{2} = Min-max method using the following formulation: 
								
								if \emph{polarity}="POS": \deqn{\frac{x-min(x)}{max(x)-min(x)}} 
                if \emph{polarity}="NEG": \deqn{\frac{max(x)-x}{max(x)-min(x)}}
								
								\item \bold{3} = Ranking method. If \emph{polarity}="POS" ranking is increasing, while if \emph{polarity}="NEG" ranking is decreasing.
						}
				}
  \item{polarity}{Polarity vector: "POS" = positive, "NEG" = negative. The polarity of a individual indicator is the sign of the relationship between the indicator and the phenomenon to be measured (e.g., in a well-being index, "GDP per capita" has 'positive' polarity and "Unemployment rate" has 'negative' polarity).}
  \item{z.mean}{If method=1, Average shifting parameter. Default is 0.}
  \item{z.std}{If method=1, Standard deviation expansion parameter. Default is 1.}
  \item{ties.method}{If method=3, A character string specifying how ties are treated, see \code{\link{rank}} for details. Default is "average".}
}

%DA FINIRE _________________________________
%\details{ 
%
%   Min-max method in case of NEG polarity
%   
%   Ranking method in case of NEG polarity  
%______________________________________________________________
%}

\value{
  \item{ci_norm}{A data.frame containing normalised score of the choosen simple indicators.}
  \item{norm_method}{Normalisation method used.}
}

\references{OECD, "\emph{Handbook on constructing composite indicators: methodology and user guide}", 2008, pag.30.}

\author{ Vidoli F.}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_mpi}}}

\examples{
data(EU_NUTS1)

# Standard z-scores normalisation #
data_norm = normalise_ci(EU_NUTS1,c(2:3),c("NEG","POS"),method=1,z.mean=0, z.std=1)
summary(data_norm$ci_norm)

# Normalisation for MPI index #
data_norm = normalise_ci(EU_NUTS1,c(2:3),c("NEG","POS"),method=1,z.mean=100, z.std=10)
summary(data_norm$ci_norm)

data_norm = normalise_ci(EU_NUTS1,c(2:3),c("NEG","POS"),method=2)
summary(data_norm$ci_norm)
}
