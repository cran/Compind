\name{ci_bod_dir}
\alias{ci_bod_dir}
\title{Directional Benefit of the Doubt (D-BoD) model}
\description{Directional Benefit of the Doubt (D-BoD) model enhance non-compensatory property by introducing directional penalties in a standard BoD model in order to consider the preference structure among simple indicators.}

\usage{ci_bod_dir(x, indic_col, dir)}

\arguments{
  \item{x}{A data.frame containing score of the simple indicators.}
   \item{indic_col}{Simple indicators column number.}
  \item{dir}{Main direction. For example you can set the average rates of substitution.}  
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_bod_dir_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="bod_dir".}
}


%\details{ 
%   INSERIRE FUNZIONE CHE CALCOLA TRAMITE ACP
%}

\references{Fusco E. (2015) "\emph{Enhancing non compensatory composite indicators: A directional proposal}", European Journal of Operational Research, 242(2), 620-630.}

\author{Vidoli F., Fusco E.}

\seealso{ \code{\link{ci_bod}}, \code{\link{ci_rbod}}}

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_bod_dir(Indic,dir=c(1,1))

data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,c(2:3),polarity = c("POS","POS"), method=2)
CI = ci_bod_dir(data_norm$ci_norm,c(1:2),dir=c(1,0.5))
}
