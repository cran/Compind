\name{ci_bod_mdir}
\alias{ci_bod_mdir}
\title{Multi-directional Benefit of the Doubt approach (MDBoD)}
\description{Multi-directional Benefit of the Doubt (MDBoD) allows to introduce the non-compensability among simple indicators in a standard BOD in an objective manner: the preference structure, i.e., the direction, is determined directly from the data and is specific for each unit.}

\usage{ci_bod_mdir(x,indic_col)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_bod_mdir_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="bod".}
  \item{ci_bod_mdir_spec}{Simple indicators specific scores.}
  \item{ci_bod_mdir_dir}{Directions for each simple indicator and unit.}
  
}

\references{

Fusco E. (2023) "\emph{Potential improvements approach in composite indicators construction: the Multi-directional Benefit of the Doubt model}", Socio-Economic Planning Sciences, vol. 85, 101447
}

\author{Fusco E.}

\seealso{ \code{\link{ci_bod_dir}}, \code{\link{mea}}}

\examples{
data(BLI_2017)
CI <- ci_bod_mdir(BLI_2017,c(2:12))
}
