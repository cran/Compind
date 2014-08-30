\name{ci_bod}
\alias{ci_bod}
\title{Benefit of the Doubt approach (BoD)}
\description{Benefit of the Doubt approach (BoD) is the application of Data Envelopment Analysis (DEA) to the field of composite indicators. It was originally proposed by Melyn and Moesen (1991) to evaluate macroeconomic performance.}
%"\emph{The Benefit of the Doubt approach is formally tantamount to the original input-oriented Data Envelopment Analysis with Constant Returns to Scale (CRS-DEA) model of Charnes et al. (1978), with all questionnaire items considered as outputs and a dummy input equal to one for all observations}", Witte and Rogge (2009).

\usage{ci_bod(x,indic_col)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
}

%\details{Please pay attention that \emph{"weights are not uniquely determined. Note also that the multiplicity of solutions is likely to depend %upon the set of constraints imposed on the weights of the maximisation problem, the wider the range of variation of weights, the lower the %possibility of obtaining a unique solution"}, OECD (2008).}

\value{An object of class "CI". This is a list containing the following elements:
  \item{ci_bod_est}{Composite indicator estimated values.}
  \item{ci_method}{Method used; for this function ci_method="bod".}
}

\references{

OECD, \emph{Handbook on constructing composite indicators: methodology and user guide}, 2008.

Melyn W. and Moesen W.W., "\emph{Towards a synthetic indicator of macroeconomic performance: unequal weighting when limited information is available}", Public Economic research Paper 17, CES, KU Leuven, 1991.

Witte, K. D., Rogge, N. "\emph{Accounting for exogenous influences in a benevolent performance evaluation of teachers}". Tech. rept. Working
Paper Series ces0913, Katholieke Universiteit Leuven, Centrum voor Economische Studien, 2009.
}

\author{Vidoli F.}

\seealso{ \code{\link{ci_bod_dir}},\code{\link{ci_rbod}} }

\examples{
i1 <- seq(0.3, 0.5, len = 100) - rnorm (100, 0.2, 0.03)
i2 <- seq(0.3, 1, len = 100)   - rnorm (100, 0.2, 0.03)
Indic = data.frame(i1, i2)
CI = ci_bod(Indic)

data(EU_NUTS1)
data_norm = normalise_ci(EU_NUTS1,c(2:3),polarity = c("POS","POS"), method=2)
CI = ci_bod(data_norm$ci_norm,c(1:2))
}
