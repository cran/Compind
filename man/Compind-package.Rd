\name{Compind-package}
\alias{Compind-package}
%\alias{Compind}
\docType{package}

\title{Composite Indicators - Compind}

\description{Compind package contains functions to enhance several approaches to the Composite Indicators (CIs) methods, focusing, in particular, on the normalisation and weighting-aggregation steps.}

\details{
\tabular{ll}{
Package: \tab Compind\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-01-01\cr
Depends: \tab Benchmarking, Hmisc, MASS, ggplot2, psych, GPArotation,  lpSolve, nonparaeff, boot\cr
License: \tab GPL-3\cr
Built: \tab R 3.0.2; ; 2014-02-06 13:14:40 UTC; unix\cr
}

Index:
\preformatted{
Compind-package         Composite Indicators - Compind
EU_2020                 Europe 2020 indicators
EU_NUTS1                EU NUTS1 Transportation data
ci_bod                  Benefit of the Doubt approach (BoD)
ci_bod_dir              Directional Benefit of the Doubt (D-BoD) model
ci_bod_var_w            Variance weighted Benefit of the Doubt approach
                        (BoD variance weighted)
ci_factor               Weighting method based on Factor Analysis
ci_mean_geom            Weighting method based on geometric aggregation
ci_mpi                  Mazziotta-Pareto Index (MPI) method
ci_rbod                 Robust Benefit of the Doubt approach (RBoD)
ci_rbod_dir             Directional Robust Benefit of the Doubt
                        approach (D-RBoD)
ci_wroclaw              Wroclaw Taxonomic Method
normalise_ci            Normalisation and polarity functions
plot_influent           Plot influents units in terms of horizontal
                        mean and variability
}

}
\author{
Francesco Vidoli, Elisa Fusco
Maintainer: Francesco Vidoli <fvidoli@gmail.com>
}

\references{

Daraio, C., Simar, L. "\emph{Introducing environmental variables in nonparametric frontier models: a probabilistic approach}", Journal of productivity analysis, 2005, 24(1), 93?121.

Fusco E., "\emph{Enhancing non compensatory composite indicators: A directional proposal}", 2013, unpublished.

OECD, "\emph{Handbook on constructing composite indicators: methodology and user guide}", 2008.

Mazziotta C., Mazziotta M., Pareto A., Vidoli F., "\emph{La sintesi di indicatori territoriali di dotazione infrastrutturale: metodi di costruzione e procedure di ponderazione a confronto}", Rivista di Economia e Statistica del territorio, n.1, 2010.

Melyn W. and Moesen W.W., "\emph{Towards a synthetic indicator of macroeconomic performance: unequal weighting when limited information is available}", Public Economic research Paper 17, CES, KU Leuven, 1991.

Simar L., Vanhems A., "\emph{Probabilistic characterization of directional distances and their robust versions}", Journal of Econometrics, 2012, 166(2), 342?354.

UNESCO, "\emph{Social indicators: problems of definition and of selection}", Paris 1974.

Vidoli F., Fusco E., Mazziotta C., "\emph{Non-compensability in composite indicators: a robust directional frontier method}", Social Indicators Research, Springer Netherlands.

Vidoli F., Mazziotta C., "\emph{Robust weighted composite indicators by means of frontier methods
with an application to European infrastructure endowment}", Statistica Applicata, Italian
Journal of Applied Statistics, 2013.

}