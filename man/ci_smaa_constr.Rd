\name{ci_smaa_constr}
\alias{ci_smaa_constr}
\title{Constrained stochastic multi-objective acceptability analysis (C-SMAA)}
\description{Stochastic multiobjective acceptability analysis (SMAA) is a multicriteria decision support technique for multiple decision makers based on exploring the weight space. Inaccurate or uncertain input data can be represented as probability distributions. In SMAA the decision makers need not express their preferences explicitly or implicitly; instead the technique analyses what kind of valuations would make each alternative the preferred one. The method produces for each alternative an acceptability index measuring the variety of different valuations that support that alternative, a central weight vector representing the typical valuations resulting in that decision, and a confidence factor measuring whether the input data is accurate enough for making an informed decision. (R Lahdelma, J. Hokkanen and P. Salminen, 1998); this function, in particular, allows to restricts the range of allowable weights within the SMAA analysis.}

\usage{ci_smaa_constr(x,indic_col,rep, label, low_w=NULL)}

\arguments{
  \item{x}{A data.frame containing simple indicators.}
  \item{indic_col}{A numeric list indicating the positions of the simple indicators.}
  \item{rep}{Number of samples.}
  \item{label}{A factor column useful to identify units.}
  \item{low_w}{Importance weights lower bound vector; default is NULL (for standard SMAA)}
}

\details{Author thanks Giuliano Resce and Raffaele Lagravinese for their help and for making available the original code of the SMAA function.\\
The lower bound vector must be set as a vector of the same size as the number of simple indicators; for example - in the presence of two indicators - if you want to constrain only one indicator, you must write: \code{low_w = c (0,0.2)}.
}


\value{An object of class "CI". This is a list containing the following elements:
      \item{ci_smaa_constr_rank_freq}{Frequence of the SMAA ranks based on the sampled alternatives' values. The rows represent the analysis units while the first column represents the number of times the unit was in first rank, the second one in second rank and so on.}
          \item{ci_smaa_constr_average_rank}{The average rank.}
  \item{ci_smaa_constr_values}{The alternative values based on a set of samples from the criteria values distribution and the samples set from the feasible weight space.}
  \item{ci_method}{Method used; for this function ci_method="smaa_const".}

}

\references{

R. Lahdelma, P. Salminen (2001) "\emph{SMAA-2: Stochastic multicriteria acceptability analysis for group decision making}", Operations Research, 49(3), pp. 444-454

S. Greco, A. Ishizaka, B. Matarazzo and G. Torrisi (2017) "\emph{Stochastic multi-attribute acceptability analysis (SMAA): an application to the ranking of Italian regions}", Regional Studies

R. Lagravinese, P. Liberati and G. Resce (2017) "\emph{Exploring health outcomes by stochastic multi-objective acceptability analysis: an application to Italian regions}", Working Papers. Collection B: Regional and sectoral economics, 1703, Universidade de Vigo, GEN - Governance and Economics research Network.

}

\author{Vidoli F.}

\seealso{ \code{\link{ci_bod}}}

\examples{


# ----- Define a function for plotting a matrix ----- #
myImagePlot <- function(x, ...){
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     title <-c()
  # check for additional function arguments
  if( length(list(...)) ){
    Lst <- list(...)
    if( !is.null(Lst$zlim) ){
       min <- Lst$zlim[1]
       max <- Lst$zlim[2]
    }
    if( !is.null(Lst$yLabels) ){
       yLabels <- c(Lst$yLabels)
    }
    if( !is.null(Lst$xLabels) ){
       xLabels <- c(Lst$xLabels)
    }
    if( !is.null(Lst$title) ){
       title <- Lst$title
    }
  }
# check for null values
if( is.null(xLabels) ){
   xLabels <- c(1:ncol(x))
}
if( is.null(yLabels) ){
   yLabels <- c(1:nrow(x))
}

layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))

 # Red and green range from 0 to 1 while Blue ranges from 1 to 0
 ColorRamp <- rgb( seq(0,1,length=256),  # Red
                   seq(0,1,length=256),  # Green
                   seq(1,0,length=256))  # Blue
 ColorLevels <- seq(min, max, length=length(ColorRamp))

 # Reverse Y axis
 reverse <- nrow(x) : 1
 yLabels <- yLabels[reverse]
 x <- x[reverse,]

 # Data Map
 par(mar = c(3,5,2.5,2))
 image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp, xlab="",
 ylab="", axes=FALSE, zlim=c(min,max))
 if( !is.null(title) ){
    title(main=title)
 }
axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
 axis(LEFT <-2, at=1:length(yLabels), labels=yLabels, las= HORIZONTAL<-1,
 cex.axis=0.7)

 # Color Scale
 par(mar = c(3,2.5,2.5,2))
 image(1, ColorLevels,
      matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1),
      col=ColorRamp,
      xlab="",ylab="",
      xaxt="n")

 layout(1)
}
# ----- END plot function ----- #


data(EU_NUTS1)

# Standard SMAA 
test <- ci_smaa_constr(EU_NUTS1,c(2,3), rep=200, label = EU_NUTS1[,1])
# source("http://www.phaget4.org/R/myImagePlot.R") 
# myImagePlot(test$ci_smaa_constr_rank_freq)
test$ci_smaa_constr_average_rank

# Constrained SMAA 
test2 <- ci_smaa_constr(EU_NUTS1,c(2,3), rep=200, label = EU_NUTS1[,1], low_w=c(0.2,0.2) )
# myImagePlot(test2$ci_smaa_constr_rank_freq)
test2$ci_smaa_constr_average_rank

}
