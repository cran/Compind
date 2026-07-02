\name{ci_bod_time_decomp}
\alias{ci_bod_time_decomp}

\title{Dynamic Benefit of the Doubt decomposition for panel data}

\description{
Computes dynamic Benefit of the Doubt (BoD) composite indicator scores for
panel data using contemporaneous, sequential and intertemporal reference
technologies. The function also computes Total Change, Catch-Up and Benchmark
Shift indices.
}

\usage{
ci_bod_time_decomp(
  data_long,
  id_col,
  time_col,
  indic_cols
)
}

\arguments{

  \item{data_long}{
  A data frame in long format containing the panel data.
  }

  \item{id_col}{
  Character string identifying the column containing the unit identifiers.
  }

  \item{time_col}{
  Character string identifying the column containing the time variable.
  }

  \item{indic_cols}{
  Character vector specifying the indicator variables used to compute the
  Benefit of the Doubt composite indicator.
  }
}

\details{

For each period the function computes three BoD efficiency scores:

\itemize{
  \item contemporaneous scores (\code{E_C});
  \item sequential scores (\code{E_S});
  \item intertemporal scores (\code{E_I}).
}

These scores are then used to compute:

\itemize{
  \item Total Change (\code{TotalChange});
  \item Catch-Up (\code{CatchUp});
  \item Benchmark Shift (\code{BenchmarkShift});
}
}

\value{

A data frame containing:

\itemize{
  \item \code{id}: unit identifier;
  \item \code{time}: time period;
  \item \code{E_C}: contemporaneous BoD score;
  \item \code{E_S}: sequential BoD score;
  \item \code{E_I}: intertemporal BoD score;
  \item \code{TotalChange}: Total Change index;
  \item \code{CatchUp}: Catch-Up index;
  \item \code{BenchmarkShift}: Benchmark Shift index.
}
}

\references{

Fusco, E. and Magrini, A. (2026) 
\emph{Dynamic Benefit of the Doubt Decomposition for Panel Data: Evidence from Sustainable Energy in the EU}.
Sustainability, 18, 5835.
}

\author{
E.Fusco and A.Magrini
}

\seealso{
\code{\link{ci_bod}}
}

\examples{

panel <- data.frame(
  id = rep(1:3, each = 4),
  time = rep(1:4, times = 3),
  I1 = c(
    0.40, 0.43, 0.47, 0.52,
    0.60, 0.63, 0.67, 0.72,
    0.80, 0.84, 0.89, 0.95
  ),
  I2 = c(
    0.50, 0.54, 0.58, 0.63,
    0.70, 0.74, 0.79, 0.84,
    0.90, 0.95, 1.00, 1.06
  )
)

res <- ci_bod_time_decomp(
  data_long = panel,
  id_col = "id",
  time_col = "time",
  indic_cols = c("I1", "I2")
)

#head(res)

}