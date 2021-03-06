% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thr_test_hom.R
\name{thr_test_hom}
\alias{thr_test_hom}
\title{Threshold Test under Homoskedasticity}
\usage{
thr_test_hom(
  df,
  yi,
  xi,
  qi,
  var.names = colnames(df),
  trim_per = 0.15,
  rep = 1000,
  cr = 0.95,
  graph = TRUE,
  quick = 1
)
}
\arguments{
\item{df}{Data frame.}

\item{yi}{Integer or character; index or column name of dependent (y)
variable in \code{df}.}

\item{xi}{Integer or character vector; indexes or column names of
independent (x) variables in \code{df}.}

\item{qi}{Integer or character; index or column name of threshold (q)
variable in \code{df}.}

\item{var.names}{Character vector; variable names with
\code{length(var.names) == ncol(df)}
corresponding to columns in \code{df} to be used in threshold regression
table. Default is \code{colnames(df)}.}

\item{trim_per}{Numeric; percentage of sample to trim from ends.
Default is \code{trim_per = .15}.}

\item{rep}{Integer; number of bootstrap replications. Default is
\code{rep = 1000}.}

\item{cr}{Numeric; confidence level used to plot the critical value in the
graph. It is not used elsewhere in the analysis. Default is
\code{cr = .95}.}

\item{graph}{Logical; graph indicator.
Set \code{TRUE} (default) to view the graph of the likelihood;
set \code{FALSE} otherwise.}

\item{quick}{Integer; indicator of method used for bootstrap.
Set \code{quick = 1} (default) for efficient, quick computation;
set \code{quick = 0} if memory is limited (perhaps for large data sets).}
}
\value{
A list with components:
    \describe{
    \item{\code{f_test}}{the value of Maximal (Quandt) F-statistic.}
    \item{\code{p_value}}{the bootstrap p-value.}
    }
}
\description{
Computes a test for a threshold in linear regression under homoskedasticity.
}
\details{
\enumerate{
    \item Do not include a constant in the independent variables;
          the function automatically adds an intercept to the regression.
    \item The function stores the sequential design matrices in order to
          speed up the bootstrap computations. It is possible that if your
          dataset is very large, this will exceed your computer RAM memory.
          If so, function will crash, and the message
          \code{Error: allocMatrix: too many elements specified}
          will be displayed. If more RAM is not available, switch to
          \code{quick = 0}. The switch \code{quick = 0} requires the bootstrap
          to re-calculate the design matrices for each bootstrap replication,
          which requires less memory, but somewhat more computer time.
    }
}
\examples{
\donttest{
## Performs part of the empirical work reported in Hansen (2000)
data <- dur_john
output <- thr_test_hom(data, 1, 2:5, 6)

output$f_test
output$p_value
}


}
\references{
Hansen, B. E. (2000). Sample splitting and threshold estimation.
    \emph{Econometrica}, 68(3):575--603.
    \url{https://doi.org/10.1111/1468-0262.00124}.
    \url{https://www.ssc.wisc.edu/~bhansen/papers/ecnmt_00.pdf}.
}
\seealso{
\code{\link{thr_test_het}} for threshold test under heteroskedasticity,
    \code{\link{thr_est}} for threshold estimation.
}
\author{
Marcel Kremer, \email{marcel.kremer@uni-due.de}

Bruce E. Hansen, \email{behansen@wisc.edu}
}
\keyword{htest}
\keyword{models}
\keyword{ts}
