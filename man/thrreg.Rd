% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/thrreg.R
\docType{package}
\name{thrreg}
\alias{thrreg}
\title{thrreg: Threshold Regression Model}
\description{
thrreg estimates the threshold regression model by Hansen (2000).
It is based on the
\href{https://www.ssc.wisc.edu/~bhansen/progs/ecnmt_00.html}{source code}
provided by Bruce E. Hansen.
}
\details{
thrreg provides three functions:
    \itemize{
    \item \code{\link{thr_test_hom}} tests for a threshold in a linear
          regression under homoskedasticity.
    \item \code{\link{thr_test_het}} tests for a threshold in a linear
          regression under heteroskedasticity.
    \item \code{\link{thr_est}} estimates the threshold and the regression
          parameters of the threshold model.
    }
}
\examples{
## See README.md on GitHub for a comprehensive example.

}
\references{
Hansen, B. E. (2000). Sample splitting and threshold estimation.
    \emph{Econometrica}, 68(3):575--603.
    \url{https://doi.org/10.1111/1468-0262.00124}.
    \url{https://www.ssc.wisc.edu/~bhansen/papers/ecnmt_00.pdf}.
}
\seealso{
Useful links:
    \itemize{
    \item \url{https://github.com/mlkremer/thrreg}
    \item Report bugs at \url{https://github.com/mlkremer/thrreg/issues}
    }
}
\author{
Marcel Kremer, \email{marcel.kremer@uni-due.de}

Bruce E. Hansen, \email{behansen@wisc.edu}
}
