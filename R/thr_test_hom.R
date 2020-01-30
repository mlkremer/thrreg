#' Threshold test under homoskedasticity
#'
#' Computes a test for a threshold in linear regression under homoskedasticity.
#'
#' @param df Data frame.
#' @param yi Integer or character; index or column name of dependent (y)
#'     variable in \code{df}.
#' @param xi Integer or character vector; indexes or column names of
#'     independent (x) variables in \code{df}.
#' @param qi Integer or character; index or column name of threshold (q)
#'     variable in \code{df}.
#' @param var.names Character vector; variable names with
#'     \code{length(var.names) == ncol(df)}
#'     corresponding to columns in \code{df} to be used in threshold regression
#'     table. Default is \code{colnames(df)}.
#' @param trim_per Numeric; percentage of sample to trim from ends.
#'     Default is \code{trim_per = .15}.
#' @param rep Integer; number of bootstrap replications. Default is \code{rep = 1000}.
#' @param cr Numeric; confidence level used to plot the critical value in the graph.
#'     It is not used elsewhere in the analysis. Default is \code{cr = .95}.
#' @param graph Logical; graph indicator.
#'     Set \code{TRUE} (default) to view the graph of the likelihood;
#'     set \code{FALSE} otherwise.
#' @param quick Integer; indicator of method used for bootstrap.
#'     Set \code{quick = 1} (default) for efficient, quick computation;
#'     set \code{quick = 0} if memory is limited (perhaps for large data sets).
#'
#' @details
#'     \enumerate{
#'     \item Do not include a constant in the independent variables;
#'           the function automatically adds an intercept to the regression.
#'     \item The function stores the sequential design matrices in order to
#'           speed up the bootstrap computations. It is possible that if your
#'           dataset is very large, this will exceed your computer RAM memory.
#'           If so, function will crash, and the message
#'           \code{Error: allocMatrix: too many elements specified}
#'           will be displayed. If more RAM is not available, switch to
#'           \code{quick = 0}. The switch \code{quick = 0} requires the bootstrap
#'           to re-calculate the design matrices for each bootstrap replication,
#'           which requires less memory, but somewhat more computer time.
#'     }
#'
#' @return A list with components:
#'     \describe{
#'     \item{\code{f_test}}{the value of Maximal (Quandt) F-statistic.}
#'     \item{\code{p_value}}{the bootstrap p-value.}
#'     }
#'
#' @inherit thrreg author references
#'
# #' @family threshold regression functions   # Automatically links to other functions in package
#'
#' @seealso
#'     \code{\link{thr_test_het}} for threshold test under heteroskedasticity,
#'     \code{\link{thr_est}} for threshold estimation.
#'
#' @keywords htest models ts
#'
#' @examples
#' \donttest{
#' ## Performs part of the empirical work reported in Hansen (2000)
#' data <- dur_john
#' output <- thr_test_hom(data, 1, 2:5, 6)
#'
#' output$f_test
#' output$p_value
#' }
#'
# #' @rdname thr_test
#'
#' @importFrom graphics legend lines plot title
#' @importFrom stats rnorm
#'
#' @export
#'
thr_test_hom <- function(df, yi, xi, qi, var.names = colnames(df), trim_per = .15,
                         rep = 1000, cr = .95, graph = TRUE, quick = 1) {

  dat <- as.matrix(df)
  if (is.character(yi)) yi <- which(colnames(df) == yi)
  if (is.character(xi)) xi <- which(colnames(df) %in% xi)
  if (is.character(qi)) qi <- which(colnames(df) == qi)

  n <- nrow(dat)
  q <- dat[,qi]
  qs <- order(q)
  y <- as.matrix(dat[qs,yi])
  x <- cbind(matrix(c(1),n,1),dat[qs,xi])
  q <- as.matrix(q[qs])
  k <- ncol(x)
  qname <- var.names[qi]
  qs <- unique(q)
  qn <- length(qs)
  qq <- matrix(c(0),qn,1)
  for (r in 1:qn) qq[r] <- colSums(q==qs[r])
  cqq <- cumsum(qq)
  sq <- (cqq>=floor(n*trim_per))*(cqq<=(floor(n*(1-trim_per))))
  qs <- as.matrix(qs[sq>0])
  cqq <- as.matrix(cqq[sq>0])
  qn <- nrow(qs)

  mi <- solve(t(x)%*%x)
  e <- y-x%*%mi%*%(t(x)%*%y)
  ee <- t(e)%*%e
  xe <- x*(e%*%matrix(c(1),1,k))
  xe <- apply(xe,2,cumsum)
  sn <- matrix(c(0),qn,1)

  if (quick == 1){
    mmistore <- matrix(c(0),k*(k+1)/2,qn)
    cqqb <- 1
    mm <- matrix(c(0),k,k)
    for (r in 1:qn) {
      cqqr <- cqq[r]
      if (cqqb==cqqr) {
        mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
      } else { mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]}
      sume <- as.matrix(xe[cqqr,])
      mmi <- solve(mm - mm%*%mi%*%mm)
      sn[r] <- t(sume)%*%mmi%*%sume
      cqqb <- cqqr+1
      ii <- 1
      for (i in 1:k) {
        mmistore[ii:(ii+i-1),r] <- mmi[i,1:i]
        ii <- ii+i
      }
    }
    si <- which.max(sn)
    smax <- sn[si]
    qmax <- qs[si]
    sig <- (ee - smax)/n
    lr <- sn/as.vector(sig)
    ftest <- smax/sig
    fboot <- matrix(c(0),rep,1)
    for (j in 1:rep) {
      y  <- rnorm(n)
      e  <- y-x%*%mi%*%(t(x)%*%y)
      ee <- t(e)%*%e
      xe <- x*(e%*%matrix(c(1),1,k))
      xe <- apply(xe,2,cumsum)
      sn <- matrix(c(0),qn,1)
      for (r in 1:qn) {
        mmi <- matrix(c(0),k,k)
        ii <- 1
        for (i in 1:k) {
          mmi[i,1:i] <- mmistore[ii:(ii+i-1),r]
          mmi[1:(i-1),i] <- mmi[i,1:(i-1)]
          ii <- ii+i
        }
        sume <- as.matrix(xe[cqq[r],])
        sn[r] <- t(sume)%*%mmi%*%sume
      }
      smax <- max(sn)
      sig <- (ee - smax)/n
      fboot[j] <- smax/sig
    }
  }

  if (quick == 0) {
    cqqb <- 1
    mm <- matrix(c(0),k,k)
    for (r in 1:qn) {
      cqqr <- cqq[r]
      if (cqqb==cqqr) {
        mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
      } else { mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]}
      sume <- as.matrix(xe[cqqr,])
      mmi <- mm - mm%*%mi%*%mm
      sn[r] <- t(sume)%*%solve(mmi)%*%sume
      cqqb <- cqqr+1
    }
    si <- which.max(sn)
    smax <- sn[si]
    qmax <- qs[si]
    sig <- (ee - smax)/n
    lr <- sn/as.vector(sig)
    ftest <- smax/sig
    fboot <- matrix(c(0),rep,1)
    for (j in 1:rep) {
      y  <- rnorm(n)
      e  <- y-x%*%mi%*%(t(x)%*%y)
      ee <- t(e)%*%e
      xe <- x*(e%*%matrix(c(1),1,k))
      xe <- apply(xe,2,cumsum)
      sn <- matrix(c(0),qn,1)
      cqqb <- 1
      mm <- matrix(c(0),k,k)
      for (r in 1:qn) {
        cqqr <- cqq[r]
        if (cqqb==cqqr) {
          mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
        } else { mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]}
        mmi <- mm - mm%*%mi%*%mm
        sume <- as.matrix(xe[cqq[r],])
        sn[r] <- t(sume)%*%solve(mmi)%*%sume
        cqqb <- cqqr+1
      }
      smax <- max(sn)
      sig <- (ee - smax)/n
      fboot[j] <- smax/sig
    }
  }

  fboot <- as.matrix(sort(fboot))
  pv <- mean(fboot >= matrix(c(1),rep,1)%*%ftest)
  crboot <- fboot[round(rep*cr)]

  if (graph==TRUE) {
    xxlim <- range(qs)
    yylim <- range(rbind(lr,crboot))
    clr <- matrix(c(1),qn,1)*crboot
    plot(qs,lr,lty=1,col=1,xlim=xxlim,ylim=yylim,type="l",ann=0)
    lines(qs,clr,lty=2,col=2)
    title(main=rbind("F Test For Threshold",
                     "Reject Linearity if F Sequence Exceeds Critical Value"),
          xlab="gamma",ylab="Fn(gamma)")
    tit <- paste(cr*100,c("% Critical"),sep="")
    legend("bottomright",c("LRn(gamma)",tit),lty=c(1,2),col=c(1,2))
  }

  cat(paste0("\\subsection{Testing for a Sample Split, Using ", qname, "}"), "\n")
  cat("\\subsubsection*{Test of Null of No Threshold Against Alternative of Threshold}", "\n")
  cat("Under Maintained Assumption of Homoskedastic Errors", "\\\\")
  cat("\\\\\n")
  cat("Number of Bootstrap Replications: ", rep, "\\\\\n")
  cat("Trimming Percentage:              ", trim_per, "\\\\")
  cat("\\\\\n")
  cat("Threshold Variable:                ", qname, "\\\\\n")
  cat("Threshold Estimate:               ", qmax, "\\\\\n")
  cat("F-test for no threshold:          ", ftest, "\\\\\n")
  cat("Bootstrap P-Value:                ", pv, "\\\\\n")
  cat("\n")

  list(f_test = as.vector(ftest), p_value = pv)
}
