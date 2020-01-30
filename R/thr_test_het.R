#' Threshold test under heteroskedasticity
#'
#' Computes a test for a threshold in linear regression under heteroskedasticity.
#'
#' @param quick Integer; indicator of method used for bootstrap.
#'     Set \code{quick = 1} for quick computation of asymptotic distribution.
#'     This method is not a proper bootstrap and may result in excess rejections.
#'     It also uses more memory.
#'     Set \code{quick = 2} (default) for a better bootstrap procedure,
#'     which also uses less memory, but is more time consuming.
#'
#' @inheritParams thr_test_hom
#'
#' @details
#'     \enumerate{
#'     \item Do not include a constant in the independent variables;
#'           the function automatically adds an intercept to the regression.
#'     \item There are two bootstrap methods which the function can use.
#'
#'           The first method, obtained by setting \code{quick = 1},
#'           is the method presented in the paper
#'           Hansen, B. E. (1996). Inference When a Nuisance Parameter is Not
#'           Identified Under the Null Hypothesis.
#'           \emph{Econometrica}, 64(2):413-430.
#'           \url{https://www.ssc.wisc.edu/~bhansen/papers/ecnmt_96.pdf},
#'           which simulates the asymptotic null distribution.
#'           A computational shortcut is also taken which speeds computational
#'           time, at the cost of greater memory usage, so may not be possible
#'           for large data sets.
#'
#'           The second method, obtained by setting \code{quick = 2},
#'           is a "fixed regressor bootstrap", which is quite close.
#'           The difference is that the bootstrap procedure calculates
#'           the variance-covariance matrix in each bootstrap replication.
#'           This results in a better finite sample approximation.
#'           The cost is greater computation time.
#'
#'           The function is set by default to use the second method
#'           (\code{quick = 2}), which has better sampling properties.
#'           If computational time is a concern, switch to the first method
#'           (\code{quick = 1}). If an "out of workspace memory" message appears,
#'           switch back to \code{quick = 2}.
#'     }
#'
#' @inherit thr_test_hom return
#' @inherit thrreg author references
#'
# #' @family threshold regression functions
#'
#' @seealso
#'     \code{\link{thr_test_hom}} for threshold test under homoskedasticity,
#'     \code{\link{thr_est}} for threshold estimation.
#'
#' @keywords htest models ts
#'
#' @examples
#' \donttest{
#' ## Performs part of the empirical work reported in Hansen (2000)
#' data <- dur_john
#' output <- thr_test_het(data, 1, 2:5, 6)
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
thr_test_het <- function(df, yi, xi, qi, var.names = colnames(df), trim_per = .15,
                         rep = 1000, cr = .95, graph = TRUE, quick = 2) {

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

  mi <- solve(t(x)%*%x, tol = 1e-1000)
  e <- y-x%*%mi%*%(t(x)%*%y)
  ee <- t(e)%*%e
  xe <- x*(e%*%matrix(c(1),1,k))
  vi <- t(xe)%*%xe
  cxe <- apply(xe,2,cumsum)
  sn <- matrix(c(0),qn,1)

  if (quick == 1) {
    mmistore <- matrix(c(0),k*(k+1)/2,qn)
    cqqb <- 1
    mm <- matrix(c(0),k,k)
    vv <- matrix(c(0),k,k)
    for (r in 1:qn) {
      cqqr <- cqq[r]
      if (cqqb==cqqr) {
        mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
        vv <- vv + as.matrix(xe[cqqb,])%*%xe[cqqb,]
      } else {
        mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]
        vv <- vv + t(xe[(cqqb:cqqr),])%*%xe[(cqqb:cqqr),]
      }
      sume <- as.matrix(cxe[cqqr,])
      mmi <- solve(vv-mm%*%mi%*%vv-vv%*%mi%*%mm+mm%*%mi%*%vi%*%mi%*%mm)
      sn[r] <- t(sume)%*%mmi%*%sume
      cqqb <- cqqr+1
      ii <- 1
      for (i in 1:k) {
        mmistore[ii:(ii+i-1),r] <- mmi[i,1:i]
        ii <- ii+i
      }
    }
    si <- which.max(sn)
    qmax <- qs[si]
    lr <- sn
    ftest <- sn[si]
    fboot <- matrix(c(0),rep,1)
    for (j in 1:rep) {
      y  <- rnorm(n)*e
      xe <- x*((y-x%*%mi%*%(t(x)%*%y))%*%matrix(c(1),1,k))
      cxe <- apply(xe,2,cumsum)
      sn <- matrix(c(0),qn,1)
      for (r in 1:qn) {
        mmi <- matrix(c(0),k,k)
        ii <- 1
        for (i in 1:k) {
          mmi[i,1:i] <- mmistore[ii:(ii+i-1),r]
          mmi[1:(i-1),i] <- mmi[i,1:(i-1)]
          ii <- ii+i
        }
        sume <- as.matrix(cxe[cqq[r],])
        sn[r] <- t(sume)%*%mmi%*%sume
      }
      fboot[j] <- max(sn)
    }
  }

  if (quick == 2) {
    cqqb <- 1
    mm <- matrix(c(0),k,k)
    vv <- matrix(c(0),k,k)
    for (r in 1:qn) {
      cqqr <- cqq[r]
      if (cqqb==cqqr) {
        mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
        vv <- vv + as.matrix(xe[cqqb,])%*%xe[cqqb,]
      } else {
        mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]
        vv <- vv + t(xe[(cqqb:cqqr),])%*%xe[(cqqb:cqqr),]
      }
      sume <- as.matrix(cxe[cqqr,])
      mmi <- vv-mm%*%mi%*%vv-vv%*%mi%*%mm+mm%*%mi%*%vi%*%mi%*%mm
      if (qr(mmi)$rank==ncol(mmi)){
        sn[r] <- t(sume)%*%solve(mmi)%*%sume
      }
      cqqb <- cqqr+1
    }
    si <- which.max(sn)
    qmax <- qs[si]
    lr <- sn
    ftest <- sn[si]
    fboot <- matrix(c(0),rep,1)
    for (j in 1:rep) {
      y  <- rnorm(n)*e
      xe <- x*((y-x%*%mi%*%(t(x)%*%y))%*%matrix(c(1),1,k))
      vi <- t(xe)%*%xe
      cxe <- apply(xe,2,cumsum)
      sn <- matrix(c(0),qn,1)
      cqqb <- 1
      mm <- matrix(c(0),k,k)
      vv <- matrix(c(0),k,k)
      for (r in 1:qn) {
        cqqr <- cqq[r]
        if (cqqb==cqqr) {
          mm <- mm + as.matrix(x[cqqb,])%*%x[cqqb,]
          vv <- vv + as.matrix(xe[cqqb,])%*%xe[cqqb,]
        } else {
          mm <- mm + t(x[(cqqb:cqqr),])%*%x[(cqqb:cqqr),]
          vv <- vv + t(xe[(cqqb:cqqr),])%*%xe[(cqqb:cqqr),]
        }
        mmi <- vv-mm%*%mi%*%vv-vv%*%mi%*%mm+mm%*%mi%*%vi%*%mi%*%mm
        sume <- as.matrix(cxe[cqqr,])
        if (qr(mmi)$rank==ncol(mmi)) {
          sn[r] <- t(sume)%*%solve(mmi)%*%sume
        }
        cqqb <- cqqr+1
      }
      fboot[j] <- max(sn)
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
  cat("Allowing Heteroskedastic Errors (White Corrected)", "\\\\")
  cat("\\\\\n")
  cat("Number of Bootstrap Replications: ", rep, "\\\\\n")
  cat("Trimming Percentage:              ", trim_per, "\\\\")
  cat("\\\\\n")
  cat("Threshold Variable:                ", qname, "\\\\\n")
  cat("Threshold Estimate:               ", qmax, "\\\\\n")
  cat("LM-test for no threshold:         ", ftest, "\\\\\n")
  cat("Bootstrap P-Value:                ", pv, "\\\\\n")
  cat("\n")

  list(f_test = ftest, p_value = pv)
}
