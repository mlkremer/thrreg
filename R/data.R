#' Durlauf and Johnson (1995) data
#'
#' Data distributed by Durlauf and Johnson (1995) to document their empirical work.
#'
#' @format A data frame with 96 rows and 7 variables:
#' \describe{
#'   \item{GDPGwth}{GDP growth as log difference between per capita GDP in 1985 and 1960.}
#'   \item{LogGDP1960}{Log per capita GDP in 1960.}
#'   \item{LogInvGDP}{Log average ratio of investment (including government investment)
#'     to GDP from 1960 to 1985.}
#'   \item{LogPopGwth}{Log average growth rate of working-age population from 1960 to 1985.}
#'   \item{LogSchool}{Log average fraction of working-age population enrolled in
#'     secondary school from 1960 to 1985.}
#'   \item{GDP1960}{Per capita GDP in 1960.}
#'   \item{Literacy}{Fraction of the population over 15 years old that is able
#'     to read and write in 1960.}
#' }
#'
#' @source
#'     Original data are reported in the appendix of Durlauf and Johnson (1995),
#'     pp. 379--380, and can be downloaded under
#'     \url{http://qed.econ.queensu.ca/jae/1995-v10.4/durlauf-johnson/}.
#'     Original data have been transformed to match Eq. (7) in Durlauf and
#'     Johnson (1995) and Section 5 in Hansen (2000), respectively.
#'     Documentation is taken from the Journal of Applied Econometrics FTP
#'     (see link above).
#'     All of the data with the exception of Literacy are from Mankiw, Romer and
#'     Weil (1992). Literacy is from the World Bank's World Development Report.
#'
#' @references
#'     Durlauf, S. N. and Johnson, P. A. (1995).
#'     Multiple regimes and cross-country growth behavior.
#'     \emph{Journal of Applied Econometrics}, 10(4):365--384.
#'     \url{https://doi.org/10.1002/jae.3950100404}.
#' @references
#'    Mankiw, N. G., Romer, D., and Weil, D. N. (1992).
#'    A contribution to the empirics of economic growth.
#'    \emph{The Quarterly Journal of Economics}, 107(2):407--437.
#'    \url{https://doi.org/10.2307/2118477}.
#'
"dur_john"

#' Original Durlauf and Johnson (1995) data
#'
#' Data distributed by Durlauf and Johnson (1995) to document their empirical work.
#'
#' @format A data frame with 121 rows and 11 variables:
#' \describe{
#'   \item{NUMBER}{Country number in Summers and Heston (1988) dataset.}
#'   \item{NONOIL}{dummy; \code{= 1} for nonoil producing countries.}
#'   \item{INTER}{dummy; \code{= 1} for countries with better quality data.}
#'   \item{OECD}{dummy; \code{= 1} for OECD countries.}
#'   \item{GDP60}{Per capita GDP in 1960.}
#'   \item{GDP85}{Per capita GDP in 1985.}
#'   \item{GDPGRO}{Average growth rate of per capita GDP from 1960 to 1985.}
#'   \item{POPGRO}{Average growth rate of working-age population from 1960 to 1985.}
#'   \item{IONY}{Average ratio of investment (including government investment)
#'     to GDP from 1960 to 1985.}
#'   \item{SCHOOL}{Average fraction of working-age population enrolled in
#'     secondary school from 1960 to 1985.}
#'   \item{LIT60}{Fraction of the population over 15 years old that is able
#'     to read and write in 1960.}
#' }
#'
#' @source
#'     \url{http://qed.econ.queensu.ca/jae/1995-v10.4/durlauf-johnson/}.
#'
#'     Original data are reported in the appendix of Durlauf and Johnson (1995),
#'     pp. 379--380. Documentation is taken from the Journal of Applied
#'     Econometrics FTP (see link above).
#'     All of the data with the exception of LIT60 are from Mankiw, Romer
#'     and Weil (1992). LIT60 is from the World Bank's World Development Report.
#'
#' @note A value of -999 indicates that the observation is missing.
