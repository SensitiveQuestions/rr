#' R Package for the Randomized Response Technique
#' 
#' \code{rr} implements methods developed by Blair, Imai, and Zhou (2015) such
#' as multivariate regression and power analysis for the randomized response
#' technique. Randomized response is a survey technique that introduces random
#' noise to reduce potential bias from non-response and social desirability
#' when asking questions about sensitive behaviors and beliefs. The current
#' version of this package conducts multivariate regression analyses for the
#' sensitive item under four standard randomized response designs: mirrored
#' question, forced response, disguised response, and unrelated question.
#' Second, it generates predicted probabilities of answering affirmatively to
#' the sensitive item for each respondent. Third, it also allows users to use
#' the sensitive item as a predictor in an outcome regression under the forced
#' response design. Additionally, it implements power analyses to help improve
#' research design. In future versions, this package will extend to new
#' modified designs that are based on less stringent assumptions than those of
#' the standard designs, specifically to allow for non-compliance and unknown
#' distribution to the unrelated question under the unrelated question design.
#' 
#' \tabular{ll}{ Package: \tab rr\cr Type: \tab Package\cr Version: \tab 1.2\cr
#' Date: \tab 2015-3-8\cr License: \tab GPL (>= 2)\cr }
#' 
#' @name rr-package
#' @useDynLib rr
#' @aliases rr-package rr
#' @docType package
#' @author Graeme Blair, Experiments in Governance and Politics, Columbia
#' University \email{graeme.blair@@gmail.com}, \url{https://graemeblair.com}
#' 
#' Kosuke Imai, Departments of Government and Statistics, Harvard University
#' \email{kimai@@harvard.edu}, \url{https://imai.fas.harvard.edu}
#' 
#' Yang-Yang Zhou, Department of Political Science, University of British Columbia
#' \email{yangyang.zhou@@ubc.ca}, \url{https://www.yangyangzhou.com}
#' 
#' Maintainer: Graeme Blair <graeme.blair@@gmail.com>
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2015) "Design
#' and Analysis of the Randomized Response Technique."  
#' \emph{Journal of the American Statistical Association.}
#' Available at \url{https://graemeblair.com/papers/randresp.pdf}.
#' @keywords package
#' @importFrom graphics abline axis lines mtext plot
#' @importFrom stats as.formula binomial coef complete.cases cov dcauchy df glm glm.control model.frame model.matrix model.matrix.default model.response na.omit pnorm qnorm quantile rnorm runif sd vcov
NULL



