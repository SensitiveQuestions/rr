

#' Nigeria Randomized Response Survey Experiment on Social Connections to Armed
#' Groups
#' 
#' This data set is a subset of the data from the randomized response technique
#' survey experiment conducted in Nigeria to study civilian contact with armed
#' groups. The survey was implemented by Blair (2014).
#' 
#' 
#' @format A data frame containing 2457 observations.  The variables are:
#' \itemize{ \item \code{Quesid}: Survey ID of civilian respondent.  \item
#' \code{rr.q1}: Randomized response survey item using the Forced Response
#' Design asking the respondent whether they hold direct social connections
#' with members of armed groups. 0 if no connection; 1 if connection.  \item
#' \code{cov.age}: Age of the respondent.  \item \code{cov.asset.index}: The
#' number of assets owned by the respondent from an index of nine assets
#' including radio, T.V., motorbike, car, mobile phone, refrigerator, goat,
#' chicken, and cow.  \item \code{cov.married}: Marital status. 0 if single; 1
#' if married.  \item \code{cov.education}: Education level of the respondent.
#' 1 if no school; 2 if started primary school; 3 if finished primary school; 4
#' if started secondary school; 5 if finished secondary school; 6 if started
#' polytechnic or college; 7 if finished polytechnic or college; 8 if started
#' university; 9 if finished university; 10 if received graduate (masters or
#' Ph.D) education.  \item \code{cov.female}: Gender. 0 if male; 1 if female.
#' \item \code{civic}: Whether or not the respondent is a member of a civic
#' group in their communities, such as youth groups , women's groups, or
#' community development committees.  }
#' @references Blair, G. (2014). "Why do civilians hold bargaining power in
#' state revenue conflicts? Evidence from Nigeria."
#' 
#' Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design and Analysis
#' of the Randomized Response Technique."  \emph{Working Paper.} Available at
#' \url{http://imai.princeton.edu/research/randresp.html}.
#' @source Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) Replication
#' data for: Design and Analysis of the Randomized Response Technique.
#' @keywords dataset
NULL





#' Nigeria Randomized Response Survey Experiment on Social Connections to Armed
#' Groups
#' 
#' This data set is a subset of the data from the randomized response technique
#' survey experiment conducted in Nigeria to study civilian contact with armed
#' groups. The survey was implemented by Blair (2014).
#' 
#' 
#' @format A data frame containing 2457 observations.  The variables are:
#' \itemize{ \item \code{Quesid}: Survey ID of civilian respondent.  \item
#' \code{rr.q1}: Randomized response survey item using the Forced Response
#' Design asking the respondent whether they hold direct social connections
#' with members of armed groups. 0 if no connection; 1 if connection.  \item
#' \code{cov.age}: Age of the respondent.  \item \code{cov.asset.index}: The
#' number of assets owned by the respondent from an index of nine assets
#' including radio, T.V., motorbike, car, mobile phone, refrigerator, goat,
#' chicken, and cow.  \item \code{cov.married}: Marital status. 0 if single; 1
#' if married.  \item \code{cov.education}: Education level of the respondent.
#' 1 if no school; 2 if started primary school; 3 if finished primary school; 4
#' if started secondary school; 5 if finished secondary school; 6 if started
#' polytechnic or college; 7 if finished polytechnic or college; 8 if started
#' university; 9 if finished university; 10 if received graduate (masters or
#' Ph.D) education.  \item \code{cov.female}: Gender. 0 if male; 1 if female.
#' \item \code{civic}: Whether or not the respondent is a member of a civic
#' group in their communities, such as youth groups , women's groups, or
#' community development committees.  }
#' @references Blair, G. (2014). "Why do civilians hold bargaining power in
#' state revenue conflicts? Evidence from Nigeria."
#' 
#' Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design and Analysis
#' of the Randomized Response Technique."  \emph{Working Paper.} Available at
#' \url{http://imai.princeton.edu/research/randresp.html}.
#' @source Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) Replication
#' data for: Design and Analysis of the Randomized Response Technique.
#' @keywords dataset
NULL





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
#' @aliases rr rr-package
#' @docType package
#' @author Graeme Blair, Experiments in Governance and Politics, Columbia
#' University \email{graeme.blair@@columbia.edu}, \url{http://graemeblair.com}
#' 
#' Kosuke Imai, Department of Politics, Princeton University
#' \email{kimai@@princeton.edu}, \url{http://imai.princeton.edu}
#' 
#' Yang-Yang Zhou, Department of Politics, Princeton University
#' \email{yz3@@princeton.edu}, \url{http://yangyangzhou.com}
#' 
#' Maintainer: Graeme Blair <graeme.blair@@columbia.edu>
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2015) "Design
#' and Analysis of the Randomized Response Technique."  \emph{Journal of the
#' American Statistical Association.} Available at
#' \url{http://graemeblair.com/papers/randresp.pdf}.
#' @keywords package
NULL



