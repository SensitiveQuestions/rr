\name{rr-package}
\alias{rr-package}
\alias{rr}
\docType{package}
\title{R Package for the Randomized Response Technique}
\description{
  \code{rr} implements methods developed by Blair, Imai, and Zhou (2015) such as     
  multivariate regression and power analysis for the randomized response technique. Randomized response is a survey technique 
  that introduces random noise to reduce potential bias from non-response and social desirability when asking questions about 
  sensitive behaviors and beliefs. The current version of this package conducts multivariate regression analyses for the 
  sensitive item under four standard randomized response designs: mirrored question, forced response, disguised response, and 
  unrelated question. Second, it generates predicted probabilities of answering affirmatively to the sensitive item for each 
  respondent. Third, it also allows users to use the sensitive item as a predictor in an outcome regression under the forced 
  response design. Additionally, it implements power analyses to help improve research design. In future versions, this 
  package will extend to new modified designs that are based on less stringent assumptions than those of the standard 
  designs, specifically to allow for non-compliance and unknown distribution to the unrelated question under the 
  unrelated question design.
}
\details{
\tabular{ll}{
Package: \tab rr\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2015-3-8\cr
License: \tab GPL (>= 2)\cr
}
}

\author{
  Graeme Blair, Experiments in Governance and Politics, Columbia University 
  \email{graeme.blair@columbia.edu}, \url{http://graemeblair.com}
  
  Kosuke Imai, Department of Politics, Princeton University
  \email{kimai@princeton.edu}, \url{http://imai.princeton.edu}
  
  Yang-Yang Zhou, Department of Politics, Princeton University
  \email{yz3@princeton.edu}, \url{http://yangyangzhou.com}
  
  Maintainer: Graeme Blair <graeme.blair@columbia.edu>
}

\references{
  Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design and Analysis of the Randomized Response Technique."   
  \emph{Working Paper.} Available at \url{http://imai.princeton.edu/research/randresp.html}.
}
\keyword{package}
