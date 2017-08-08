#' Power Analysis for Randomized Response
#' 
#' \code{power.rr.test} is used to conduct power analysis for randomized
#' response survey designs.
#' 
#' This function allows users to conduct power analysis for randomized response
#' survey designs, both for the standard designs ("forced-known", "mirrored",
#' "disguised", "unrelated-known") and modified designs ("forced-unknown", and
#' "unrelated -unknown").
#' 
#' @param p The probability of receiving the sensitive question (Mirrored
#' Question Design, Unrelated Question Design); the probability of answering
#' truthfully (Forced Response Design); the probability of selecting a red card
#' from the 'yes' stack (Disguised Response Design).
#' @param p0 The probability of forced 'no' (Forced Response Design).
#' @param p1 The probability of forced 'yes' (Forced Response Design).
#' @param q The probability of answering 'yes' to the unrelated question, which
#' is assumed to be independent of covariates (Unrelated Question Design).
#' @param design Call of design (including modified designs) used:
#' "forced-known", "mirrored", "disguised", "unrelated-known",
#' "forced-unknown", and "unrelated-unknown".
#' @param n Number of observations. Exactly one of 'n' or 'power' must be NULL.
#' @param r For the modified designs only (i.e. "forced-unknown" for Forced
#' Response with Unknown Probability and "unrelated-unknown" for Unrelated
#' Question with Unknown Probability), \code{r} is the proportion of
#' respondents allocated to the first group, which is the group that is
#' directed to answer the sensitive question truthfully with probability
#' \code{p} as opposed to the second group which is directed to answer the
#' sensitive question truthfully with probability \code{1-p}.
#' @param presp For a one sample test, the probability of possessing the
#' sensitive trait under the alternative hypothesis.
#' @param presp.null For a one sample test, the probability of possessing the
#' sensitive trait under the null hypothesis. The default is \code{NULL}
#' meaning zero probability of possessing the sensitive trait.
#' @param sig.level Significance level (Type I error probability).
#' @param prespT For a two sample test, the probability of the treated group
#' possessing the sensitive trait under the alternative hypothesis.
#' @param prespC For a two sample test, the probability of the control group
#' possessing the sensitive trait under the alternative hypothesis.
#' @param prespT.null For a two sample test, the probability of the treated
#' group possessing the sensitive trait under the null hypothesis. The default
#' is \code{NULL} meaning there is no difference between the treated and
#' control groups, specifically that \code{prespT.null} is the same as
#' \code{prespC.null}, the probability of the control group possessing the
#' sensitive trait under the null hypothesis.
#' @param prespC.null For a two sample test, the probability of the control
#' group possessing the sensitive trait under the null hypothesis.
#' @param power Power of test (Type II error probability). Exactly one of 'n'
#' or 'power' must be NULL.
#' @param type One or two sample test. For a two sample test, the alternative
#' and null hypotheses refer to the difference between the two samples of the
#' probabilities of possessing the sensitive trait.
#' @param alternative One or two sided test.
#' @param solve.tolerance When standard errors are calculated, this option
#' specifies the tolerance of the matrix inversion operation solve.
#' @return \code{power.rr.test} contains the following components (the
#' inclusion of some components such as the design parameters are dependent
#' upon the design used):
#' 
#' \item{n}{Point estimates for the effects of covariates on the randomized
#' response item.} \item{r}{Standard errors for estimates of the effects of
#' covariates on the randomized response item.} \item{presp}{For a one sample
#' test, the probability of possessing the sensitive trait under the
#' alternative hypothesis.  For a two sample test, the difference between the
#' probabilities of possessing the sensitive trait for the treated and control
#' groups under the alternative hypothesis.} \item{presp.null}{For a one sample
#' test, the probability of possessing the sensitive trait under the null
#' hypothesis. For a two sample test, the difference between the probabilities
#' of possessing the sensitive trait for the treated and control groups under
#' the null hypothesis.} \item{sig.level}{Significance level (Type I error
#' probability).} \item{power}{Power of test (Type II error probability).}
#' \item{type}{One or two sample test.} \item{alternative}{One or two sided
#' test.}
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2015) "Design
#' and Analysis of the Randomized Response Technique."  \emph{Journal of the
#' American Statistical Association.} Available at
#' \url{http://graemeblair.com/papers/randresp.pdf}.
#' @keywords analysis power
#' @examples
#' 
#' 
#' 
#' ## Calculate the power to detect a sensitive item proportion of .2
#' ## with the forced design with known probabilities of 2/3 in truth-telling group,
#' ## 1/6 forced to say "yes" and 1/6 forced to say "no" and sample size of 200.
#' 
#' power.rr.test(p = 2/3, p1 = 1/6, p0 = 1/6, n = 200, 
#'              presp = .2, presp.null = 0,
#'              design = "forced-known", sig.level = .01,
#'              type = "one.sample", alternative = "one.sided")
#' 				       
#' \dontrun{
#' 
#' ## Find power varying the number of respondents from 250 to 2500 and 
#' ## the population proportion of respondents possessing the sensitive 
#' ## trait from 0 to .15
#' 
#' presp.seq <- seq(from = 0, to = .15, by = .0025)
#' n.seq <- c(250, 500, 1000, 2000, 2500)
#' power <- list()
#' for(n in n.seq) {
#'     power[[n]] <- rep(NA, length(presp.seq))
#'     for(i in 1:length(presp.seq))
#'         power[[n]][i] <- power.rr.test(p = 2/3, p1 = 1/6, p0 = 1/6, n = n, 
#'                                        presp = presp.seq[i], presp.null = 0,
#'                                        design = "forced-known", sig.level = .01, 
#'                                        type = "one.sample",
#'                                        alternative = "one.sided")$power
#'     }
#'     
#' ## Replicates the results for Figure 2 in Blair, Imai, and Zhou (2014)
#' }
#' 
#' 
#' @export power.rr.test
power.rr.test <- function(p, p0, p1, q, design, n = NULL, r, presp, presp.null = NULL, sig.level, 
                          prespT, prespC, prespT.null = NULL, prespC.null,
                          power = NULL, type = c("one.sample", "two.sample"), 
                          alternative = c("one.sided", "two.sided"),
                          solve.tolerance = .Machine$double.eps
                          ){

se.f.rr <- function(p, p0, p1, q, design, n, r, presp){
  
  if(missing(presp))
    stop("Missing probability of latent response, see documentation")
  if(missing(n))
    stop("Missing sample size, see documentation")
  
  # Get c, d
  if (design != "forced-unknown" & design != "unrelated-unknown") {
    cd <- rrcd(p, p0, p1, q, design)
    c <- cd$c
    d <- cd$d
  } else {
    if(missing(r)){
      r <- 1/2
    }
    if (any((r > 1) | (r < 0))) {
      stop("r must be a probability")
    }
  }
    
  if(design == "forced-known" | design == "unrelated-known"){
    
    information.rr <- c^2  * ((c * presp + d)^(-1) + (1 - c * presp - d)^(-1))
    vcov <- solve(information.rr, tol = solve.tolerance) / n
    se <- sqrt(diag(vcov))
    
  } else if(design == "mirrored" | design == "disguised"){
    
    information.rr <- c^2 * (((c * presp + d)^(-2) * ((1 - 2*d) * presp  + d)) +
          ((1 - c * presp - d)^(-2) * ((2*d - 1) * presp + (1 - d))))
    vcov <- solve(information.rr, tol = solve.tolerance) / n
    se <- sqrt(diag(vcov))
    
  } else if (design == "unrelated-unknown") {
    
    if(missing(p) & !missing(q))
      stop(paste("for the", design, "design, you must specify p"))
    if(missing(q) & !missing(p))
      stop(paste("for the", design, "design, you must specify q"))
    if(missing(p) & missing(q))
      stop(paste("for the", design, "design, you must specify p and q"))
    
    if (any((p > 1) | (p < 0))) {
      stop("p must be a probability")
    }
    
    if (any((q > 1) | (q < 0))) {
      stop("q must be a probability")
    }
    
    information.rr <- ((r*(p^2)) / ((p*presp + (1 - p)*q)*(1 - p*presp - (1 - p)*q))) +
                        (((1-r)*((1 - p)^2)) / (((1 - p)*presp + p*q)*(1 - (1 - p)*presp - p*q)))
    
    vcov <- solve(information.rr, tol = solve.tolerance) / n
    se <- sqrt(diag(vcov))
    
  } else if (design == "forced-unknown") {
    
    if(missing(p) & !missing(p1))
      stop("for the forced design, you must specify p")
    if(missing(p1) & !missing(p))
      stop("for the forced design, you must specify p1")
    if(missing(p) & missing(p1))
      stop("for the forced design, you must specify p and p1")
    
    if (any((p > 1) | (p < 0))) {
      stop("p must be a probability")
    }
    if (any((p1 > 1) | (p1 < 0))) {
      stop("p1 must be a probability")
    }
    if(missing(p0)){
      p0 <- 0
    }
    
    if (any((p0 > 1) | (p0 < 0))) {
      stop("p0 must be a probability")
    }
    if (any( round(p + p0 + p1, 10) !=1)) {
      stop("p, p0, and p1 must sum to 1")
    }
    
    if(missing(q)){
      q <- 1
    }
    if (any((q > 1) | (q < 0))) {
      stop("q must be a probability")
    }
    
    information.rr <- ((r*(p^2)) / ((p*presp + (1 - p)*q)*(1 - p*presp - (1 - p)*q))) +
      (((1-r)*((1 - p)^2)) / (((1 - p)*presp + p*q)*(1 - (1 - p)*presp - p*q)))
    
    vcov <- solve(information.rr, tol = solve.tolerance) / n
    se <- sqrt(diag(vcov))
    
  }
  
  return.object <- list(vcov = vcov,
                        se = se)
  return(return.object)
  
}

    
  if(missing(sig.level))
    stop("Missing significance level (Type I error probability)")
  if(missing(type))
    stop("Missing specification, one.sample or two.sample")
  if(missing(alternative))
    stop("Missing specification, one.sided or two.sided test")
  if((missing(n) & missing(power)) | (is.null(n) & is.null(power)))
    stop("Exactly one of 'n' or 'power' must be NULL")
  if((!missing(n) & !missing(power)) | (is.null(n) & is.null(power)))
    stop("Exactly one of 'n' or 'power' must be NULL") 
  if (any((sig.level > 1) | (sig.level < 0))) {
    stop("Significance level must be a probability")
  }
  if (any((power > 1) | (power < 0))) {
    stop("Power must be a probability")
  }
  
  if (design != "forced-unknown" & design != "unrelated-unknown") {
    r <- NULL
  } else {
    if(missing(r)){
        r <- 1/2
      }
    if (any((r > 1) | (r < 0))) {
      stop("r must be a probability")
      }
  }
  
  ## Calculating power varying n and presp
  if(is.null(power)){
  # one sample
    if(type == "one.sample"){
      
      # for one sample, presp.null is 0 unless otherwise specified 
      if(is.null(presp.null)){
        presp.null <- 0
      }
      
      if (any((presp > 1) | (presp < 0))) {
        stop("Proportion of respondents with sensitive trait must be a probability")
      }
      if (any((presp.null > 1) | (presp.null < 0))) {
        stop("Null proportion of respondents with sensitive trait must be a probability")
      }
      
      se.rr.alt <- se.f.rr(p, p0, p1, q, design, n, r, presp)
      se <- se.rr.alt$se
    
      se.rr.null <- se.f.rr(p, p0, p1, q, design, n, r, presp.null)
      se.null <- se.rr.null$se
    
    # one.sided
      if(alternative == "one.sided") {
      
        thld.null <- presp.null + qnorm(1 - (sig.level)) * se.null
        power <- pnorm(thld.null, presp, se, lower.tail = FALSE)
      
    # two.sided
      } else if(alternative == "two.sided"){
      
        thld.low <- presp.null - qnorm(1 - (sig.level/2)) * se.null
        thld.high <- presp.null + qnorm(1 - (sig.level/2)) * se.null
        power.low <- pnorm(thld.low, presp, se, lower.tail = TRUE)
        power.high <- pnorm(thld.high, presp, se, lower.tail = FALSE)
        power <- power.low + power.high
      
      }
    } else if(type == "two.sample"){
      
      # for two sample, prespT.null is prespC.null unless otherwise specified
      if(is.null(prespT.null)){
        prespT.null <- prespC.null
      }
      
      if (any((prespT > 1) | (prespT < 0) | (prespC > 1) | (prespC < 0))) {
        stop("Proportion of respondents with sensitive trait must be a probability")
      }
      if (any((prespT.null > 1) | (prespT.null < 0) | (prespC.null > 1) | (prespC.null < 0))) {
        stop("Null proportion of respondents with sensitive trait must be a probability")
      }
      
      presp <- prespT - prespC
      se.rr.alt.T <- se.f.rr(p, p0, p1, q, design, n, r, prespT)$se
      se.rr.alt.C <- se.f.rr(p, p0, p1, q, design, n, r, prespC)$se
      se <- sqrt(se.rr.alt.T^2 + se.rr.alt.C^2)
            
      presp.null <- prespT.null - prespC.null
      se.rr.null.T <- se.f.rr(p, p0, p1, q, design, n, r, prespT.null)$se
      se.rr.null.C <- se.f.rr(p, p0, p1, q, design, n, r, prespC.null)$se
      se.null <- sqrt(se.rr.null.T^2 + se.rr.null.C^2)
      
      # one.sided
      if(alternative == "one.sided") {
        
        thld.null <- presp.null + qnorm(1 - (sig.level)) * se.null
        power <- pnorm(thld.null, presp, se, lower.tail = FALSE)
        
        # two.sided
      } else if(alternative == "two.sided"){
        
        thld.low <- presp.null - qnorm(1 - (sig.level/2)) * se.null
        thld.high <- presp.null + qnorm(1 - (sig.level/2)) * se.null
        power.low <- pnorm(thld.low, presp, se, lower.tail = TRUE)
        power.high <- pnorm(thld.high, presp, se, lower.tail = FALSE)
        power <- power.low + power.high
        
      }      
    }
  }
  
  ## Calculating n given power and presp
  if(is.null(n)){
    
    # Parameters for design types
  if(design != "forced-unknown" & design != "unrelated-unknown") {
    N <- 1:100000
  } else {
    N <- 2:100000
  }
  
    # one sample
    if(type == "one.sample"){
      
      # for one sample, presp.null is 0 unless otherwise specified 
      if(is.null(presp.null)){
        presp.null <- 0
      }
      
      if (any((presp > 1) | (presp < 0))) {
        stop("Proportion of respondents with sensitive trait must be a probability")
      }
      if (any((presp.null > 1) | (presp.null < 0))) {
        stop("Null proportion of respondents with sensitive trait must be a probability")
      }
      
      # one.sided
      if(alternative == "one.sided") {
        
        for(i in 1:length(N)){
        
          sample <- N[i]
          se.rr.alt <- se.f.rr(p, p0, p1, q, design, sample, r, presp)
          se <- se.rr.alt$se
          se.rr.null <- se.f.rr(p, p0, p1, q, design, sample, r, presp.null)
          se.null <- se.rr.null$se
          thld.null <- presp.null + qnorm(1 - (sig.level)) * se.null
          power.n <- pnorm(thld.null, presp, se, lower.tail = FALSE)
          
          if(abs(power - power.n) < .01){
            n <- N[i] 
            break
          }
          
          if(sample == length(N)){
            stop("Your n exceeds 100,000")
          }
        }
        # two.sided
      } else if(alternative == "two.sided"){
        
        for(i in 1:length(N)){
          
          sample <- N[i]
          se.rr.alt <- se.f.rr(p, p0, p1, q, design, sample, r, presp)
          se <- se.rr.alt$se
          se.rr.null <- se.f.rr(p, p0, p1, q, design, sample, r, presp.null)
          se.null <- se.rr.null$se
          thld.low <- presp.null - qnorm(1 - (sig.level/2)) * se.null
          thld.high <- presp.null + qnorm(1 - (sig.level/2)) * se.null
          power.low <- pnorm(thld.low, presp, se)
          power.high <- pnorm(thld.high, presp, se, lower.tail = FALSE)
          power.n <- power.low + power.high
        
          if(abs(power - power.n) < .01){
            n <- N[i] 
            break
          }
          
          if(sample == length(N)){
            stop("Your n exceeds 100,000")
          }
        }
      }
    } else if(type == "two.sample"){
      
      # for two sample, prespT.null is prespC.null unless otherwise specified
      if(is.null(prespT.null)){
        prespT.null <- prespC.null
      }
      
      if (any((prespT > 1) | (prespT < 0) | (prespC > 1) | (prespC < 0))) {
        stop("Proportion of respondents with sensitive trait must be a probability")
      }
      if (any((prespT.null > 1) | (prespT.null < 0) | (prespC.null > 1) | (prespC.null < 0))) {
        stop("Null proportion of respondents with sensitive trait must be a probability")
      }
      
      # one.sided
      if(alternative == "one.sided") {
        
        for(i in 1:length(N)){
          
          sample <- N[i]
          presp <- prespT - prespC
          se.rr.alt.T <- se.f.rr(p, p0, p1, q, design, sample, r, prespT)$se
          se.rr.alt.C <- se.f.rr(p, p0, p1, q, design, sample, r, prespC)$se
          se <- sqrt(se.rr.alt.T^2 + se.rr.alt.C^2)
          
          presp.null <- prespT.null - prespC.null
          se.rr.null.T <- se.f.rr(p, p0, p1, q, design, sample, r, prespT.null)$se
          se.rr.null.C <- se.f.rr(p, p0, p1, q, design, sample, r, prespC.null)$se
          se.null <- sqrt(se.rr.null.T^2 + se.rr.null.C^2)
          
          thld.null <- presp.null + qnorm(1 - (sig.level)) * se.null
          power.n <- pnorm(thld.null, presp, se, lower.tail = FALSE)
          
          if(abs(power - power.n) < .01){
            n <- N[i] 
            break
          }
          
          if(sample == length(N)){
            stop("Your n exceeds 100,000")
          }
        }
        # two.sided
      } else if(alternative == "two.sided"){
        
        for(i in 1:length(N)){
          
          sample <- N[i]
          presp <- prespT - prespC
          se.rr.alt.T <- se.f.rr(p, p0, p1, q, design, sample, r, prespT)$se
          se.rr.alt.C <- se.f.rr(p, p0, p1, q, design, sample, r, prespC)$se
          se <- sqrt(se.rr.alt.T^2 + se.rr.alt.C^2)
          
          presp.null <- prespT.null - prespC.null
          se.rr.null.T <- se.f.rr(p, p0, p1, q, design, sample, r, prespT.null)$se
          se.rr.null.C <- se.f.rr(p, p0, p1, q, design, sample, r, prespC.null)$se
          se.null <- sqrt(se.rr.null.T^2 + se.rr.null.C^2)
          
          thld.low <- presp.null - qnorm(1 - (sig.level/2)) * se.null
          thld.high <- presp.null + qnorm(1 - (sig.level/2)) * se.null
          power.low <- pnorm(thld.low, presp, se)
          power.high <- pnorm(thld.high, presp, se, lower.tail = FALSE)
          power.n <- power.low + power.high
          
          if(abs(power - power.n) < .01){
            n <- N[i] 
            break
          }
          
          if(sample == length(N)){
            stop("Your n exceeds 100,000")
           }
          }
        }
      }     
    }     
  
  return.object <- list(design = design,
                        n = n,
                        r = r,
                        presp = presp,
                        presp.null = presp.null,
                        sig.level = sig.level,
                        power = power,
                        type = type,
                        alternative = alternative)
  
  return(return.object)
}





#' Power Analysis Plot for Randomized Response
#' 
#' \code{power.rr.plot} generates a power analysis plot for randomized response
#' survey designs.
#' 
#' This function generates a power analysis plot for randomized response survey
#' designs, both for the standard designs ("forced-known", "mirrored",
#' "disguised", "unrelated-known") and modified designs ("forced-unknown", and
#' "unrelated -unknown"). The x-axis shows the population proportions with the
#' sensitive trait; the y-axis shows the statistical power; and different
#' sample sizes are shown as different lines in grayscale.
#' 
#' @param p The probability of receiving the sensitive question (Mirrored
#' Question Design, Unrelated Question Design); the probability of answering
#' truthfully (Forced Response Design); the probability of selecting a red card
#' from the 'yes' stack (Disguised Response Design).
#' @param p0 The probability of forced 'no' (Forced Response Design).
#' @param p1 The probability of forced 'yes' (Forced Response Design).
#' @param q The probability of answering 'yes' to the unrelated question, which
#' is assumed to be independent of covariates (Unrelated Question Design).
#' @param design Call of design (including modified designs) used:
#' "forced-known", "mirrored", "disguised", "unrelated-known",
#' "forced-unknown", and "unrelated-unknown".
#' @param n.seq A sequence of number of observations or sample sizes.
#' @param r For the modified designs only (i.e. "forced-unknown" for Forced
#' Response with Unknown Probability and "unrelated-unknown" for Unrelated
#' Question with Unknown Probability), \code{r} is the proportion of
#' respondents allocated to the first group, which is the group that is
#' directed to answer the sensitive question truthfully with probability
#' \code{p} as opposed to the second group which is directed to answer the
#' sensitive question truthfully with probability \code{1-p}.
#' @param presp.seq For a one sample test, a sequence of probabilities of
#' possessing the sensitive trait under the alternative hypothesis.
#' @param presp.null For a one sample test, the probability of possessing the
#' sensitive trait under the null hypothesis. The default is \code{NULL}
#' meaning zero probability of possessing the sensitive trait.
#' @param sig.level Significance level (Type I error probability).
#' @param prespT.seq For a two sample test, a sequence of probabilities of the
#' treated group possessing the sensitive trait under the alternative
#' hypothesis.
#' @param prespC.seq For a two sample test, a sequence of probabitilies of the
#' control group possessing the sensitive trait under the alternative
#' hypothesis.
#' @param prespT.null For a two sample test, the probability of the treated
#' group possessing the sensitive trait under the null hypothesis. The default
#' is \code{NULL} meaning there is no difference between the treated and
#' control groups, specifically that \code{prespT.null} is the same as
#' \code{prespC.null}, the probability of the control group possessing the
#' sensitive trait under the null hypothesis.
#' @param prespC.null For a two sample test, the probability of the control
#' group possessing the sensitive trait under the null hypothesis.
#' @param type One or two sample test. For a two sample test, the alternative
#' and null hypotheses refer to the difference between the two samples of the
#' probabilities of possessing the sensitive trait.
#' @param alternative One or two sided test.
#' @param solve.tolerance When standard errors are calculated, this option
#' specifies the tolerance of the matrix inversion operation solve.
#' @param legend Indicator of whether to include a legend of sample sizes. The
#' default is \code{TRUE}.
#' @param legend.x Placement on the x-axis of the legend. The default is
#' \code{"bottomright"}.
#' @param legend.y Placement on the y-axis of the legend.
#' @param par Option to set or query graphical parameters within the function.
#' The default is \code{TRUE}.
#' @param ... Additional arguments to be passed to \code{par()}
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design
#' and Analysis of the Randomized Response Technique."  \emph{Working Paper.}
#' Available at \url{http://imai.princeton.edu/research/randresp.html}.
#' @keywords analysis power
#' @examples
#' 
#' 
#' ## Generate a power plot for the forced design with known 
#' ## probabilities of 2/3 in truth-telling group, 1/6 forced to say "yes" 
#' ## and 1/6 forced to say "no", varying the number of respondents from 
#' ## 250 to 2500 and the population proportion of respondents 
#' ## possessing the sensitive trait from 0 to .15.
#' 
#' 
#' presp.seq <- seq(from = 0, to = .15, by = .0025)
#' n.seq <- c(250, 500, 1000, 2000, 2500)
#' power.rr.plot(p = 2/3, p1 = 1/6, p0 = 1/6, n.seq = n.seq, 
#'               presp.seq = presp.seq, presp.null = 0,
#'               design = "forced-known", sig.level = .01, 
#'               type = "one.sample",
#'               alternative = "one.sided", legend = TRUE)
#' 				       
#'     
#' ## Replicates the results for Figure 2 in Blair, Imai, and Zhou (2014)
#' 
#' 
#' @export power.rr.plot
power.rr.plot <- function(p, p0, p1, q, design, n.seq, r, presp.seq, presp.null = NULL, sig.level, 
                          prespT.seq, prespC.seq, prespT.null = NULL, prespC.null,
                          type = c("one.sample", "two.sample"), 
                          alternative = c("one.sided", "two.sided"),
                          solve.tolerance = .Machine$double.eps, legend = TRUE, legend.x = "bottomright", 
                          legend.y, par = TRUE, ...
){
  if(missing(n.seq))
    stop("Please provide a vector of sample sizes")
  
  if(missing(presp.seq))
    stop("Please provide a vector of proportions with the sensitive trait")
  
  power <- list()
  
  if(missing(presp.null))
    presp.null <- NULL
  
  if(missing(prespT.null))
    prespT.null <-NULL
  
  if(alternative == "one.sided"){
    alt <- "One Sided"
  }else{
    alt <- "Two Sided"
  }
  
  if(type == "one.sample"){
  for(n in n.seq) {
    power[[n]] <- rep(NA, length(presp.seq))
    for(i in 1:length(presp.seq))
      power[[n]][i] <- power.rr.test(p = p, p0 = p0, p1 = p1, q = q, n = n, 
                                     presp = presp.seq[i], presp.null = presp.null,
                                     design = design, sig.level = sig.level, 
                                     type = type,
                                     alternative = alternative)$power
  }
  }else{
    if(length(prespT.seq) != length(prespC.seq))
      stop("length of vector prespT.seq is not equal to length of vector prespC.seq")
    
    presp.seq <- prespT.seq
    
    for(n in n.seq) {
      power[[n]] <- rep(NA, length(presp.seq))
      for(i in 1:length(presp.seq))
        power[[n]][i] <- power.rr.test(p = p, p0 = p0, p1 = p1, q = q, n = n, 
                                       prespT = presp.seq[i], 
                                       prespC = prespC.seq[i], 
                                       prespT.null = prespT.null, 
                                       prespC.null = prespC.null,
                                       design = design, sig.level = sig.level, 
                                       type = type,
                                       alternative = alternative)$power
  }
  }
  
  if(par == TRUE)
  par(oma = c(0, 0, 0, 0), mar = c(3.6, 3.6, 0, 0), las = 1, mgp = c(2, .7, 0), tck = -.01, cex = 0.8)  
  
  plot(0, 1, type = "n", xlim = c(0, max(presp.seq)), ylim = c(0, 1.05), axes = F, 
       xlab = "", ylab = "")
  
  abline(h = 0, lty = "dashed", col = "darkgray")
  
  col.seq <- rev(gray((n.seq)/(3000+150)))
  
  for(n in 1:length(n.seq))
    lines(presp.seq, power[[n.seq[n]]], lwd = 1.5, col = col.seq[n])
  
  axis(1)
  axis(2)
  mtext("Population Proportion of Respondents", side = 1, las = 0, line = 2.5, cex = 0.8)
  mtext(paste("Statistical Power for", alt, "Test"), side = 2, las = 0, line = 2.5, cex = 0.8)
  
  if(legend == TRUE){
    
    legend(x = legend.x, y = legend.y, legend = n.seq, title = "Sample size", lwd = 1, col = col.seq, bty = "n", cex = 0.8)

  }
}
