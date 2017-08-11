logistic <- function(x) exp(x)/(1+exp(x))



#' Randomized Response as a Regression Predictor
#' 
#' \code{rrreg.predictor} is used to jointly model the randomized response item
#' as both outcome and predictor for an additional outcome given a set of
#' covariates.
#' 
#' This function allows users to perform multivariate regression analysis with
#' the randomized response item as a predictor for a separate outcome of
#' interest. It does so by jointly modeling the randomized response item as
#' both outcome and predictor for an additional outcome given the same set of
#' covariates. Four standard designs are accepted by this function: mirrored
#' question, forced response, disguised response, and unrelated question.
#' 
#' @usage rrreg.predictor(formula, p, p0, p1, q, design, data, rr.item,
#' model.outcome = "logistic", fit.sens = "bayesglm", fit.outcome = "bayesglm",
#' bstart = NULL, tstart = NULL, parstart = TRUE, maxIter = 10000, verbose =
#' FALSE, optim = FALSE, em.converge = 10^(-4), glmMaxIter = 20000, estconv =
#' TRUE, solve.tolerance = .Machine$double.eps)
#' @param formula An object of class "formula": a symbolic description of the
#' model to be fitted with the randomized response item as one of the
#' covariates.
#' @param p The probability of receiving the sensitive question (Mirrored
#' Question Design, Unrelated Question Design); the probability of answering
#' truthfully (Forced Response Design); the probability of selecting a red card
#' from the 'yes' stack (Disguised Response Design).
#' @param p0 The probability of forced 'no' (Forced Response Design).
#' @param p1 The probability of forced 'yes' (Forced Response Design).
#' @param q The probability of answering 'yes' to the unrelated question, which
#' is assumed to be independent of covariates (Unrelated Question Design).
#' @param design One of the four standard designs: "forced-known", "mirrored",
#' "disguised", or "unrelated-known".
#' @param data A data frame containing the variables in the model. Observations
#' with missingness are list-wise deleted.
#' @param rr.item A string containing the name of the randomized response item
#' variable in the data frame.
#' @param model.outcome Currently the function only allows for logistic
#' regression, meaning the outcome variable must be binary.
#' @param fit.sens Indicator for whether to use Bayesian generalized linear
#' modeling (bayesglm) in the Maximization step for the
#' Expectation-Maximization (EM) algorithm to generate coefficients for the
#' randomized response item as the outcome.  Default is \code{"bayesglm"};
#' otherwise input \code{"glm"}.
#' @param fit.outcome Indicator for whether to use Bayesian generalized linear
#' modeling (bayesglm) in the Maximization step for the EM algorithm to
#' generate coefficients for the outcome variable given in the formula with the
#' randomized response item as a covariate. Default is \code{"bayesglm"};
#' otherwise input \code{"glm"}.
#' @param bstart Optional starting values of coefficient estimates for the
#' randomized response item as outcome for the EM algorithm.
#' @param tstart Optional starting values of coefficient estimates for the
#' outcome variable given in the formula for the EM algorithm.
#' @param parstart Option to use the function \code{rrreg} to generate starting
#' values of coefficient estimates for the randomized response item as outcome
#' for the EM algorithm. The default is \code{TRUE}, but if starting estimates
#' are inputted by the user in \code{bstart}, this option is overidden.
#' @param maxIter Maximum number of iterations for the Expectation-Maximization
#' algorithm. The default is \code{10000}.
#' @param verbose A logical value indicating whether model diagnostics counting
#' the number of EM iterations are printed out.  The default is \code{FALSE}.
#' @param optim A logical value indicating whether to use the quasi-Newton
#' "BFGS" method to calculate the variance-covariance matrix and standard
#' errors. The default is \code{FALSE}.
#' @param em.converge A value specifying the satisfactory degree of convergence
#' under the EM algorithm. The default is \code{10^(-4)}.
#' @param glmMaxIter A value specifying the maximum number of iterations to run
#' the EM algorithm. The default is \code{20000} .
#' @param estconv Option to base convergence on the absolute value of the
#' difference between subsequent coefficients generated through the EM
#' algorithm rather than the subsequent log-likelihoods. The default is
#' \code{TRUE}.
#' @param solve.tolerance When standard errors are calculated, this option
#' specifies the tolerance of the matrix inversion operation solve.
#' @return \code{rrreg.predictor} returns an object of class "rrpredreg"
#' associated with the randomized response item as predictor.  The object
#' \code{rrpredreg} is a list that contains the following components (the
#' inclusion of some components such as the design parameters are dependent
#' upon the design used):
#' 
#' \item{est.t}{Point estimates for the effects of the randomized response item
#' as predictor and other covariates on the separate outcome variable specified
#' in the formula.} \item{se.t}{Standard errors for estimates of the effects of
#' the randomized response item as predictor and other covariates on the
#' separate outcome variable specified in formula.} \item{est.b}{Point
#' estimates for the effects of covariates on the randomized response item.}
#' \item{vcov}{Variance-covariance matrix for estimates of the effects of the
#' randomized response item as predictor and other covariates on the separate
#' outcome variable specified in formula as well as for estimates of the
#' effects of covariates on the randomized response item.} \item{se.b}{Standard
#' errors for estimates of the effects of covariates on the randomized response
#' item.} \item{data}{The \code{data} argument.} \item{coef.names}{Variable
#' names as defined in the data frame.} \item{x}{The model matrix of
#' covariates.} \item{y}{The randomized response vector.} \item{o}{The separate
#' outcome of interest vector.} \item{design}{Call of standard design used:
#' "forced-known", "mirrored", "disguised", or "unrelated-known".} \item{p}{The
#' \code{p} argument.} \item{p0}{The \code{p0} argument.} \item{p1}{The
#' \code{p1} argument.} \item{q}{The \code{q} argument.} \item{call}{The
#' matched call.}
#' @seealso \code{\link{rrreg}} for multivariate regression.
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design
#' and Analysis of the Randomized Response Technique."  \emph{Working Paper.}
#' Available at \url{http://imai.princeton.edu/research/randresp.html}.
#' @keywords regression predictor joint model
#' @examples
#' 
#' data(nigeria)
#' 
#' ## Define design parameters
#' 
#' set.seed(44)
#' 
#' p <- 2/3  # probability of answering honestly in Forced Response Design
#' p1 <- 1/6 # probability of forced 'yes'
#' p0 <- 1/6 # probability of forced 'no'
#' 
#' ## Fit joint model of responses to an outcome regression of joining a civic 
#' ## group and the randomized response item of having a militant social connection
#' \dontrun{
#' rr.q1.pred.obj <- 
#'     rrreg.predictor(civic ~ cov.asset.index + cov.married + I(cov.age/10) + 
#'               I((cov.age/10)^2) + cov.education + cov.female 
#'               + rr.q1, rr.item = "rr.q1", parstart = FALSE, estconv = TRUE,
#'               data = nigeria, verbose = FALSE, optim = TRUE,
#'               p = p, p1 = p1, p0 = p0, design = "forced-known")
#' 
#' summary(rr.q1.pred.obj)
#' }
#' ## Replicates Table 4 in Blair, Imai, and Zhou (2014)
#' 
#' @importFrom arm bayesglm
#' @export
rrreg.predictor <- function(formula, p, p0, p1, q, design, data, rr.item, model.outcome = "logistic",
                            fit.sens = "bayesglm", fit.outcome = "bayesglm",
                            bstart = NULL, tstart = NULL, parstart = TRUE, maxIter = 10000, verbose = FALSE, 
                            optim = FALSE, em.converge = 10^(-4), glmMaxIter = 20000,
                            estconv = TRUE, solve.tolerance = .Machine$double.eps) {
  
  # Set up data objects
  df.na <- model.frame(formula, data)
  df <- df.na[complete.cases(df.na),]
  
  xy <- model.matrix.default(formula, df)
  x1 <- xy[, colnames(xy) != paste(rr.item)] 
  
  bscale <- c(1, unname(apply((as.matrix(x1[,-1])), 2, sd)))
  tscale <- c(bscale, 1)
  scale <- c(bscale, tscale)
  
  x <- sweep(x1, 2, bscale, `/`) #rescale/standardize each covar, not including intercept
  o <- model.response(df) # outcome variable
  y <- df[,paste(rr.item)] # sensitive item
  n <- length(o)
  
  xg1 <- cbind(x, 1) # data inputs for the outcome model given x and z = 1, 0, y
  xg0 <- cbind(x, 0)
  xgy <- cbind(x, y)
  
  dat <- as.data.frame(xgy)
  
  # Get c, d
  #cd <- rrcd(p, p0, p1, q, design)
  #c <- cd$c
  #d <- cd$d
  
  # For optim, must merge bpar and tpar to be one vector of estimates to optimize, within the function, we separate them
  obs.llik <- function(par, y, x, o, p, p1, xg1, xg0, fit.sens, fit.outcome) {
      
      bpar <- par[1:ncol(x)]
      tpar <- par[ncol(x)+1:ncol(xg1)]
      fx <- logistic(x %*% bpar)
      gx1 <- (o == 1)*logistic(xg1 %*% tpar) + (o == 0)*(1 - logistic(xg1 %*% tpar))
      gx0 <- (o == 1)*logistic(xg0 %*% tpar) + (o == 0)*(1 - logistic(xg0 %*% tpar))

      bayes <- 0
      if(fit.sens == "bayesglm")
          bayes <- bayes + sum(dcauchy(x = tpar, scale = rep(2.5, length(tpar)), log = TRUE))
      if(fit.outcome == "bayesglm")
          bayes <- bayes + sum(dcauchy(x = bpar, scale = rep(2.5, length(bpar)), log = TRUE))
      
      if((1 - p - p1) == 0){ ## when probability of forced "no" is 0
        forceno <- 1
      } else {
        forceno <- log(1 - p - p1)
      }
      
      if(p1 == 0){ ## when probability of forced "yes" is 0
        forceyes <- 1
      } else {
        forceyes <- log(p1)
      }
            
     llik <- sum(log( exp(log(fx) + log(gx1) + y * log(p + p1) + (1-y) * forceno) +
                       exp(log(1 - fx) + log(gx0) + y * forceyes + (1-y) * log(1 - p1)) ) ) + bayes
     
    return(llik)
  }
  
  ### NOW THE EM ALGORITHM ###
  
  ## set some starting values for bpar and tpar ##
  if(!missing(bstart) & (length(bstart) != ncol(x)))
    stop("The starting values are not the right dimensions. bstart should be a vector of length the number of variables,
         excluding the randomized response item and including the intercept if applicable.")
  
  if(!missing(tstart) & (length(tstart) != ncol(xy)))
    stop("The starting values are not the right dimensions. tstart should be a vector of length the number of variables,
         including the intercept if applicable.")
  
  
  if(!is.null(bstart)){ 
    bpar <- bstart
  } else {
    if(parstart == TRUE) {
      bpar <- rrreg(y ~ x - 1, p = p, p1 = p1, p0 = p0, data = dat, design = "forced-known")$est
    } else {
      bpar <- rnorm(ncol(x))
    }
  }
  
  if(!is.null(tstart)){
    tpar <- tstart
  } else {
    if(parstart == TRUE) {
      bfit <- rrreg(y ~ x - 1, p = p, p1 = p1, p0 = p0, data = dat, design = "forced-known")
      propz <- predict.rrreg(bfit, avg = FALSE, given.y = TRUE, quasi.bayes = "single")$est
      predz <- ifelse(as.numeric(unlist(propz)) < .5, 0, 1)
      xgz <- cbind(x, predz)
      tpar <- glm(cbind(o, 1-o) ~ xgz - 1, family = binomial("logit"))$coef  
    } else {
      tpar <- rnorm(ncol(xg1))
    }
  }
  
  ## start off with an infinitely negative log likelihood
  pllik <- -Inf
  
  ## calculate the log likelihood at the starting values
  par <- c(bpar, tpar)
  llik <- obs.llik(par = par, y = y, x = x, o = o, p = p, p1 = p1, xg1 = xg1, xg0 = xg0, fit.sens = fit.sens, fit.outcome = fit.outcome)
  
  ## set a counter to zero, which will iterate
  
  counter <- 0
  
  ## begin the while loop, which goes until the log likelihood it calculates
  ## is very very close to the log likelihood at the last iteration (the
  ## difference is less than 10^(-5) different
  
  while (
    if((counter > 0) & (estconv == TRUE)){
      (max(abs(par - ppar)) > em.converge & (counter < maxIter))
    } else {
    (((llik - pllik) > em.converge) & (counter < maxIter))
    }) {
    
    # weight: *note that this is the weight for forced-known design
    fx <- logistic(x %*% bpar)  
    gx1 <- (o == 1)*logistic(xg1 %*% tpar) + (o == 0)*(1 - logistic(xg1 %*% tpar))
    gx0 <- (o == 1)*logistic(xg0 %*% tpar) + (o == 0)*(1 - logistic(xg0 %*% tpar))
    gxy <- (o == 1)*logistic(xgy %*% tpar) + (o == 0)*(1 - logistic(xgy %*% tpar))
    
    wr <- exp(log(p) + log(gxy) + y * log(fx) + (1-y) * log(1-fx)) / 
      (exp(log(p) + log(gxy) + y * log(fx) + (1-y) * log(1-fx)) +
       exp(y*log(p1) + (1-y) * log(p0) + log(gx1 * fx + gx0 * (1-fx))))
    
    ## M step to obtain updated coefficients for bpar and tpar
    
    if(fit.sens == "glm") {
        zfit <- glm(cbind(y, 1-y) ~ x - 1, family = binomial("logit"), weights = wr, 
                    control = list(maxit = glmMaxIter))
    } else if(fit.sens == "bayesglm") {
        zfit <- bayesglm(cbind(y, 1-y) ~ x - 1, family = binomial("logit"), weights = wr,
                         control = glm.control(maxit = maxIter), prior.scale.for.intercept = 2.5)
    }
    
    if(fit.outcome == "glm") {
        vfit <- glm(cbind(o, 1-o) ~ xgy - 1, family = binomial("logit"), weights = wr,
                    control = list(maxit = glmMaxIter)) 
    } else if(fit.outcome == "bayesglm") {
        vfit <- bayesglm(cbind(o, 1-o) ~ xgy - 1, family = binomial("logit"), weights = wr,
                         control = glm.control(maxit = maxIter), prior.scale.for.intercept = 2.5)
    }
    
    ppar <- par
    pllik <- llik
    
    bpar <- coef(zfit)
    tpar <- coef(vfit)
    par <- c(bpar, tpar)
    
    if(verbose) {
        cat(paste("llik", counter, round(llik, 4), "\n"))
        cat(paste("par", counter, max(abs(par-ppar)), "\n"))
    }
    
    ## calculate the new log likelihood with the parameters you just
    ## estimated in the M step

    llik <- obs.llik(par = par, y = y, x = x, o = o, p = p, p1 = p1, xg1 = xg1, xg0 = xg0, fit.sens = fit.sens, fit.outcome = fit.outcome)
    
    counter <- counter + 1
    
    # Warning if llik is not monotonically increasing
    if (llik < pllik) {
        if (verbose) {
            cat(paste("llik", counter, round(llik, 4), "\n"))
            cat(paste("par", counter, max(abs(par-ppar)), "\n"))
        }
        warning("log-likelihood did not always increase monotonically.")
    }
    
    if(counter == (maxIter-1))
      warning("number of iterations exceeded maximum in ML.")
    
  } ## end of while loop
  
  # Rescale coefficients
  bpar <- bpar/bscale
  tpar <- tpar/tscale
  
  # Now get standard errors
  if(optim) {
    MLEfit <- optim(par = par, fn = obs.llik, method = "BFGS",
                    y = y, x = x, o = o, xg1 = xg1, xg0 = xg0,
                    p = p, p1 = p1, fit.sens = fit.sens, fit.outcome = fit.outcome,
                    hessian = TRUE,
                    control = list(maxit = 0, fnscale = -1))
    
    vcov <- solve(-MLEfit$hessian)
    
    #rescale vcov by dividing the diag by var(x_i) and the off diag by sd(x_i)sd(x_j)
    vcov <- vcov/(scale %*% t(scale))
    
    se <- sqrt(diag(vcov))
    
    se.b <- se[1:ncol(x)]
    se.t <- se[ncol(x)+1:ncol(xg1)]
    
  } else {
    bpar.us <- par[1:ncol(x)] #take unscaled bpar and tpar
    tpar.us <- par[ncol(x)+1:ncol(xg1)]
    fx <- logistic(x %*% bpar.us)
    gx1 <- (o == 1)*logistic(xg1 %*% tpar.us) + (o == 0)*(1 - logistic(xg1 %*% tpar.us))
    gx0 <- (o == 1)*logistic(xg0 %*% tpar.us) + (o == 0)*(1 - logistic(xg0 %*% tpar.us))
    
    if((1 - p - p1) == 0){ ## when probability of forced "no" is 0
      forceno <- 1
    } else {
      forceno <- (1 - p - p1)
    }
    
    if(p1 == 0){ ## when probability of forced "yes" is 0
      forceyes <- 1
    } else {
      forceyes <- p1
    }
    
    score.rrpred <- cbind( (as.vector(fx * gx1 * (p + p1)^y * forceno^(1 - y)) + 
                                 as.vector((1 - fx) * gx0 * forceyes^y * (1 - p1)^(1 - y)))^(-1) *
                           ((as.vector(fx * gx1 * (p + p1)^y * forceno^(1 - y) * (1 - fx)) * x) +
                           (as.vector(-fx * gx0 * forceyes^y * (1 - p1)^(1 - y) * (1 - fx)) * x)), 
                           (as.vector(fx * gx1 * (p + p1)^y * forceno^(1 - y)) + 
                                 as.vector((1 - fx) * gx0 * forceyes^y * (1 - p1)^(1 - y)))^(-1) *
                           ((as.vector(fx * gx1 * (p + p1)^y * forceno^(1 - y) * (1 - gx1)) * xg1) +
                           (as.vector((1 - fx) * gx0 * forceyes^y * (1 - p1)^(1 - y) * (1 - gx0)) * xg0))
                                )
    
    information.rrpred <- crossprod(score.rrpred) / nrow(x)
    vcov <- solve(information.rrpred, tol = solve.tolerance) / nrow(x)
    vcov <- vcov/(scale %*% t(scale))
    
    se <- sqrt(diag(vcov))
    
    se.b <- se[1:ncol(x)]
    se.t <- se[ncol(x)+1:ncol(xg1)]
      
  }

  return.object.rrreg.predictor <- list(est.t = tpar,
                                  se.t = se.t,
                                  est.b = bpar,
                                  se.b = se.b,
                                  vcov = vcov,
                                  data = df,
                                  coef.names = c(colnames(x), paste(rr.item)),
                                  x = x1, #unscaled data
                                  y = y,
                                  o = o,
                                  design = design,
                                  call = match.call()                              
                                  )  
  
  if(design == "forced-known") {
    return.object.rrreg.predictor$p <- p
    return.object.rrreg.predictor$p1 <- p1
    return.object.rrreg.predictor$p0 <- p0
  }
  
  if(design == "mirrored" | design == "disguised") {
    return.object.rrreg.predictor$p <- p
  }
  
  if(design == "unrelated-known") {
    return.object.rrreg.predictor$p <- p
    return.object.rrreg.predictor$q <- q
  }
  
  class(return.object.rrreg.predictor) <- "rrreg.predictor"
  
  return(return.object.rrreg.predictor)
  
}


#' @export
vcov.rrreg.predictor <- function(object, ...){
  
  vcov <- object$vcov
  
  rownames(vcov) <- colnames(vcov) <- c(colnames(object$x), object$coef.names)
  
  return(vcov)
  
}

#' @export
coef.rrreg.predictor <- function(object, ...){
  
  coef <- object$est.t
  
  names(coef) <- object$coef.names
  
  return(coef)
  
}

#' @export
summary.rrreg.predictor <- function(object, ...) {
  structure(object, class = c("summary.rrreg.predictor", class(object)))
}

#' @export
print.rrreg.predictor <- function(x, ...){
  
  cat("\nRandomized Response as a Regression Predictor \n\nCall: ")
  
  dput(x$call)
  
  cat("\n")
  
  tb <- matrix(NA, ncol = 4, nrow = length(x$est.t))
  colnames(tb) <- c("est.t", "se.t", "est.b", "se.b")
  rownames(tb) <- x$coef.names
  
  tb[,1] <- x$est.t
  tb[,2] <- x$se.t
  tb[,3] <- c(x$est.b, NA)
  tb[,4] <- c(x$se.b, NA)
  
  print(as.matrix(round(tb,5)))
  
  cat("\n")
  
  summarize.design(x$design, x$p, x$p0, x$p1, x$q)
  
  cat("\n")
  
  invisible(x)
  
}

#' @export
print.summary.rrreg.predictor <- function(x, ...){
  
  cat("\nRandomized Response as a Regression Predictor \n\nCall: ")
  
  dput(x$call)
  
  cat("\n")
  
  tb <- matrix(NA, ncol = 4, nrow = length(x$est.t))
  colnames(tb) <- c("est.t", "se.t", "est.b", "se.b")
  rownames(tb) <- x$coef.names
  
  tb[,1] <- x$est.t
  tb[,2] <- x$se.t
  tb[,3] <- c(x$est.b, NA)
  tb[,4] <- c(x$se.b, NA)
  
  print(as.matrix(round(tb,5)))
  
  cat("\n")
  
  summarize.design(x$design, x$p, x$p0, x$p1, x$q)
  
  cat("\n")
  
  invisible(x)
  
}



#' Predicted Probabilities for Randomized Response as a Regression Predictor
#' 
#' \code{predict.rrreg.predictor} is used to generate predicted probabilities
#' from a multivariate regression object of survey data using the randomized
#' response item as a predictor for an additional outcome.
#' 
#' This function allows users to generate predicted probabilities for the
#' additional outcome variables with the randomized response item as a
#' covariate given an object of class "rrreg.predictor" from the
#' \code{rrreg.predictor()} function. Four standard designs are accepted by
#' this function: mirrored question, forced response, disguised response, and
#' unrelated question. The design, already specified in the "rrreg.predictor"
#' object, is then directly inputted into this function.
#' 
#' @usage predict.rrreg.predictor(object, fix.z = NULL, alpha = .05,
#' n.sims = 1000, avg = FALSE, newdata = NULL, quasi.bayes = FALSE, keep.draws
#' = FALSE, ...)
#' @param object An object of class "rrreg.predictor" generated by the
#' \code{rrreg.predictor()} function.
#' @param fix.z An optional value or vector of values between 0 and 1 that the
#' user inputs as the proportion of respondents with the sensitive trait or
#' probability that each respondent has the sensitive trait, respectively. If
#' the user inputs a vector of values, the vector must be the length of the
#' data from the "rrreg.predictor" object. Default is \code{NULL} in which case
#' predicted probabilities are generated for the randomized response item.
#' @param alpha Confidence level for the hypothesis test to generate upper and
#' lower confidence intervals. Default is \code{.05}.
#' @param n.sims Number of sampled draws for quasi-bayesian predicted
#' probability estimation. Default is \code{1000}.
#' @param avg Whether to output the mean of the predicted probabilities and
#' uncertainty estimates. Default is \code{FALSE}.
#' @param newdata Optional new data frame of covariates provided by the user.
#' Otherwise, the original data frame from the "rreg" object is used.
#' @param quasi.bayes Option to use Monte Carlo simulations to generate
#' uncertainty estimates for predicted probabilities. Default is \code{FALSE}
#' meaning no uncertainty estimates are outputted.
#' @param keep.draws Option to return the Monte Carlos draws of the quantity of
#' interest, for use in calculating differences for example.
#' @param ... Further arguments to be passed to
#' \code{predict.rrreg.predictor()} command.
#' @return \code{predict.rrreg.predictor} returns predicted probabilities
#' either for each observation in the data frame or the average over all
#' observations. The output is a list that contains the following components:
#' 
#' \item{est}{Predicted probabilities of the additional outcome variable given
#' the randomized response item as a predictor generated either using fitted
#' values or quasi-Bayesian simulations. If \code{avg} is set to \code{TRUE},
#' the output will only include the mean estimate.} \item{se}{Standard errors
#' for the predicted probabilities of the additional outcome variable given the
#' randomized response item as a predictor generated using Monte Carlo
#' simulations. If \code{quasi.bayes} is set to \code{FALSE}, no standard
#' errors will be outputted.} \item{ci.lower}{Estimates for the lower
#' confidence interval. If \code{quasi.bayes} is set to \code{FALSE}, no
#' confidence interval estimate will be outputted.} \item{ci.upper}{Estimates
#' for the upper confidence interval. If \code{quasi.bayes} is set to
#' \code{FALSE}, no confidence interval estimate will be outputted.}
#' \item{qoi.draws}{Monte Carlos draws of the quantity of interest, returned
#' only if \code{keep.draws} is set to \code{TRUE}.}
#' @seealso \code{\link{rrreg.predictor}} to conduct multivariate regression
#' analyses with the randomized response as predictor in order to generate
#' predicted probabilities.
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design
#' and Analysis of the Randomized Response Technique."  \emph{Working Paper.}
#' Available at \url{http://imai.princeton.edu/research/randresp.html}.
#' @keywords predicted probabilities fitted values
#' @examples
#' 
#' \dontrun{
#' data(nigeria)
#' 
#' ## Define design parameters
#' 
#' set.seed(44)
#' 
#' p <- 2/3  # probability of answering honestly in Forced Response Design
#' p1 <- 1/6 # probability of forced 'yes'
#' p0 <- 1/6 # probability of forced 'no'
#' 
#' ## Fit joint model of responses to an outcome regression of joining a civic 
#' ## group and the randomized response item of having a militant social connection
#' 
#' rr.q1.pred.obj <- 
#'     rrreg.predictor(civic ~ cov.asset.index + cov.married + I(cov.age/10) + 
#'               I((cov.age/10)^2) + cov.education + cov.female 
#'               + rr.q1, rr.item = "rr.q1", parstart = FALSE, estconv = TRUE,
#'               data = nigeria, verbose = FALSE, optim = TRUE,
#'               p = p, p1 = p1, p0 = p0, design = "forced-known")
#' 
#' ## Generate predicted probabilities for the likelihood of joining 
#' ## a civic group across respondents using quasi-Bayesian simulations. 
#' 
#' rr.q1.rrreg.predictor.pred <- predict(rr.q1.pred.obj, 
#'                                  avg = TRUE, quasi.bayes = TRUE, 
#'                                  n.sims = 10000)
#' 
#' }
#' @importFrom MASS mvrnorm
#' @method predict rrreg.predictor
predict.rrreg.predictor <- function(object, fix.z = NULL, alpha = .05, n.sims = 1000, 
                           avg = FALSE, newdata = NULL, quasi.bayes = FALSE, keep.draws = FALSE, ...) {
  
  if(missing(fix.z)){
  
  if(missing(newdata)) {
    #generate predicted probs for z
    propz <- (object$p * object$y + object$p0 * (1-object$y) + object$p1 * object$y) * logistic(object$x %*% object$est.b) /
      (object$p * (object$y * logistic(object$x %*% object$est.b) + (1-object$y) * (1-logistic(object$x %*% object$est.b))) +
         object$p1 * object$y + object$p0 * (1-object$y))
    
    xvar <- cbind(object$x, propz)
  } else {
    if(nrow(newdata)==0)
      stop("No data in the provided data frame.")
    
    #generate predicted probs for z
    covars <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata)
    
    newx <- covars[,!colnames(covars) %in% paste(object$call$rr.item)]
    
    newy <- covars[,colnames(covars) %in% paste(object$call$rr.item)]
    
    propz <- (object$p * newy + object$p0 * (1-newy) + object$p1 * newy) * logistic(newx %*% object$est.b) /
      (object$p * (newy * logistic(newx %*% object$est.b) + (1-newy) * (1-logistic(newx %*% object$est.b))) +
         object$p1 * newy + object$p0 * (1-newy))
    
    xvar <- cbind(newx, propz)
  }
  
  } else{
    if(length(fix.z) == 1){
      if(any((fix.z > 1) | (fix.z < 0))){
        stop("fix.z must be a probability (of respondents possessing the sensitive trait)")
      }
      
      propz <- rep(fix.z, length(object$y))
      
    } else if(length(fix.z) == length(object$y)){
      if(any((fix.z > 1) | (fix.z < 0))){
        stop("fix.z must be a probability (of respondents possessing the sensitive trait)")
      }
      propz <- fix.z
      
    } else{
      stop("fix.z must be either a single value or a vector of 
           length of the number of observations of the data from the 'rrreg.predictor' object")
    }
    
    if(missing(newdata)) {
      #generate predicted probs for z
      xvar <- cbind(object$x, propz)
    } else {
      if(nrow(newdata)==0)
        stop("No data in the provided data frame.")
      
      covars <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata)
      
      newx <- covars[,!colnames(covars) %in% paste(object$call$rr.item)]
      
      xvar <- cbind(newx, propz)
    }
  } 
  
  if(quasi.bayes == FALSE) {
    
      pred <- logistic(xvar %*% object$est.t)
    
    if(avg == TRUE)
      pred <- t(as.matrix(apply(pred, 2, mean)))
    
    return.object <- list(est = pred)
    
  } else {
    
    draws <- mvrnorm(n = n.sims, mu = c(object$est.b, object$est.t), Sigma = object$vcov)
    
      pred <- logistic(xvar %*% t(draws[,(length(object$est.b)+1):(length(c(object$est.b, object$est.t)))]))
    
    if(avg == TRUE)
      pred <- t(as.matrix(apply(pred, 2, mean)))
    
    est <- apply(pred, 1, mean)
    lower.ci <- apply(pred, 1, function(x, a = alpha) quantile(x, a/2))  
    upper.ci <- apply(pred, 1, function(x, a = alpha) quantile(x, 1-a/2)) 
    se <- apply(pred, 1, sd)
    
    return.object <- list(est = est, 
                          se = se, 
                          ci.lower = lower.ci, 
                          ci.upper = upper.ci)

    if(keep.draws == TRUE)
        return.object$qoi.draws <- pred
        
    }
    
    
   return(return.object) 
    
}
  

