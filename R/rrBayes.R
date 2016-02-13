logistic <- function(x) exp(x)/(1+exp(x))



#' Bayesian Randomized Response Regression
#' 
#' Function to conduct multivariate regression analyses of survey data with the
#' randomized response technique using Bayesian MCMC.
#' 
#' This function allows the user to perform regression analysis on data from
#' the randomized response technique using a Bayesian MCMC algorithm.
#' 
#' The Metropolis algorithm for the Bayesian MCMC estimators in this function
#' must be tuned to work correctly. The \code{beta.tune} and, for the mixed
#' effects model \code{Psi.tune}, are required, and the values, one for each
#' estimated parameter, will need to be manipulated. The output of the
#' \code{rrreg.bayes} function displays the acceptance ratios from the
#' Metropolis algorithm. If these values are far from 0.4, the tuning
#' parameters should be changed until the ratios approach 0.4.
#' 
#' Convergence is at times difficult to achieve, so we recommend running
#' multiple chains from overdispersed starting values by, for example, running
#' an MLE using the rrreg() function, and then generating a set of
#' overdispersed starting values using those estimates and their estimated
#' variance-covariance matrix. An example is provided below for each of the
#' possible designs. Running \code{summary()} after such a procedure will
#' output the Gelman-Rubin convergence statistics in addition to the estimates.
#' If the G-R statistics are all below 1.1, the model is said to have
#' converged.
#' 
#' @usage rrreg.bayes(formula, p, p0, p1, design, data, group.mixed,
#' formula.mixed = ~1, verbose = FALSE, n.draws = 10000, burnin = 5000, thin =
#' 1, beta.start, beta.mu0, beta.A0, beta.tune, Psi.start, Psi.df, Psi.scale,
#' Psi.tune)
#' @param formula An object of class "formula": a symbolic description of the
#' model to be fitted.
#' @param p The probability of receiving the sensitive question (Mirrored
#' Question Design, Unrelated Question Design); the probability of answering
#' truthfully (Forced Response Design); the probability of selecting a red card
#' from the 'yes' stack (Disguised Response Design).
#' @param p0 The probability of forced 'no' (Forced Response Design).
#' @param p1 The probability of forced 'yes' (Forced Response Design).
#' @param design Character indicating the design. Currently only "forced-known"
#' is supported.
#' @param data A data frame containing the variables in the model.
#' @param group.mixed A string indicating the variable name of a numerical
#' group indicator specifying which group each individual belongs to for a
#' mixed effects model.
#' @param formula.mixed To specify a mixed effects model, include this formula
#' object for the group-level fit. ~1 allows intercepts to vary, and including
#' covariates in the formula allows the slopes to vary also.
#' @param verbose A logical value indicating whether model diagnostics are
#' printed out during fitting.
#' @param n.draws Number of MCMC iterations.
#' @param burnin The number of initial MCMC iterations that are discarded.
#' @param thin The interval of thinning between consecutive retained iterations
#' (1 for no thinning).
#' @param beta.start Optional starting values for the sensitive item fit. This
#' should be a vector of length the number of covariates.
#' @param beta.mu0 Optional vector of prior means for the sensitive item fit
#' parameters, a vector of length the number of covariates.
#' @param beta.A0 Optional matrix of prior precisions for the sensitive item
#' fit parameters, a matrix of dimension the number of covariates.
#' @param beta.tune A required vector of tuning parameters for the Metropolis
#' algorithm for the sensitive item fit. This must be set and refined by the
#' user until the acceptance ratios are approximately .4 (reported in the
#' output).
#' @param Psi.start Optional starting values for the variance of the random
#' effects in the mixed effects models. This should be a scalar.
#' @param Psi.df Optional prior degrees of freedom parameter for the variance
#' of the random effects in the mixed effects models.
#' @param Psi.scale Optional prior scale parameter for the variance of the
#' random effects in the mixed effects models.
#' @param Psi.tune A required vector of tuning parameters for the Metropolis
#' algorithm for variance of the random effects in the mixed effects models.
#' This must be set and refined by the user until the acceptance ratios are
#' approximately .4 (reported in the output).
#' @return \code{rrreg.bayes} returns an object of class "rrreg.bayes".  The
#' function \code{summary} is used to obtain a table of the results.
#' 
#' \item{beta}{The coefficients for the sensitive item fit. An object of class
#' "mcmc" that can be analyzed using the \code{coda} package.} \item{data}{The
#' \code{data} argument.} \item{coef.names}{Variable names as defined in the
#' data frame.} \item{x}{The model matrix of covariates.} \item{y}{The
#' randomized response vector.} \item{design}{Call of standard design used:
#' "forced-known", "mirrored", "disguised", or "unrelated-known".} \item{p}{The
#' \code{p} argument.} \item{p0}{The \code{p0} argument.} \item{p1}{The
#' \code{p1} argument.} \item{mixed}{Indicator for whether a mixed effects
#' model was run.} \item{call}{the matched call.}
#' 
#' If a mixed-effects model is used, then several additional objects are
#' included: \item{Psi}{The coefficients for the group-level fit. An object of
#' class "mcmc" that can be analyzed using the \code{coda} package.}
#' \item{gamma}{The random effects estimates. An object of class "mcmc" that
#' can be analyzed using the \code{coda} package.}
#' \item{coef.names.mixed}{Variable names for the predictors for the
#' second-level model} \item{z}{The predictors for the second-level model.}
#' \item{groups}{A vector of group indicators.}
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design
#' and Analysis of the Randomized Response Technique."  \emph{Working Paper.}
#' Available at \url{http://imai.princeton.edu/research/randresp.html}.
#' @examples
#' 
#'  \dontrun{ 
#'  
#' data(nigeria)
#'  
#' ## Define design parameters
#' p <- 2/3  # probability of answering honestly in Forced Response Design
#' p1 <- 1/6 # probability of forced 'yes'
#' p0 <- 1/6 # probability of forced 'no'
#' 
#' ## run three chains with overdispersed starting values
#' 
#' set.seed(1)
#' 
#' ## starting values constructed from MLE model
#' mle.estimates <- rrreg(rr.q1 ~ cov.asset.index + cov.married + 
#'                          I(cov.age/10) + I((cov.age/10)^2) + cov.education + cov.female, 
#'                          data = nigeria, 
#'                       p = p, p1 = p1, p0 = p0,
#'                       design = "forced-known")
#' 
#' library(MASS)
#' draws <- mvrnorm(n = 3, mu = coef(mle.estimates), 
#'   Sigma = vcov(mle.estimates) * 9)
#' 
#' ## run three chains
#' bayes.1 <- rrreg.bayes(rr.q1 ~ cov.asset.index + cov.married + 
#'                          I(cov.age/10) + I((cov.age/10)^2) + cov.education + cov.female,   
#'                       data = nigeria, p = p, p1 = p1, p0 = p0,
#'                       beta.tune = .0001, beta.start = draws[1,],
#'                       design = "forced-known")
#' 
#' bayes.2 <- rrreg.bayes(rr.q1 ~ cov.asset.index + cov.married + 
#'                          I(cov.age/10) + I((cov.age/10)^2) + cov.education + cov.female,   
#'                       data = nigeria, p = p, p1 = p1, p0 = p0,
#'                       beta.tune = .0001, beta.start = draws[2,],
#'                       design = "forced-known")
#' 
#' bayes.3 <- rrreg.bayes(rr.q1 ~ cov.asset.index + cov.married + 
#'                          I(cov.age/10) + I((cov.age/10)^2) + cov.education + cov.female,   
#'                       data = nigeria, p = p, p1 = p1, p0 = p0,
#'                       beta.tune = .0001, beta.start = draws[3,],
#'                       design = "forced-known")
#'                       
#' bayes <- as.list(bayes.1, bayes.2, bayes.3)
#' 
#' summary(bayes)
#' 
#'  }
#' @importFrom coda mcmc
#' 
#' @export
rrreg.bayes <- function(formula, p, p0, p1, design, data, 
                        group.mixed, formula.mixed = ~1,
                        verbose = FALSE, n.draws = 10000,
                        burnin = 5000, thin = 1,
                        beta.start, beta.mu0, beta.A0, beta.tune,
                        Psi.start, Psi.df, Psi.scale, Psi.tune) {
  
  ictreg.call <- match.call()
  
  # set up data frame, with support for standard and modified responses
  mf <- match.call(expand.dots = FALSE)
  
  # make all other call elements null in mf <- NULL in next line
  mf$p <- mf$p0 <- mf$p1 <- mf$q <- mf$design <- mf$group.mixed <- mf$formula.mixed <- mf$start <- mf$maxIter <- mf$verbose <- mf$optim <- mf$em.converge <- mf$n.draws <- mf$beta.start <- mf$beta.mu0 <- mf$beta.A0 <- mf$beta.tune <- mf$Psi.start <- mf$Psi.df <- mf$Psi.scale <- mf$Psi.tune <- mf$burnin <- mf$thin <- NULL
    
  mf[[1]] <- as.name("model.frame")
  mf$na.action <- 'na.pass'
  mf <- eval.parent(mf)
  
  # check bayes inputs
  if(n.draws - burnin < 0)
    stop("The input burnin must be set to a number smaller than n.draws.")
  burnin <- burnin + 1
  if(burnin < 1)
    stop("The input burnin must be greater than or equal to zero.")
  if(thin > (n.draws - burnin))
    stop("The input thin must be greater than the number of draws after burnin (n.draws - burnin).")
  
  # define design, response data frames
  x <- model.matrix.default(attr(mf, "terms"), mf)
  y <- model.response(mf)
  
  # get mixed effects group-level predictors
  mixed <- missing("group.mixed") == FALSE
    
  if (mixed == TRUE)
    z <- model.matrix(formula.mixed, model.frame(formula.mixed, data, na.action = 'na.pass'))
  
  # list-wise missing deletion
  na.x <- apply(is.na(x), 1, sum)
  na.y <- is.na(y)
  
  if (mixed == TRUE) {
    na.z <- apply(is.na(z), 1, sum)
    na.cond <- na.x==0 & na.y==0 & na.z==0
  } else {
    na.cond <- na.x==0 & na.y==0
  }
    
  ## list wise delete
  y <- y[na.cond == TRUE]
  x <- x[na.cond == TRUE, , drop = FALSE]
  if (mixed == TRUE)
    z <- z[na.cond == TRUE, , drop = FALSE]
    
  ## group indicator for mixed effects regression
  if (mixed == TRUE) {
    if(!(group.mixed %in% colnames(data)))
      stop("The covariate named in group.mixed cannot be found in the provided data frame.")
    grp <- data[na.cond == TRUE, paste(group.mixed)]
    if(class(grp) == "character")
      grp <- as.factor(grp)
    
    grp.labels <- sort(unique(grp))
    z.grp <- unique(cbind(grp, z))
    if(nrow(z.grp) > length(unique(grp)))
      warning("Some or all of the variables specified in formula.mixed are not constant within groups.")
  }
    
  n <- length(y)
  
  ## set up starting values
  
  nPar <- ncol(x)
  
  if (missing("beta.start")) {
    beta.start <- runif(nPar)
  }
  if (missing("beta.mu0")) {
    beta.mu0 <- rep(0, nPar)
  }
  
  if (missing("beta.A0")) {
    beta.A0 <- diag(nPar) * .01
  }
  
  if (missing("beta.tune")) {
    stop("The Metropolis tuning input object beta.tune is required.")
  } else if (class(beta.tune) == "numeric" & length(beta.tune) == 1 & nPar > 1) {
    beta.tune <- diag(nPar) * beta.tune
  } else if (class(beta.tune) == "numeric" & length(beta.tune) == nPar) {
    beta.tune <- diag(beta.tune)
  } else if (class(beta.tune) == "numeric" | (class(beta.tune) == "matrix" & 
                                                (ncol(beta.tune) != nPar | nrow(beta.tune) != nPar))) {
    stop("The input beta.tune is malformed. It should be a scalar, a vector of length the number of predictors, or a diagonal matrix of dimensions the number of predictors.")
  }
    
  if (mixed == TRUE) {
            
    if (missing("Psi.start")) {
      Psi.start <- diag(10, ncol(z)) 
    } else if (class(Psi.start) == "numeric" & 
                 (length(Psi.start) == 1 | length(Psi.start) == ncol(z))) {
      Psi.start <- diag(ncol(z)) * Psi.start
    } else if (class(Psi.start) == "numeric" | (class(Psi.start) == "matrix" & 
                                                  (nrow(Psi.start) != ncol(z) | ncol(Psi.start) != ncol(z)))) {
      stop("The input Psi.start is malformed. It should either be a square matrix of dimensions the number of columns of Z, a scalar, or a vector of length number of columns of Z.")
    }
    
    if (missing("Psi.df")) {
      Psi.df <- ncol(z) + 2
    }
    if (missing("Psi.scale")) {
      Psi.scale <- diag(2, ncol(z)) 
    } else if (class(Psi.scale) == "numeric" & 
                 (length(Psi.scale) == 1 | length(Psi.scale) == ncol(z))) {
      Psi.scale <- diag(ncol(z)) * Psi.scale
    } else if (class(Psi.scale) == "numeric" | (class(Psi.scale) == "matrix" & 
                 (nrow(Psi.scale) != ncol(z) | ncol(Psi.scale) != ncol(z)))) {
      stop("The input Psi.scale is malformed. It should either be a square matrix of dimensions the number of columns of Z, a scalar, or a vector of length number of columns of Z.")
    }

    if (missing("Psi.tune")) {
      stop("The Metropolis tuning input object psi.tune is required.")
    } else if (class(Psi.tune) == "numeric" & length(Psi.tune) == 1 & length(unique(grp)) > 1) {
      Psi.tune <- rep(Psi.tune, length(unique(grp)))
    } else if (class(Psi.tune) == "matrix" | (class(Psi.tune) == "numeric" & (length(Psi.tune) != length(unique(grp))))) {
      stop("The input Psi.tune is malformed. It should be either a scalar or a vector of length the number of groups.")
    }
  }
  
  if(mixed == FALSE) {
    bayes.fit <- rrbayes.fit(Y = y, X = x, p = p, p1 = p1, n.draws = n.draws,
                             beta.start = beta.start, beta.mu0 = beta.mu0, 
                             beta.A0 = beta.A0, beta.prop = beta.tune, verbose = verbose)
    
    beta.mcmc <- mcmc(data = bayes.fit, start = burnin, thin = thin, end = nrow(bayes.fit))
    
    return.object <- list(beta = beta.mcmc,
                          data = df,
                          coef.names = colnames(x),
                          x = x,
                          y = y,
                          design = design,
                          p = p,
                          p1 = p1,
                          p0 = p0,
                          mixed = mixed,
                          call = match.call())
    
  } else {
    bayes.fit <- rrbayesmixed.fit(Y = y, X = x, Z = z, 
                                  p = p, p1 = p1, grp = as.numeric(grp),
                                  n.draws = n.draws, 
                                  beta.start = beta.start, Psi = Psi.start,
                                  beta.mu0 = beta.mu0, beta.A0 = beta.A0, Psi.df = Psi.df,
                                  Psi.scale = Psi.scale, 
                                  beta.tune = beta.tune, Psi.tune = Psi.tune, verbose = verbose)
        
    beta.mcmc <- mcmc(data = bayes.fit$beta, start = burnin, 
                      thin = thin, end = nrow(bayes.fit$beta))
    gamma.mcmc <- mcmc(data = bayes.fit$gamma, start = burnin, 
                       thin = thin, end = nrow(bayes.fit$gamma))
    Psi.mcmc <- mcmc(data = bayes.fit$Psi, start = burnin, 
                     thin = thin, end = nrow(bayes.fit$Psi))
    
    return.object <- list(beta = beta.mcmc, 
                          gamma = gamma.mcmc, 
                          Psi = Psi.mcmc,
                          data = df,
                          coef.names = colnames(x),
                          coef.names.mixed = colnames(z),
                          x = x,
                          y = y,
                          z = z,
                          groups = grp,
                          group.names = grp.labels,
                          design = design,
                          p = p,
                          p1 = p1,
                          p0 = p0,
                          mixed = mixed,
                          call = match.call())
    
  }
  
  class(return.object) <- "rrreg.bayes"
  
  return(return.object)
  
}

#' @export
coef.rrreg.bayes <- function(object, ranef = FALSE, ...) {
  
  beta.coef <- apply(object$beta, 2, mean)
  names(beta.coef) <- object$coef.names
  
  return.object <- list(beta = beta.coef)
  
  if(object$mixed == TRUE) {

    if (ranef == TRUE) {
      
      psi.label <- c()
      for(i in 1:length(object$coef.names.mixed)){
        for(j in 1:i){
          psi.label <- c(psi.label, paste(object$coef.names.mixed[i], object$coef.names.mixed[j]))
        }
      }
      
      Psi.coef <- apply(object$Psi, 2, mean)
      names(Psi.coef) <- psi.label
      
      return.object$Psi <- Psi.coef
      
      if(length(object$coef.names.mixed) > 1)
        gamma.label <- paste(rep(object$coef.names.mixed, length(object$group.names)), "group",
                                 rep(object$group.names, each = length(object$coef.names.mixed)))
      else
        gamma.label <- object$group.names
      
      gamma.coef <- apply(object$gamma, 2, mean)
      names(gamma.coef) <- gamma.label
      
      return.object$gamma <- gamma.coef
    }
     
  }
    
  return.object
  
}

#' @importFrom coda as.mcmc
#' @export
coef.rrreg.bayes.list <- function(object, ranef = FALSE, ...) {
  
  object$beta <- as.mcmc(do.call(rbind, as.list(object$beta)))
  
  if(object$mixed == TRUE){
    object$Psi <- as.mcmc(do.call(rbind, as.list(object$Psi)))
  
    if (ranef == TRUE)
      object$gamma <- as.mcmc(do.call(rbind, as.list(object$gamma)))
  }
  
  class(object) <- "rrreg.bayes"
  
  coef(object, ranef = ranef, ... = ...)
  
}


#' @method sd rrreg.bayes
sd.rrreg.bayes <- function(object, ranef = FALSE, ...) {
  
  beta.coef <- apply(object$beta, 2, sd)
  names(beta.coef) <- object$coef.names
  
  return.object <- list(beta = beta.coef)
  
  if(object$mixed == TRUE) {
    
    psi.label <- c()
    for(i in 1:length(object$coef.names.mixed)){
      for(j in 1:i){
        psi.label <- c(psi.label, paste(object$coef.names.mixed[i], object$coef.names.mixed[j]))
      }
    }
    
    Psi.coef <- apply(object$Psi, 2, sd)
    names(Psi.coef) <- psi.label
    
    return.object$Psi <- Psi.coef
    
    if (ranef == TRUE) {
      if(length(object$coef.names.mixed) > 1)
        gamma.label <- paste(rep(object$coef.names.mixed, length(object$group.names)), "group",
                             rep(object$group.names, each = length(object$coef.names.mixed)))
      else
        gamma.label <- object$group.names
      
      gamma.coef <- apply(object$gamma, 2, sd)
      names(gamma.coef) <- gamma.label
      
      return.object$gamma <- gamma.coef
    }
  }
  
  return.object
  
}

#' @method sd rrreg.bayes.list
sd.rrreg.bayes.list <- function(object, ranef = FALSE, ...) {
  
  object$beta <- as.mcmc(do.call(rbind, as.list(object$beta)))
  
  if(object$mixed == TRUE) {
    object$Psi <- as.mcmc(do.call(rbind, as.list(object$Psi)))
    
    if(ranef == TRUE)
      object$gamma <- as.mcmc(do.call(rbind, as.list(object$gamma)))
  }
  
  class(object) <- "rrreg.bayes"
  
  sd.rrreg.bayes(object, ranef = ranef, ... = ...)
  
}

#' @export
vcov.rrreg.bayes <- function(object, ...) {  
  if(object$mixed == FALSE) {
    cov(cbind(object$beta))
  } else {
    cov(cbind(object$beta, object$Psi, object$gamma))
  }
}

#' @export
vcov.rrreg.bayes.list <- function(object, ...) {
  
  object$beta <- as.mcmc(do.call(rbind, as.list(object$beta)))
  
  if(object$mixed == TRUE) {
    object$Psi <- as.mcmc(do.call(rbind, as.list(object$Psi)))
    object$gamma <- as.mcmc(do.call(rbind, as.list(object$gamma)))
  }
  class(object) <- "rrreg.bayes"
    
  vcov(object, ... = ...)
  
}

#' @export
print.rrreg.bayes <- print.rrreg.bayes.list <- function(x, ...) {
  
  cat("\nRandomized Response Technique Bayesian Regression \n\nCall: ")
  
  dput(x$call)
  
  cat("\nCoefficient estimates\n")
  
  print(coef(x))
    
  summarize.design(x$design, x$p, x$p0, x$p1, x$q)
  
  cat("\n\n")
  
  invisible(x)
  
}

#' @importFrom coda as.mcmc.list
#' @export
as.list.rrreg.bayes <- function(...) {
  
  x <- list(...)
  
  beta.list <- list()
  for (i in 1:length(x))
    beta.list[[i]] <- x[[i]]$beta
  
  beta.list <- as.mcmc.list(beta.list)
  
  if (x[[1]]$mixed == TRUE) {
    gamma.list <- list()
    for (i in 1:length(x))
      gamma.list[[i]] <- x[[i]]$gamma
    
    gamma.list <- as.mcmc.list(gamma.list)
    
    Psi.list <- list()
    for (i in 1:length(x))
      Psi.list[[i]] <- x[[i]]$Psi
    
    Psi.list <- as.mcmc.list(Psi.list)
  }
  
  return.object <- x[[1]]
  return.object$beta <- beta.list
  if (x[[1]]$mixed == TRUE) {
    return.object$gamma <- gamma.list
    return.object$Psi <- Psi.list
  }
  
  class(return.object) <- "rrreg.bayes.list"
  
  return.object
  
}

#' @export
summary.rrreg.bayes <- function(object, ...) {
  structure(object, class = c("summary.rrreg.bayes", class(object)))
}

#' @export
print.summary.rrreg.bayes <- function(x, ...) {
  
  cat("\nRandomized Response Bayesian Regression \n\nCall: ")
  
  dput(x$call)
  
  cat("\nIndividual-level predictors \n")
  print(matrix(c(round(cbind(coef(x)$beta, sd.rrreg.bayes(x)$beta),5)), nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
               dimnames = list(x$coef.names, c("Est.", "S.E."))))
    
  if(x$mixed == TRUE) {
    cat("\nRandom effects \n")
    
    if(length(x$coef.names.mixed) > 1)
      gamma.label <- paste(rep(x$coef.names.mixed, length(x$group.names)), "group",
                           rep(x$group.names, each = length(x$coef.names.mixed)))
    else
      gamma.label <- x$group.names
    
    print(matrix(c(round(cbind(coef(x, ranef = T)$gamma, 
                               sd.rrreg.bayes(x, ranef = T)$gamma),5)), 
                 nrow = length(gamma.label), ncol = 2, byrow = FALSE,
                 dimnames = list(gamma.label, c("Est.", "S.E."))))
    
    cat("\nVariance of the random effects\n")
    
    if(length(coef(x)$Psi) == 1) {
      psi.label <- "Psi"
    } else {
      psi.label <- c()
      for(i in 1:length(x$coef.names.mixed)){
        for(j in 1:i){
          psi.label <- c(psi.label, paste(x$coef.names.mixed[i], x$coef.names.mixed[j]))
        }
      }
    }
      
    print(matrix(c(round(cbind(coef(x)$Psi, sd.rrreg.bayes(x)$Psi),5)),
                 nrow = length(psi.label), ncol = 2, byrow = FALSE,
                 dimnames = list(psi.label, c("Est.", "S.E."))))
    
  }
  
  cat("\n")
  
  summarize.design(x$design, x$p, x$p0, x$p1, x$q)
  
  cat("\n\n")
  
  invisible(x)
  
}

#' @export
summary.rrreg.bayes.list <- function(object, ...) {
  structure(object, class = c("summary.rrreg.bayes.list", class(object)))
}


#' @importFrom coda gelman.diag
#' @export
print.summary.rrreg.bayes.list <- function(x, ...) {
  
  cat("\nRandomized Response Bayesian Regression \n\nCall: ")
  
  dput(x$call)
  
  cat("\nSummary from",length(x$beta),"chains")
  
  cat("\n\nIndividual-level predictors\n")
  
  print(matrix(c(round(cbind(coef(x)$beta, sd.rrreg.bayes.list(x)$beta),5)),
               nrow = length(x$coef.names), ncol = 2, byrow = FALSE,
               dimnames = list(x$coef.names, c("Est.", "S.E."))))
  
  cat("\nGelman-Rubin statistics:\n")
  
  gelmanrubin <- round(gelman.diag(x$beta)$psrf[,1],4)
  names(gelmanrubin) <- x$coef.names
  
  print(gelmanrubin)
  
  if(x$mixed == TRUE) {
    
    cat("\n\nRandom effects\n")
  
    if(length(x$coef.names.mixed) > 1)
      gamma.label <- paste(rep(x$coef.names.mixed, length(x$group.names), "group",
                               rep(x$group.names), each = length(x$coef.names.mixed)))
    else
      gamma.label <- x$group.names
    
    print(matrix(c(round(cbind(coef(x, ranef = T)$gamma, 
                               sd.rrreg.bayes.list(x, ranef = T)$gamma),5)), 
                 nrow = length(gamma.label), ncol = 2, byrow = FALSE,
                 dimnames = list(gamma.label, c("Est.", "S.E."))))
    
    cat("\nGelman-Rubin statistics:\n")
    
    gelmanrubin <- round(gelman.diag(x$gamma)$psrf[,1],4)
    names(gelmanrubin) <- unique(x$groups)
    
    print(gelmanrubin)
    
    cat("\nVariance of the random effects\n")
    
    if(length(coef(x)$Psi) == 1) {
      psi.label <- "Psi"
    } else {
      psi.label <- c()
      for(i in 1:length(x$coef.names.mixed)){
        for(j in 1:i){
          psi.label <- c(psi.label, paste(x$coef.names.mixed[i], x$coef.names.mixed[j]))
        }
      }
    }
    
    print(matrix(c(round(cbind(coef(x)$Psi, sd.rrreg.bayes.list(x)$Psi),5)),
                 nrow = length(psi.label), ncol = 2, byrow = FALSE,
                 dimnames = list(psi.label, c("Est.", "S.E."))))
    
    cat("\nGelman-Rubin statistics:\n")
    
    gelmanrubin <- round(gelman.diag(x$Psi)$psrf[,1],4)
    names(gelmanrubin) <- x$coef.names.mixed
    
    print(gelmanrubin)
    
  }
  
  cat("\n")
  
  summarize.design(x$design, x$p, x$p0, x$p1, x$q)
  
  cat("\n\n")
  
  invisible(x)
  
}


rrbayes.fit <- function(Y, X, p, p1, n.draws, beta.start, beta.mu0,
                        beta.A0, beta.prop, verbose) {
  
  n <- length(Y)
  k <- ncol(X)
  res <- .C("R2rrLogit", as.integer(Y), as.double(X),
            as.double(beta.start), as.double(p), as.double(p1),
            as.integer(n), as.integer(k), as.double(beta.mu0),
            as.double(beta.A0), as.double(beta.prop),
            as.integer(n.draws), as.integer(1),
            verbose = verbose,
            store = double(k*n.draws),
            PACKAGE = "rr")$store
  
  res <- matrix(res, byrow = TRUE, ncol = k)
  
  class(res) <- "rrBayes"
  return(res)
  
}


rrbayesmixed.fit <- function(Y, X, Z, p, p1, grp,
                             n.draws, beta.start, Psi,
                             beta.mu0, beta.A0, Psi.df,
                             Psi.scale, beta.tune, Psi.tune, verbose) {
  
  n <- length(Y)
  k <- ncol(X)
  m <- ncol(Z)
  n.grp <- length(table(grp))
  
  res <- .C("R2rrLogitMixed", as.integer(Y), as.double(X), as.double(t(Z)),
            as.double(p), as.double(p1), as.integer(grp-1), as.double(beta.start),
            as.double(Psi), as.integer(n), as.integer(k), as.integer(m),
            as.integer(n.grp), as.integer(max(table(grp))), as.double(beta.mu0),
            as.double(beta.A0), as.integer(Psi.df), as.double(Psi.scale),
            as.double(beta.tune), as.double(Psi.tune), as.integer(n.draws),
            integer(1), integer(n.grp),
            verbose = verbose,
            betaStore = double(n.draws*k), gammaStore = double(n.draws*m*n.grp),
            PsiStore = double(n.draws*m*(m+1)/2),
            PACKAGE = "rr")
  
  res <- list(beta = matrix(res$betaStore, byrow = TRUE, ncol = k),
              gamma = matrix(res$gammaStore, byrow = TRUE, ncol = m*n.grp),
              Psi = matrix(res$PsiStore, byrow = TRUE, ncol = m*(m+1)/2))
  return(res)
  
}
