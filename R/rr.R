logistic <- function(x) exp(x)/(1+exp(x))

rrreg <- function(formula, p, p0, p1, q, design, data, start = NULL, 
                  maxIter = 10000, verbose = FALSE, 
                  optim = FALSE, em.converge = 10^(-8), 
                  glmMaxIter = 10000, solve.tolerance = .Machine$double.eps) {
  
    df <- model.frame(formula, data, na.action = na.omit)
    x1 <- model.matrix.default(formula, df)
    bscale <- c(1, unname(apply((as.matrix(x1[,-1])), 2, sd)))
    x <- sweep(x1, 2, bscale, `/`) #rescale/standardize each covar, not including intercept
    y <- model.response(df)
    n <- length(y)

    # Get c, d
    cd <- rrcd(p, p0, p1, q, design)
    c <- cd$c
    d <- cd$d

    if(missing(design)) {
        stop("Missing design specification, see documentation")
    } else if(design != "forced-known" & design != "mirrored" & design != "disguised" & design != "unrelated-known"
              ) {
        stop(paste("The design you specified,", design, "is not implemented in this version of the software."))
    }

    obs.llik <- function(par, y, x, c, d, gx = NULL, bayesglm = FALSE) {
        if (is.null(gx))
            gx <- logistic(x %*% par)
        llik <- sum( y * log( c * gx + d) + (1-y) * log(1 - c * gx - d))
        ##if(bayesglm = TRUE)
        ##    llik <- llik + sum(dcauchy(x = par, scale = rep(2.5, length(par), log = TRUE))
        return(llik)
    }
    
    ## NOW THE ALGORITHM

    ## set some starting values
    if(!missing(start) & (length(start) != ncol(x)))
        stop("The starting values are not the right dimensions. start should be a vector of length the number of variables, including the intercept if applicable.")
        
    if(!is.null(start))
        par <- start
    else
        par <- rnorm(ncol(x))
    gx <- logistic(x %*% par)

    ## start off with an infinitely negative log likelihood
    pllik <- -Inf

    ## calculate the log likelihood at the starting values
    llik <- obs.llik(par, y = y, x = x, c = c, d = d, gx = gx)

    ## set a counter to zero, which you will iterate
    ## (this just allows to you stop after a maximum # of iterations if you want
    ## e.g. 10k

    counter <- 0

    ## begin the while loop, which goes until the log likelihood it calculates
    ## is very very close to the log likelihood at the last iteration (the
    ## difference is less than 10^(-8) different

    while (((llik - pllik) > em.converge) & (counter < maxIter)) {

        ## first the E step is run each time, from the parameters at the end of the
        ## last iteration, or the starting values in the first iteration
        w1 <- (c*gx) /  (c*gx + d)
        w0 <- c*(1-gx) / (1-c*gx - d)

        w <- rep(NA, n)
        w[y == 0] <- w0[y==0]
        w[y == 1] <- w1[y==1]

        ##if(bayesglm == FALSE) {
        lfit <- glm(cbind(y, 1-y) ~ x - 1, family = binomial("logit"), weights = w, start = par,
                    control = list(maxit = glmMaxIter))
        ##} else {
        ##    lfit <- bayesglm(as.formula(paste("cbind(", y.var, ", 1-", y.var, ") ~ ",
        ##                                           paste(x.vars.ceiling, collapse=" + "))),
        ##                          weights = dtmpC$w, family = binomial(logit),
        ##                          start = coef.qufit.start, data = dtmpC,
        ##                          control = glm.control(maxit = maxIter), scaled = F)
        ##}
        par <- coef(lfit)
        ## update gx
        gx <- logistic(x %*% par)

        pllik <- llik

        if(verbose==T)
            cat(paste(counter, round(llik, 4), "\n"))

        ## calculate the new log likelihood with the parameters you just
        ## estimated in the M step
        llik <- obs.llik(par, y = y, x = x, c = c, d = d, gx = gx)

        counter <- counter + 1

        ## stop for obvious reasons
        if (llik < pllik)
            stop("log-likelihood is not monotonically increasing.")

        if(counter == (maxIter-1))
            warning("number of iterations exceeded maximum in ML")

    } ## end of while loop
    
    # Rescale coefficients
    bpar <- c(par/bscale)
    
        if(optim) {
        MLEfit <- optim(par = par, fn = obs.llik, method = "BFGS",
                        y = y, x = x,
                        hessian = TRUE,
                        control = list(maxit = 0, fnscale = -1),
                        c = c, d = d)
        vcov <- solve(-MLEfit$hessian, tol = solve.tolerance)
        vcov <- vcov/(bscale %*% t(bscale))
        
        se <- sqrt(diag(vcov))
    } else { 
         score.rr <- as.vector(y * ( c * gx + d )^(-1) * c * gx * (1-gx)) * x +
           as.vector((1-y) * ( 1 - c * gx - d )^(-1) * (-c) * gx * (1-gx)) * x

         information.rr <- crossprod(score.rr) / nrow(x)
         vcov <- solve(information.rr, tol = solve.tolerance) / nrow(x)
         vcov <- vcov/(bscale %*% t(bscale))
         
         se <- sqrt(diag(vcov))
    }

    return.object <- list(est = bpar, 
                          vcov = vcov,
                          se = se,
                          data = df,
                          coef.names = colnames(x),
                          x = x1, #unscaled data
                          y = y,
                          design = design,
                          call = match.call())

    if(design == "forced-known") {
      return.object$p <- p
      return.object$p1 <- p1
      return.object$p0 <- p0
    }
    
    if(design == "mirrored" | design == "disguised") {
      return.object$p <- p
    }
    
    if(design == "unrelated-known") {
      return.object$p <- p
      return.object$q <- q
    }
    
    class(return.object) <- "rrreg"

    return(return.object)

}

summarize.design <- function(design, p, p0, p1, q) {
    
    design.short <- strsplit(design, "-")[[1]][1]
    if(design.short == "forced")
        param <- paste("p = ", round(p,2), ", p0 = ", round(p0,2), ", and p1 = ", round(p1, 2), sep = "")
    
    return(cat("Randomized response ", design.short, " design with ", param, ".", sep = ""))
    
}
    

rrcd <- function(p, p0, p1, q, design) {
   if(missing(design)) {
        stop("Missing design specification, see documentation")
    } else if(design != "forced-known" & design != "mirrored" & design != "disguised" & design != "unrelated-known") {
        stop(paste("The design you specified,", design, "is not implemented in this version of the software."))
    }
  
  if (design == "forced-known") {
    if(missing(p) & !missing(p1))
      stop("for the forced design with known probabilities, you must specify p")
    if(missing(p1) & !missing(p))
      stop("for the forced design with known probabilities, you must specify p1")
    if(missing(p) & missing(p1))
      stop("for the forced design with known probabilities, you must specify p and p1")
    
    if (any((p > 1) | (p < 0))) {
      stop("p must be a probability")
    }
    if (any((p1 > 1) | (p1 < 0))) {
      stop("p1 must be a probability")
    }
    if (any((p0 > 1) | (p0 < 0))) {
      stop("p0 must be a probability")
    }
    if (any( round(p + p0 + p1, 10) !=1)) {
      stop("p, p0, and p1 must sum to 1")
    }
    
    c <- p
    d <- p1
    
  } else if (design == "mirrored" | design == "disguised") {
   
    if(missing(p))
      stop(paste("for the", design, "design, you must specify p"))
    if (any((p > 1) | (p < 0)) | p == .5) {
      stop("p must be a probability that is not .5")
    }
    
    c <- 2 * p - 1
    d <- 1 - p
    
  } else if (design == "unrelated-known") {
    
    warning("q, the probability of answering 'yes' to the unrelated question, is assumed to be independent of covariates")
    
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
    
    c <- p
    d <- (1 - p) * q
  
  } 
  
  return(list(c = c,
              d = d))
}


vcov.rrreg <- function(object, ...){

  vcov <- object$vcov

  rownames(vcov) <- colnames(vcov) <- object$coef.names

  return(vcov)

}

coef.rrreg <- function(object, ...){

  coef <- object$est

  names(coef) <- object$coef.names

  return(coef)

}

predict.rrreg <- function(object, given.y = FALSE, alpha = .05, 
                          n.sims = 1000, avg = FALSE, newdata = NULL, 
                          quasi.bayes = FALSE, keep.draws = FALSE, ...) {
  
  z.post <- function(object = NULL, x = NULL, y = NULL, coef = NULL) {
    
    if(is.null(coef))
      coef <- object$est
    
    if(object$design == "forced-known") {

        pred <- ((object$p + object$p1)^y * object$p0^(1 - y) * logistic(x %*% coef)) /
            (object$p * logistic(x %*% coef)^y * (1 - logistic(x %*% coef))^(1 - y) +
                 object$p1^y * object$p0^(1 - y))
      
    } else if (object$design == "disguised") {

        pred <- ( logistic(x %*% coef) * object$p^y * (1 - object$p)^(1 - y) ) /
            ( logistic(x %*% coef) * object$p^y * (1 - object$p)^(1 - y) +
                 (1 - logistic(x %*% coef)) * (1 - object$p)^y * object$p^(1 - y) )

    } else if (object$design == "mirrored") {

        pred <- ( object$p^y * (1 - object$p)^(1 - y) * logistic(x %*% coef) ) /
            ( object$p^y * (1 - object$p)^(1 - y) * logistic(x %*% coef) +
                 object$p^(1 - y) * (1 - object$p)^y * (1 - logistic(x %*% coef)) )
          
    } else if (object$design == "unrelated-known") {

        pred <- ( (object$p * y + (1 - object$p) * object$q^y * (1 - object$q)^(1 - y))
                 * logistic(x %*% coef) ) /
        ( (object$p * logistic(x %*% coef)^y * (1 - logistic(x %*% coef))^(1 - y)) +
           ((1 - object$p) * object$q^y * (1 - object$q)^(1 - y)) )

    }
    
    return(pred)
    
  }
  
  if(missing(newdata)) {
    xvar <- object$x
    y <- object$y
  } else {
    if(nrow(newdata)==0)
      stop("No data in the provided data frame.")
    xvar <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata)
    
    if(given.y == TRUE){
    newdata1 <- model.frame(object$call$formula, newdata, na.action = na.omit)
    xvar <- model.matrix(as.formula(paste("~", c(object$call$formula[[3]]))), newdata1)
    y <- model.response(newdata1)
    }
  }
  
  if(quasi.bayes == FALSE) {
    
    if(given.y == TRUE) {
      pred <- z.post(object, x = xvar, y = y)
    } else {
      pred <- logistic(xvar %*% coef(object))
    }
    
    if(avg == TRUE)
      pred <- t(as.matrix(apply(pred, 2, mean)))
    
    return.object <- list(est = pred)
    
  } else {
    
    draws <- mvrnorm(n = n.sims, mu = object$est, Sigma = vcov(object))
    
    if(given.y == TRUE) {
      pred <- z.post(object, x = xvar, y = y, coef = t(draws))
    } else {
      pred <- logistic(xvar %*% t(draws))
    }
    
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


summary.rrreg <- function(object, ...) {
  structure(object, class = c("summary.rrreg", class(object)))
}

print.rrreg <- function(x, ...){

  cat("\nRandomized Response Technique Regression \n\nCall: ")

  dput(x$call)

  cat("\n")

  tb <- matrix(NA, ncol = 2, nrow = length(x$est))
  colnames(tb) <- c("est", "se")
  rownames(tb) <- colnames(x$x)

  tb[,1] <- x$est
  tb[,2] <- x$se

  print(as.matrix(round(tb,5)))

  cat("\n")

  summarize.design(x$design, x$p, x$p0, x$p1, x$q)

  cat("\n")

  invisible(x)

}

print.summary.rrreg <- function(x, ...){

  cat("\nRandomized Response Technique Regression \n\nCall: ")

  dput(x$call)

  cat("\n")

  tb <- matrix(NA, ncol = 2, nrow = length(x$est))
  colnames(tb) <- c("Est.", "S.E.")
  rownames(tb) <- colnames(x$x)

  tb[,1] <- x$est
  tb[,2] <- x$se

  print(as.matrix(round(tb,5)))

  cat("\n")

  summarize.design(x$design, x$p, x$p0, x$p1, x$q)

  cat("\n")

  invisible(x)

}

