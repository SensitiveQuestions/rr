logistic <- function(x) exp(x)/(1+exp(x))

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


vcov.rrreg.predictor <- function(object, ...){
  
  vcov <- object$vcov
  
  rownames(vcov) <- colnames(vcov) <- c(colnames(object$x), object$coef.names)
  
  return(vcov)
  
}

coef.rrreg.predictor <- function(object, ...){
  
  coef <- object$est.t
  
  names(coef) <- object$coef.names
  
  return(coef)
  
}

summary.rrreg.predictor <- function(object, ...) {
  structure(object, class = c("summary.rrreg.predictor", class(object)))
}

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
  

