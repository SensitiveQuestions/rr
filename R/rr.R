logistic <- function(x) exp(x)/(1+exp(x))

#' Randomized Response Regression
#' 
#' \code{rrreg} is used to conduct multivariate regression analyses of survey
#' data using randomized response methods.
#' 
#' This function allows users to perform multivariate regression analysis on
#' data from the randomized response technique.  Four standard designs are
#' accepted by this function: mirrored question, forced response, disguised
#' response, and unrelated question. The method implemented by this function is
#' the Maximum Likelihood (ML) estimation for the Expectation-Maximization (EM)
#' algorithm.
#' 
#' @usage rrreg(formula, p, p0, p1, q, design, data, start = NULL, maxIter =
#' 10000, verbose = FALSE, optim = FALSE, em.converge = 10^(-8), glmMaxIter =
#' 10000, solve.tolerance = .Machine$double.eps)
#' @param formula An object of class "formula": a symbolic description of the
#' model to be fitted.
#' @param p The probability of receiving the sensitive question (Mirrored
#' Question Design, Unrelated Question Design); the probability of answering
#' truthfully (Forced Response Design); the probability of selecting a red card
#' from the 'yes' stack (Disguised Response Design). For "mirrored" and
#' "disguised" designs, p cannot equal .5.
#' @param p0 The probability of forced 'no' (Forced Response Design).
#' @param p1 The probability of forced 'yes' (Forced Response Design).
#' @param q The probability of answering 'yes' to the unrelated question, which
#' is assumed to be independent of covariates (Unrelated Question Design).
#' @param design One of the four standard designs: "forced-known", "mirrored",
#' "disguised", or "unrelated-known".
#' @param data A data frame containing the variables in the model.
#' @param start Optional starting values of coefficient estimates for the
#' Expectation-Maximization (EM) algorithm.
#' @param maxIter Maximum number of iterations for the Expectation-Maximization
#' algorithm. The default is \code{10000}.
#' @param verbose A logical value indicating whether model diagnostics counting
#' the number of EM iterations are printed out.  The default is \code{FALSE}.
#' @param optim A logical value indicating whether to use the quasi-Newton
#' "BFGS" method to calculate the variance-covariance matrix and standard
#' errors. The default is \code{FALSE}.
#' @param em.converge A value specifying the satisfactory degree of convergence
#' under the EM algorithm. The default is \code{10^(-8)}.
#' @param glmMaxIter A value specifying the maximum number of iterations to run
#' the EM algorithm. The default is \code{10000}.
#' @param solve.tolerance When standard errors are calculated, this option
#' specifies the tolerance of the matrix inversion operation solve.
#' @return \code{rrreg} returns an object of class "rrreg".  The function
#' \code{summary} is used to obtain a table of the results.  The object
#' \code{rrreg} is a list that contains the following components (the inclusion
#' of some components such as the design parameters are dependent upon the
#' design used):
#' 
#' \item{est}{Point estimates for the effects of covariates on the randomized
#' response item.} \item{vcov}{Variance-covariance matrix for the effects of
#' covariates on the randomized response item.} \item{se}{Standard errors for
#' estimates of the effects of covariates on the randomized response item.}
#' \item{data}{The \code{data} argument.} \item{coef.names}{Variable names as
#' defined in the data frame.} \item{x}{The model matrix of covariates.}
#' \item{y}{The randomized response vector.} \item{design}{Call of standard
#' design used: "forced-known", "mirrored", "disguised", or "unrelated-known".}
#' \item{p}{The \code{p} argument.} \item{p0}{The \code{p0} argument.}
#' \item{p1}{The \code{p1} argument.} \item{q}{The \code{q} argument.}
#' \item{call}{The matched call.}
#' @seealso \code{\link{predict.rrreg}} for predicted probabilities.
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design
#' and Analysis of the Randomized Response Technique." Working Paper. Available
#' at \url{http://imai.princeton.edu/research/randresp.html}.
#' @keywords regression
#' @examples
#' 
#' \dontrun{
#' data(nigeria)
#' 
#' set.seed(1)
#' 
#' ## Define design parameters
#' p <- 2/3  # probability of answering honestly in Forced Response Design
#' p1 <- 1/6 # probability of forced 'yes'
#' p0 <- 1/6 # probability of forced 'no'
#' 
#' ## Fit linear regression on the randomized response item of whether 
#' ## citizen respondents had direct social contacts to armed groups
#' 
#' rr.q1.reg.obj <- rrreg(rr.q1 ~ cov.asset.index + cov.married + 
#'                     I(cov.age/10) + I((cov.age/10)^2) + cov.education + cov.female,   
#'                     data = nigeria, p = p, p1 = p1, p0 = p0, 
#'                     design = "forced-known")
#'   
#' summary(rr.q1.reg.obj)
#' 
#' ## Replicates Table 3 in Blair, Imai, and Zhou (2014)
#' }
#' 
#' @export
rrreg <- function(formula, p, p0, p1, q, design, data, start = NULL, 
                  maxIter = 10000, verbose = FALSE, 
                  optim = FALSE, em.converge = 10^(-8), 
                  glmMaxIter = 10000, solve.tolerance = .Machine$double.eps, 
                  h = NULL, group = NULL, use.gmm = FALSE, matrixMethod = "efficient") {
  
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

    if (use.gmm == FALSE) {
      
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
  
           information.rr <- crossprod(score.rr)/nrow(x)
           vcov <- solve(information.rr, tol = solve.tolerance)/nrow(x)
           vcov <- vcov/(bscale %*% t(bscale))
           
           se <- sqrt(diag(vcov))
      }
      
    } else {

      bpar <- start

    }
    
    ## auxiliary data functionality checks 
    if (!is.null(h) & is.null(group)) {
      stop("Need to specify character vector of group assignments. See documentation")
    } else if (is.null(h) & !is.null(group)) {
      stop("Need to specify named vector of group moments. See documentation")
    } 

    aux.check <- !(is.null(h) | is.null(group))

    ## AUXILIARY DATA FUNCTIONALITY  
    if (aux.check) {     

        score <- function (beta, y, x, w = NULL, c, d ) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            coef <- (c * y / (c * logistic(x %*% beta) + d) - 
                c * (1 - y) / (1 - c * logistic(x %*% beta) - d)) * c(dlogistic.coef(x %*% beta))
            w.coef <- c(coef) * w
            colSums(w.coef * x)/sum(w)
        }

        jbn <- function (beta, y, x, w = NULL, c, d ) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            coef1 <- (y * c/(c * logistic(x %*% beta) + d) -
                (1 - y) * c /(1 - c * logistic(x %*% beta) - d)) *
                    exp(x %*% beta) * (1 - exp(2 * x %*% beta)) / (1 + exp(x %*% beta))^4
            coef2 <- (y * c^2 / (c * logistic(x %*% beta) + d)^2 + 
                (1 - y) * c^2/(1 - c * logistic(x %*% beta) - d)^2) *
                    exp(2 * x %*% beta)/(1 + exp(x %*% beta))^4
            t(w * c(coef1 - coef2) * x) %*% x / sum(w)
        }

        # gmm objective function
        gmm <- function (beta, y, x, w = NULL, c, d, W = NULL, h, group) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            m1 <- score(beta, y, x, w, c, d)
            g <- c()
            group.labels <- names(h)
            for (i in 1:length(h)){
                wg <- w[group == group.labels[i]]
                g[i] <- sum(c(h[i]) - logistic(x[group == group.labels[i], , drop = FALSE] %*% beta))/sum(w)
            }
            if (missing(W)) W <- diag(length(c(m1, g)))
            return(t(c(m1, g)) %*% W %*% c(m1, g))
        }

        # gmm gradient
        grad <- function (beta, y, x, w = NULL, c, d , W = NULL, h, group) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            score <- score(beta, y, x, w, c, d)
            jbn <- jbn(beta, y, x, w, c, d)
            g <- c()
            group.labels <- names(h)
            for (i in 1:length(h)){
                wg <- w[group == group.labels[i]]
                g[i] <- sum(c(h[i]) - logistic(x[group == group.labels[i], , drop = FALSE] %*% beta))/sum(w)
                grad.iter <- colSums(-c(dlogistic.coef(x[group == group.labels[i], , drop = FALSE] %*% beta)) * x[group == group.labels[i], , drop = FALSE])/sum(w)
                jbn <- rbind(jbn, grad.iter)
            } 
            if (missing(W)) W <- diag(length(c(score, g)))
            t(2 * c(score, g)) %*% W %*% jbn
        }

        # gmm weighting matrix
        weight.matrix <- function (beta, y, x, w = NULL, c, d, h, group, matrixMethod) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            coef <- (c * y / (c * logistic(x %*% beta) + d) - 
                c * (1 - y) / (1 - c * logistic(x %*% beta) - d)) * c(dlogistic.coef(x %*% beta))
            w.coef <- c(coef) * w
            dg <- c()
            group.labels <- names(h)
            for (i in 1:length(group.labels)){
                wg <- w[group == group.labels[i]]
                dg[i] <- sum(c(h[i]) - logistic(x[group == group.labels[i], , drop = FALSE] %*% beta)^2)/sum(w)
            }
            matrix1 <- t(w.coef * x) %*% (w.coef * x)/sum(w)
            matrix2 <- diag(dg, nrow = length(dg))
            efficient.W <- solve(adiag(matrix1, matrix2))
            if (matrixMethod == "efficient" | matrixMethod == "cue") {
              efficient.W
            } else if (matrixMethod == "princomp") {
              decomp <- eigen(efficient.W)
              decomp$values * t(decomp$vectors) %*% decomp$vectors
            }
        }

        # continuously updating objective function
        gmm.cue <- function (beta, y, x, w = NULL, c, d, h, group, matrixMethod) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            m1 <- score(beta, y, x, w, c, d)
            g <- c()
            group.labels <- names(h)
            for (i in 1:length(h)){
                wg <- w[group == group.labels[i]]
                g[i] <- sum(c(h[i]) - logistic(x[group == group.labels[i], , drop = FALSE] %*% beta))/sum(w)
            }
            W <- weight.matrix(beta = beta, y = y, x = x, w = w, c = c, d = d, 
              h = h, group = group, matrixMethod = matrixMethod)
            t(c(m1, g)) %*% W %*% c(m1, g)
        }


        # coefs with auxiliary information
        rr.true.coefs <- function (beta, y, x, w = NULL, c, d , W = NULL, h = NULL, group = NULL, matrixMethod) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            if (missing(W) & missing(h)) {
                W <- diag(length(ncol(x)))
            } else if (missing(W)) {
                W <- diag(length(c(ncol(x), h)))
            }
            if (matrixMethod == "efficient" | matrixMethod == "princomp") {
              weight.matrix <- weight.matrix(beta = beta, y = y, x = x, h = h, group = group, c = c, d = d, matrixMethod = matrixMethod)
              fit1 <- optim(par = beta, fn = gmm, gr = grad, method = "BFGS", 
                control = list(maxit = 500), y = y, x = x, h = h, group = group, 
                c = c, d = d, W = weight.matrix)
              newpar <- fit1$par
              weight.matrix <- weight.matrix(beta = newpar, y = y, x = x, h = h, group = group, c = c, d = d, matrixMethod = matrixMethod)
              fit2 <- optim(par = newpar, fn = gmm, gr = grad, method = "BFGS", 
                control = list(maxit = 500), y = y, x = x, h = h, group = group, 
                c = c, d = d, W = weight.matrix)
              fit2             
            } else {
              fit <- optim(par = beta, fn = gmm.cue, method = "BFGS", 
                control = list(maxit = 500), y = y, x = x, h = h, group = group, 
                c = c, d = d, matrixMethod = matrixMethod)
              fit
            }
        }

        # vcov with auxiliary information
        rr.true.vcov <- function(beta, y, x, w, c, d, h = NULL, group = NULL, matrixMethod) {
            n <- length(y)
            if (missing(w)) w <- rep(1, n)
            jbn <- jbn(beta, y, x, w, c, d)
            group.labels <- names(h)
            for (i in 1:length(group.labels)) {
                wg <- w[group == group.labels[i]]
                grad.iter <- t(w[group == group.labels[i]]) %*% (-c(dlogistic.coef(x[group == group.labels[i], , drop = FALSE] %*% beta)) * x[group == group.labels[i], , drop = FALSE])/sum(w)
                jbn <- rbind(jbn, grad.iter)
            }
            weight.matrix <- weight.matrix(beta = beta, y = y, x = x, h = h, group = group, c = c, d = d, matrixMethod = matrixMethod)
            if (matrixMethod != "princomp") {
              solve(t(jbn) %*% weight.matrix %*% jbn)/(sum(w))
            } else {
              V <- solve(weight.matrix(beta = beta, y = y, x = x, h = h, group = group, c = c, d = d, matrixMethod = "efficient"))
              solve(t(jbn) %*% weight.matrix %*% jbn) %*% t(jbn) %*% weight.matrix %*% V %*% t(weight.matrix) %*% jbn %*% solve(t(jbn) %*% weight.matrix %*% jbn)/sum(w)
            }
        }



        true.fit <- rr.true.coefs(beta = bpar, y = y, x = x1, c = c, d = d, h = h, group = group, matrixMethod = matrixMethod)
        bpar <- true.fit$par
        vcov <- rr.true.vcov(beta = bpar, y = y, x = x1, c = c, d = d, h = h, group = group, matrixMethod = matrixMethod)
        se <- sqrt(diag(vcov))

        # Sargan-Hansen overidentification test
        J.stat <- true.fit$val * n
        overid.p <- round(1 - pchisq(J.stat, df = length(h)), 4)

    }

    return.object <- list(est = bpar, 
                          vcov = vcov,
                          se = se,
                          data = df,
                          coef.names = colnames(x),
                          x = x1, #unscaled data
                          y = y,
                          design = design,
                          call = match.call(), 
                          aux = aux.check)

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

    if (aux.check) {
      return.object$nh <- length(h)
      return.object$wm <- ifelse(matrixMethod == "cue", "continuously updating", 
        ifelse(matrixMethod == "princomp", "principal components", "efficient"))
      return.object$J.stat <- round(J.stat, 4)
      return.object$overid.p <- overid.p
    }
    
    class(return.object) <- "rrreg"

    return(return.object)

}

#' @export
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

#' @export
vcov.rrreg <- function(object, ...){

  vcov <- object$vcov

  rownames(vcov) <- colnames(vcov) <- object$coef.names

  return(vcov)

}

#' @export
coef.rrreg <- function(object, ...){

  coef <- object$est

  names(coef) <- object$coef.names

  return(coef)

}



#' Predicted Probabilities for Randomized Response Regression
#' 
#' \code{predict.rrreg} is used to generate predicted probabilities from a
#' multivariate regression object of survey data using randomized response
#' methods.
#' 
#' This function allows users to generate predicted probabilities for the
#' randomized response item given an object of class "rrreg" from the
#' \code{rrreg()} function. Four standard designs are accepted by this
#' function: mirrored question, forced response, disguised response, and
#' unrelated question. The design, already specified in the "rrreg" object, is
#' then directly inputted into this function.
#' 
#' @usage predict.rrreg(object, given.y = FALSE, alpha = .05, n.sims =
#' 1000, avg = FALSE, newdata = NULL, quasi.bayes = FALSE, keep.draws = FALSE,
#' ...)
#' @param object An object of class "rrreg" generated by the \code{rrreg()}
#' function.
#' @param given.y Indicator of whether to use "y" the response vector to
#' calculate the posterior prediction of latent responses. Default is
#' \code{FALSE}, which simply generates fitted values using the logistic
#' regression.
#' @param alpha Confidence level for the hypothesis test to generate upper and
#' lower confidence intervals. Default is \code{ .05}.
#' @param n.sims Number of sampled draws for quasi-bayesian predicted
#' probability estimation. Default is \code{1000}.
#' @param avg Whether to output the mean of the predicted probabilities and
#' uncertainty estimates. Default is \code{FALSE}.
#' @param newdata Optional new data frame of covariates provided by the user.
#' Otherwise, the original data frame from the "rreg" object is used.
#' @param quasi.bayes Option to use Monte Carlo simulations to generate
#' uncertainty estimates for predicted probabilities. Default is \code{FALSE}.
#' @param keep.draws Option to return the Monte Carlos draws of the quantity of
#' interest, for use in calculating differences for example.
#' @param ... Further arguments to be passed to \code{predict.rrreg()} command.
#' @return \code{predict.rrreg} returns predicted probabilities either for each
#' observation in the data frame or the average over all observations. The
#' output is a list that contains the following components:
#' 
#' \item{est}{Predicted probabilities for the randomized response item
#' generated either using fitted values, posterior predictions, or
#' quasi-Bayesian simulations. If \code{avg} is set to \code{TRUE}, the output
#' will only include the mean estimate.} \item{se}{Standard errors for the
#' predicted probabilities of the randomized response item generated using
#' Monte Carlo simulations. If \code{quasi.bayes} is set to \code{FALSE}, no
#' standard errors will be outputted.} \item{ci.lower}{Estimates for the lower
#' confidence interval. If \code{quasi.bayes} is set to \code{FALSE}, no
#' confidence interval estimate will be outputted.} \item{ci.upper}{Estimates
#' for the upper confidence interval. If \code{quasi.bayes} is set to
#' \code{FALSE}, no confidence interval estimate will be outputted.}
#' \item{qoi.draws}{Monte Carlos draws of the quantity of interest, returned
#' only if \code{keep.draws} is set to \code{TRUE}.}
#' @seealso \code{\link{rrreg}} to conduct multivariate regression analyses in
#' order to generate predicted probabilities for the randomized response item.
#' @references Blair, Graeme, Kosuke Imai and Yang-Yang Zhou. (2014) "Design
#' and Analysis of the Randomized Response Technique."  \emph{Working Paper.}
#' Available at \url{http://imai.princeton.edu/research/randresp.html}.
#' @keywords predicted probabilities fitted values
#' @examples
#' 
#' \dontrun{
#' data(nigeria)
#' 
#' set.seed(1)
#' 
#' ## Define design parameters
#' p <- 2/3  # probability of answering honestly in Forced Response Design
#' p1 <- 1/6 # probability of forced 'yes'
#' p0 <- 1/6 # probability of forced 'no'
#' 
#' ## Fit linear regression on the randomized response item of 
#' ## whether citizen respondents had direct social contacts to armed groups
#' 
#' rr.q1.reg.obj <- rrreg(rr.q1 ~ cov.asset.index + cov.married + I(cov.age/10) + 
#'                       I((cov.age/10)^2) + cov.education + cov.female,   
#'                       data = nigeria, p = p, p1 = p1, p0 = p0, 
#'                       design = "forced-known")
#' 
#' ## Generate the mean predicted probability of having social contacts to 
#' ## armed groups across respondents using quasi-Bayesian simulations. 
#' 
#' rr.q1.reg.pred <- predict(rr.q1.reg.obj, given.y = FALSE, 
#'                                 avg = TRUE, quasi.bayes = TRUE, 
#'                                 n.sims = 10000)
#' 
#' ## Replicates Table 3 in Blair, Imai, and Zhou (2014)
#' }
#' 
#' @importFrom MASS mvrnorm
#' 
#' @method predict rrreg 
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

  if (x$aux) cat("Incorporating ", x$nh, " auxiliary moment(s). Weighting method: ", x$wm, ".\n", 
    "The overidentification test statistic was: ", x$J.stat, " (p < ", x$overid.p, ")", ".\n", sep = "")

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

  if (x$aux) cat("Incorporating ", x$nh, " auxiliary moment(s). Weighting method: ", x$wm, ".\n", 
    "The overidentification test statistic was: ", x$J.stat, " (p < ", x$overid.p, ")", ".\n", sep = "")

  invisible(x)

}

