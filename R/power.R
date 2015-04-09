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