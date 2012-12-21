##====Beta-binomial Count Simulations====

## I use cat settings.R betabinsim.R > bbsim.R in the .csh file to get full sim file
getwd()
ls() 
##----Global parameter settings----
sims     <- 10000
depthMin <- 5000
depthMax <- 10000

##----Necessities----
require(MASS)
require(betareg)
require(aod)
source("~/Metagenomics/metastats.R")

##----Beta-binomial Count Generator----
# Generates binomial counts where the probability of success is the reparameterized beta density
rbetabin <- function(n, depth, mu, theta){ # theta is precision ( or 1/dispersion)
    y <- rep(NA, n)
    
    shape1 <- mu * theta
    shape2 <- (1 - mu) * theta
    
    y <- rbinom(n, depth, rbeta(n, shape1, shape2))
    return(y)
}

##----Random uniform discrete generator for generating depths----
runifdisc <- function(n, min=0, max=1) sample(min:max, n, replace=T)



##----Results List----
results.bb <- list(
    details   = data.frame(sims, n, depthMin, depthMax, mu, a),
    precision = matrix(rep(NA, 2 * sims), nrow = sims, dimnames = list(1:sims, c("group1", "group2"))),
    counts    = matrix(rep(NA, 2 * n * sims), nrow = sims),
    depths    = matrix(rep(NA, 2 * n * sims), nrow = sims),
    ttest     = matrix(rep(NA, 2 * sims), nrow = sims, dimnames = list(1:sims, c("p.tt", "s.tt"))),
    metastat  = rep(NA, sims),
    poisson   = data.frame(matrix(rep(NA, 3 * sims), nrow = sims, 
                                  dimnames = list(1:sims, c("p.po", "s.po", "c.po")))),
    negbin    = data.frame(matrix(rep(NA, 3 * sims), nrow = sims, 
                                  dimnames = list(1:sims, c("p.nb", "s.nb", "c.nb")))),
    betaBin   = data.frame(matrix(rep(NA, 7 * sims), nrow = sims, 
                                  dimnames = list(1:sims, c("l.bb", "p.bbv", "s.bbv", "c.bbv",
                                                            "p.bbf", "s.bbf", "c.bbf")))),
    betaReg   = data.frame(matrix(rep(NA, 7 * sims), nrow = sims, 
                                  dimnames = list(1:sims, c("l.br", "p.brv", "s.brv", "c.brv",
                                                            "p.brf", "s.brf", "c.brf")))),
    perf      = list(
        p.na      = rep(NA, 10),
        false.pos = rep(NA, 8),
        true.neg  = rep(NA, 8),
        false.neg = rep(NA, 8),
        true.pos  = rep(NA, 8)
    )     
    
)

set.seed(2001)
for(i in 1:sims){
    ## Sampling depth of each subject
    depth1 <- runifdisc(n, depthMin, depthMax) # group 1 depths
    depth2 <- runifdisc(n, depthMin, depthMax) # group 2 depths
    
    ## Variable precision between groups
    theta1 <- 1/runif(1, 1e-08, 0.05) -> results.bb$precision[i, "group1"]
    theta2 <- 1/runif(1, 1e-08, 0.05) -> results.bb$precision[i, "group2"]
    
    ## Generate beta-binomial counts
    if(i <= round(sims/2)){
        y1 <- rbetabin(n, depth1, mu, theta1)     
        y2 <- rbetabin(n, depth2, mu, theta2) # both groups have equal means
    } 
    else{
        y1 <- rbetabin(n, depth1, mu, theta1)     
        y2 <- rbetabin(n, depth2, a * mu, theta2) # group 2 mu is a * group 1 mu
    }
    
    ## t-test
    tt   <- t.test(y1/depth1, y2/depth2)
    results.bb$ttest[i, "p.tt"] <- tt$p.value    
    results.bb$ttest[i, "s.tt"] <- tt$statistic  
    
    ## Metastats
    # Actually contructs a 2 row feature matrix and does 2 tests
    # I only care about the first test. Need to do both to easily implement Metastats code
    aa     <- matrix(c(y1, y2, depth1 - y1, depth2 - y2), nrow = 2, byrow = TRUE)
    jobj   <- list(matrix = aa, taxa = c("1", "2"))
    result <- detect_differentially_abundant_features(jobj, n + 1, "t.diffAb", B = 1000)
    results.bb$metastat[i] <- result[1] 
    
    ## Helpful for upcoming regression models
    group <- factor(rep(1:2, each = n))
    y     <- c(y1, y2)         -> results.bb$counts[i, ]
    depth <- c(depth1, depth2) -> results.bb$depths[i, ]
    
    ## Poisson regression
    # Returns a p-value and ChiSq-stat for the beta1 parameter in the poisson model
    # Although it's not using the SFR package, this is exactly what SFR does for this type of comparison
    # Also returns a logical for whether or not the model converged
    mp   <- glm(y ~ offset(log(depth)) + group, family = "poisson")
    bb   <- anova.glm(mp, test = "Chisq")
    results.bb$poisson[i, "p.po"] <- bb["group", "Pr(>Chi)"]
    results.bb$poisson[i, "s.po"] <- bb["group", "Deviance"]
    results.bb$poisson[i, "c.po"] <- mp$converged
    
    ## Negative binomial regression
    # Returns a p-value and ChiSq-stat for the beta1 parameter in the negative binomial model
    # Also returns a logical for whether or not the model converged
    mn   <- glm.nb(y ~ offset(log(depth)) + group)
    cc   <- anova.glm(mn, test = "Chisq")
    results.bb$negbin[i, "p.nb"] <- cc["group", "Pr(>Chi)"]
    results.bb$negbin[i, "s.nb"] <- cc["group", "Deviance"]
    results.bb$negbin[i, "c.nb"] <- mn$converged
    
    ## Beta binomial regression
    # Assumes a fixed dispersion across groups
    mbbv <- betabin(cbind(y, depth - y) ~ group, ~ group, data = data.frame(y, depth, group))
    mbbf <- betabin(cbind(y, depth - y) ~ group, ~ 1, data = data.frame(y, depth, group))
    # p-value for likelihood ratio test (fixed or variable model?)
    results.bb$betaBin[i, "l.bb"] <- pchisq(-2 * (logLik(mbbf) - logLik(mbbv))[1], df = 1, lower.tail = FALSE)
    # variable and fixed model results
    results.bb$betaBin[i, "p.bbv"] <- summary(mbbv)@Coef["group2", "Pr(> |z|)"]
    results.bb$betaBin[i, "p.bbf"] <- summary(mbbf)@Coef["group2", "Pr(> |z|)"]
    results.bb$betaBin[i, "s.bbv"] <- summary(mbbv)@Coef["group2", "z value"]
    results.bb$betaBin[i, "s.bbf"] <- summary(mbbf)@Coef["group2", "z value"]
    results.bb$betaBin[i, "c.bbv"] <- mbbv@code
    results.bb$betaBin[i, "c.bbf"] <- mbbf@code
    
    ## Beta regression (variable and fixed dispersions)
    prop <- y/depth
    prop[prop==0] <- .1/depth[prop==0] #zero-count correction (.1 is ideal amount)
    mbrv   <- betareg(prop ~ group | group, control = betareg.control(maxit = 10000))
    mbrf   <- betareg(prop ~ group, link.phi = "log", control = betareg.control(maxit = 25000))
    # p-value for likelihood ratio test (fixed or variable model?)
    results.bb$betaReg[i, "l.br"] <- pchisq(-2 * (logLik(mbrf) - logLik(mbrv))[1], df = 1, lower.tail = FALSE)
    # variable and fixed model results
    results.bb$betaReg[i, "p.brv"] <- summary(mbrv)$coefficients$mean["group2", "Pr(>|z|)"]
    results.bb$betaReg[i, "p.brf"] <- summary(mbrf)$coefficients$mean["group2", "Pr(>|z|)"]
    results.bb$betaReg[i, "s.brv"] <- summary(mbrv)$coefficients$mean["group2", "z value"]
    results.bb$betaReg[i, "s.brf"] <- summary(mbrf)$coefficients$mean["group2", "z value"]
    results.bb$betaReg[i, "c.brv"] <- mbrv$converged
    results.bb$betaReg[i, "c.brf"] <- mbrf$converged    
}

## write out workspace and p-value results
pvalues <- data.frame(p.tt  = results.bb$ttest[, "p.tt"],
                      p.ms  = results.bb$metastat,
                      p.po  = results.bb$poisson[, "p.po"],
                      p.nb  = results.bb$negbin[, "p.nb"],
                      l.bb  = results.bb$betaBin[, "l.bb"],
                      p.bbv = results.bb$betaBin[, "p.bbv"],
                      p.bbf = results.bb$betaBin[, "p.bbf"],
                      l.br  = results.bb$betaReg[, "l.br"],
                      p.brv = results.bb$betaReg[, "p.brv"],
                      p.brf = results.bb$betaReg[, "p.brf"]
                      )
write.table(pvalues, paste(sim.id, "pvalues.txt", sep = ""), sep = "\t", row.names = FALSE)

## How many methods have NAs in their p-values
results.bb$perf$p.na <- apply(pvalues, 2, FUN = function(x) sum(is.na(x)))


## At an alpha=0.05 level how did the methods do? (Note that NAs are removed)
cutoff <- round(sims/2) # cutoff between sames and differences
# False Postives; How many known 'sames' (group1 = group2) did methods get wrong
results.bb$perf$false.pos <- apply(pvalues[1:cutoff, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] < 0.05))
# True Negatives; How many known 'sames' (group1 = group2) did methods get right
results.bb$perf$true.neg <- apply(pvalues[1:cutoff, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] >= 0.05))
# False Negatives; How many known 'differences' (group1 != group2) did methods get wrong
results.bb$perf$false.neg <- apply(pvalues[(cutoff + 1):sims, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] >= 0.05))
# True Positives; How many known 'differences' (group1 != group2) dit methods get right
results.bb$perf$true.pos <- apply(pvalues[(cutoff + 1):sims, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] < 0.05))

## give results specific sim.id name
assign(paste("results", sim.id, sep = "."), results.bb)

## save workspace image containing results
save(list = paste("results", sim.id, sep = "."), file = paste(sim.id, "RData", sep = "."))

## cleanup
rm(list =ls())
