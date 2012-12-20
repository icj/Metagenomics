##====Negative-binomial Count Simulations====

## I use cat settings.R negbinsim.R > nbsim.R in the .csh file to get full sim file
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
#source("~/Dropbox/Metagenomics/metastats.R")
source("~/Metagenomics/metastats.R")

##----Random uniform discrete generator for generating depths----
runifdisc <- function(n, min=0, max=1) sample(min:max, n, replace=T)

##----Results List----
results.nb <- list(
    details   = data.frame(sims, n, depthMin, depthMax, mu, a),
    dispersion = matrix(rep(NA, 2 * sims), nrow = sims, dimnames = list(1:sims, c("group1", "group2"))),
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

set.seed(1812)
for(i in 1:sims){
    ## Sampling depth of each subject
    depth1 <- runifdisc(n, depthMin, depthMax) # group 1 depths
    depth2 <- runifdisc(n, depthMin, depthMax) # group 2 depths
    
    ## Variable dispersion between groups
    ## R models dispersion as 1/phi
    ## I fit neg-bin models to real data and got 25th and 75th percentiles dispersions to be .15 and .62
    theta1 <- 1/runif(1, .15, .62) -> results.nb$dispersion[i, "group1"]
    theta2 <- 1/runif(1, .15, .62) -> results.nb$dispersion[i, "group2"]
    
    ## Generate negative-binomial counts
    ## while loops take care of the occasional case where y > depth
    ## the case is not occasional for the mu = 0.2, a = 4 case, which give
    ## 80% mean proportions for group2 in the known difference simulations
    if(i <= round(sims/2)){
        ii <- 1
        while(ii == 1){
            y1 <- rnbinom(n, mu = mu * depth1, size = theta1)
            if(!any(y1 >= depth1)) ii <- 2
        }
        jj <- 1
        while(jj == 1){
            y2 <- rnbinom(n, mu = mu * depth2, size = theta2) # both groups have equal means
            if(!any(y2 >= depth2)) jj <- 2
        }
    } else{
        ii <- 1
        while(ii == 1){
            y1 <- rnbinom(n, mu = mu * depth1, size = theta1)
            if(!any(y1 >= depth1)) ii <- 2
        }
        jj <- 1
        while(jj == 1){
            y2 <- rnbinom(n, mu = a * mu * depth2, size = theta2) # group 2 mu is a * group 1 mu
            if(!any(y2 >= depth2)) jj <- 2
        }  
    }
    
    ## t-test
    tt   <- t.test(y1/depth1, y2/depth2)
    results.nb$ttest[i, "p.tt"] <- tt$p.value    
    results.nb$ttest[i, "s.tt"] <- tt$statistic  
    
    ## Metastats
    # Actually contructs a 2 row feature matrix and does 2 tests
    # I only care about the first test. Need to do both to easily implement Metastats code
    aa     <- matrix(c(y1, y2, depth1 - y1, depth2 - y2), nrow = 2, byrow = TRUE)
    jobj   <- list(matrix = aa, taxa = c("1", "2"))
    result <- detect_differentially_abundant_features(jobj, n + 1, "t.diffAb", B = 1000)
    results.nb$metastat[i] <- result[1] 
    
    ## Helpful for upcoming regression models
    group <- factor(rep(1:2, each = n))
    y     <- c(y1, y2)         -> results.nb$counts[i, ]
    depth <- c(depth1, depth2) -> results.nb$depths[i, ]
    
    ## Poisson regression
    # Returns a p-value and ChiSq-stat for the beta1 parameter in the poisson model
    # Although it's not using the SFR package, this is exactly what SFR does for this type of comparison
    # Also returns a logical for whether or not the model converged
    mp   <- glm(y ~ offset(log(depth)) + group, family = "poisson")
    bb   <- anova.glm(mp, test = "Chisq")
    results.nb$poisson[i, "p.po"] <- bb["group", "Pr(>Chi)"]
    results.nb$poisson[i, "s.po"] <- bb["group", "Deviance"]
    results.nb$poisson[i, "c.po"] <- mp$converged
    
    ## Negative binomial regression
    # Returns a p-value and ChiSq-stat for the beta1 parameter in the negative binomial model
    # Also returns a logical for whether or not the model converged
    mn   <- glm.nb(y ~ offset(log(depth)) + group)
    cc   <- anova.glm(mn, test = "Chisq")
    results.nb$negbin[i, "p.nb"] <- cc["group", "Pr(>Chi)"]
    results.nb$negbin[i, "s.nb"] <- cc["group", "Deviance"]
    results.nb$negbin[i, "c.nb"] <- mn$converged
    
    ## Beta binomial regression
    # Assumes a fixed dispersion across groups
    mbbv <- betabin(cbind(y, depth - y) ~ group, ~ group, data = data.frame(y, depth, group))
    mbbf <- betabin(cbind(y, depth - y) ~ group, ~ 1, data = data.frame(y, depth, group))
    # p-value for likelihood ratio test (fixed or variable model?)
    results.nb$betaBin[i, "l.bb"] <- pchisq(-2 * (logLik(mbbf) - logLik(mbbv))[1], df = 1, lower.tail = FALSE)
    # variable and fixed model results
    results.nb$betaBin[i, "p.bbv"] <- summary(mbbv)@Coef["group2", "Pr(> |z|)"]
    results.nb$betaBin[i, "p.bbf"] <- summary(mbbf)@Coef["group2", "Pr(> |z|)"]
    results.nb$betaBin[i, "s.bbv"] <- summary(mbbv)@Coef["group2", "z value"]
    results.nb$betaBin[i, "s.bbf"] <- summary(mbbf)@Coef["group2", "z value"]
    results.nb$betaBin[i, "c.bbv"] <- mbbv@code
    results.nb$betaBin[i, "c.bbf"] <- mbbf@code
    
    ## Beta regression (variable and fixed dispersions)
    prop <- y/depth
    prop[prop==0] <- .1/depth[prop==0] #zero-count correction (.1 is ideal amount)
    mbrv   <- betareg(prop ~ group | group, control = betareg.control(maxit = 10000))
    mbrf   <- betareg(prop ~ group, link.phi = "log", control = betareg.control(maxit = 25000))
    # p-value for likelihood ratio test (fixed or variable model?)
    results.nb$betaReg[i, "l.br"] <- pchisq(-2 * (logLik(mbrf) - logLik(mbrv))[1], df = 1, lower.tail = FALSE)
    # variable and fixed model results
    results.nb$betaReg[i, "p.brv"] <- summary(mbrv)$coefficients$mean["group2", "Pr(>|z|)"]
    results.nb$betaReg[i, "p.brf"] <- summary(mbrf)$coefficients$mean["group2", "Pr(>|z|)"]
    results.nb$betaReg[i, "s.brv"] <- summary(mbrv)$coefficients$mean["group2", "z value"]
    results.nb$betaReg[i, "s.brf"] <- summary(mbrf)$coefficients$mean["group2", "z value"]
    results.nb$betaReg[i, "c.brv"] <- mbrv$converged
    results.nb$betaReg[i, "c.brf"] <- mbrf$converged    
}

## write out workspace and p-value results
pvalues <- data.frame(p.tt  = results.nb$ttest[, "p.tt"],
                      p.ms  = results.nb$metastat,
                      p.po  = results.nb$poisson[, "p.po"],
                      p.nb  = results.nb$negbin[, "p.nb"],
                      l.bb  = results.nb$betaBin[, "l.bb"],
                      p.bbv = results.nb$betaBin[, "p.bbv"],
                      p.bbf = results.nb$betaBin[, "p.bbf"],
                      l.br  = results.nb$betaReg[, "l.br"],
                      p.brv = results.nb$betaReg[, "p.brv"],
                      p.brf = results.nb$betaReg[, "p.brf"]
)


## How many methods have NAs in their p-values
results.nb$perf$p.na <- apply(pvalues, 2, FUN = function(x) sum(is.na(x)))


## At an alpha=0.05 level how did the methods do? (Note that NAs are removed)
cutoff <- round(sims/2) # cutoff between sames and differences
# False Postives; How many known 'sames' (group1 = group2) did methods get wrong
results.nb$perf$false.pos <- apply(pvalues[1:cutoff, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] < 0.05))
# True Negatives; How many known 'sames' (group1 = group2) did methods get right
results.nb$perf$true.neg <- apply(pvalues[1:cutoff, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] >= 0.05))
# False Negatives; How many known 'differences' (group1 != group2) did methods get wrong
results.nb$perf$false.neg <- apply(pvalues[(cutoff + 1):sims, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] >= 0.05))
# True Positives; How many known 'differences' (group1 != group2) did methods get right
results.nb$perf$true.pos <- apply(pvalues[(cutoff + 1):sims, -c(5, 8)], 2, FUN = function(x) sum(x[!is.na(x)] < 0.05))

## give results specific sim.id name
assign(paste("results", sim.id, sep = "."), results.nb)

## save workspace image containing results
save(list = paste("results", sim.id, sep = "."), file = paste(sim.id, "RData", sep = "."))

## write out pvalue tables
write.table(pvalues, paste(sim.id, "pvalues.txt", sep = ""), sep = "\t", row.names = FALSE)

## cleanup
rm(list =ls())
