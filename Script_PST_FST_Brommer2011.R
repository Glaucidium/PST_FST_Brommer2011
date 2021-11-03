##############################
# PST-FST robustness test by Brommer (2011)
##############################

library(VCA)
library(boot)
library(ggplot2)
library(stats)

############
# Pstformula
############

# data = dataframe with dependent (trait) and independent variables.
# d = index for bootstrap function (fixed)
# ch = c/h^2 variable

Pstformula <- function (data, d, ch){  
  DBind <- data[d, ]
  anova.res <- anovaVCA(trait ~ sex + pop + who, Data = DBind)$aov.tab # change name of traits of interest acording to the infile data (covariables)
  pstbych <- ch * anova.res['pop','VC'] / 
    (ch * anova.res['pop','VC'] + (2 * anova.res['error','VC'])) # PST formula for a fixed ch
  return(pstbych)
}

#############
# PstBrommer
#############

# data = dataframe with dependent (trait) and independent variables
# a = to use in the seq function, the lower c/h2 value of interest (from)
# b = to use in the seq function, the higher c/h2 value of interest (to)
# c = to use in the seq function, the increment of the sequence (by)
# boot.rep = number of bootstrap replicates.


PstBrommer <- function(data, a = 0.01, b = 2, c = 0.01, boot.rep) {
  chvector <- seq(a, b, by =  c) # to determine vector of possible ch values
  pst <- as.numeric(rep(NA, length(chvector))) # empty vector to store PST results
  CIup <- as.numeric(rep(NA, length(chvector))) # empty vector to store ci 95% PST
  CIlow <- as.numeric(rep(NA, length(chvector))) # empty vector to store ci 5% PST
  for (i in 1:length(chvector)){
    ch <- chvector[i]
    boot.out <- boot(data = data, statistic = Pstformula, R = boot.rep, ch = ch, parallel = "snow")
    pst[i] <- boot.out$t0 # boot resultado pst
    CIup[i] <- boot.ci(boot.out, type = 'perc')$percent[[5]]# boot.ci 95% upper interval
    CIlow[i] <- boot.ci(boot.out, type = 'perc')$percent[[4]]# boot.ci 5% lower interval
    print(i)
  }
  pstdf <- data.frame(chvector, pst, CIlow, CIup)
  return(pstdf)
}

#################
## Example of use
#################

# example dataframe with one trait of interest and covariables. Remember to check ANOVA assuptions!
df.pop1 <- data.frame (pop = as.factor(morpho.data$Clas), 
                       sex = as.factor(morpho.data$Sex), 
                       trait = as.numeric(morpho.data$PC1), 
                       who = as.factor(morpho.data$Who))

# PST value for a fixed c/h2. 
boot.out <- boot(data = df.pop1, statistic = Pstformula, R = 100000, ch = 1, parallel = "snow") # PST for c/h2 = 1
boot.out$t0
boot.ci(boot.out, type = 'perc')$percent[[5]]
boot.ci(boot.out, type = 'perc')$percent[[4]]

# PST for each c/h2 of interest. In this example, c/h2 will take values from 0.01 to 2 each 0.02
pstdf <- PstBrommer (data = df.pop1, boot.rep = 50, a = 0.01, b = 2, c = 0.02 )

#########
# plot
#########

plot(pstdf$ch, pstdf$pst, type = "n", 
     ylim = c(0, 1), xlim = c(0, 2), 
     xlab = as.expression(expression("c/"*"h"^"2")), 
     ylab = as.expression(expression("P"[ST] *" or F"[ST])),
     xaxs="r", yaxs="r",
     xaxt="n",
     main = 'DDRAD data pop1 PC1',
     bty = "n")
axis (side = 1, at = seq(0, 2, 0.2))

# PST lines
lines(smooth.spline(pstdf$ch, pstdf$CIlow, spar = 0.60),
      col = "blue", lty = "dotted", lwd = 1.5)
lines(pstdf$ch, pstdf$pst, col = "blue", lwd = 2)
lines(smooth.spline(pstdf$ch, pstdf$CIup, spar = 0.60),
      col = "blue", lty = "dotted", lwd = 1.5)

# FST lines -  mean and CI
abline (h = 0.086, col = "red", lwd =2) # mean 
abline (h = 0.088, col = "red", lty = "dashed", lwd = 1.5) # 95% CI
abline (h = 0.083, col = "red", lty = "dashed", lwd = 1.5) # 5% CI

ch2critic <- approxfun(pstdf$CIlow, pstdf$ch) # creates a formula to determine c/h2 critic
ch2critic(0.088) # Place between parenthesis the upper CI value of FST

