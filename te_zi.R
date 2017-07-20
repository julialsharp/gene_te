
##
## Zero Inflated Models
##
source("te_common.R")
library(pscl)
library(boot)

zdata <- te_data2 %>%  
  dplyr::select (asd, devreg, cns, housekeeping, wgc, te)

# Negative Binomial
zinb_fit <- zeroinfl(te ~ asd + devreg + cns + housekeeping + wgc, data = zdata, dist = "negbin", link = "logit")
summary(zinb_fit)

# Poisson
zip_fit <- zeroinfl(te ~ .| asd + devreg + cns + housekeeping + wgc, data = zdata)
summary(zip_fit)

# null model
zip_null <- update(zip_fit, . ~ 1)

# df = 4 for predictors in the full model
pchisq(2 * (logLik(zip_fit) - logLik(zip_null)), df = 4, lower.tail = FALSE)

# Standard Poisson model using glm
glp_fit <- glm(te ~ asd + devreg + cns + housekeeping + wgc, family = poisson, data = zdata)
summary(glp_fit)

vuong(glp_fit, zip_fit)

# print in a copyable format
dput(coef(zip_fit, "count"))
dput(coef(zip_fit, "zero"))

# bootstrap to get confidence intervals
f <- function(data, i) {
  library(pscl)
  m <- zeroinfl(te ~ .| asd + devreg + cns + housekeeping + wgc, data = data[i, ],
                start = list(count = c(4.62158803741013, 0.720112468782827, 0.117378116123572,
                                       0.433387832391776, -0.41464572078391, -0.0330734424583709), 
                             zero = c(-2.15043568742934, -0.685442348507311, -0.0108790570643863, 
                                      0.536171475367259, -0.776906673590443, 0.1979986930577)))
  as.vector(t(do.call(rbind, coef(summary(m)))[, 1:2]))
}

set.seed(100)

# number of repetitions should be larger than the number of data rows
# https://stat.ethz.ch/pipermail/r-help/2011-February/269006.html
res <- boot(zdata, f, R = 25000, parallel = "snow", ncpus = 20)

## print results
res

## basic parameter estimates with percentile and bias adjusted CIs
parms <- t(sapply(c(1, 3, 5, 7, 9, 11), 
                  function(i) {
                    out <- boot.ci(res, index = c(i, i + 1), type = c("perc", "bca"))
                    with(out, c(Est = t0, pLL = percent[4], pUL = percent[5],
                                bcaLL = bca[4], bcaLL = bca[5]))
                  }))

## add row names
row.names(parms) <- names(coef(zip_fit))
## print results
parms

## compare with normal based approximation
confint(zip_fit)

## exponentiated parameter estimates with percentile and bias adjusted CIs
exp_parms <- t(sapply(c(1, 3, 5, 7, 9), 
                     function(i) {
                       out <- boot.ci(res, index = c(i, i + 1), type = c("perc", "bca"), h = exp)
                       with(out, c(Est = t0, pLL = percent[4], pUL = percent[5],
                                   bcaLL = bca[4], bcaLL = bca[5]))
                     }))

## add row names
row.names(exp_parms) <- names(coef(zip_fit))
## print results
exp_parms


