source("te_common.R")

##
## Basic models
##
base_model <- function(r, d) {
  f_str <- paste(r, '~ asd + devreg + cns + housekeeping + wgc')
  fit <- lm(as.formula(f_str), data = d)
  p <- autoplot(fit, which = 1:6, colour = 'dodgerblue3',
                smooth.colour = 'black', smooth.linetype = 'dashed',
                ad.colour = 'blue', 
                label.size = 3, label.n = 5, label.colour = 'blue',
                ncol = 3)
  print(p)
  print(paste("========= Model: ", f_str, " ========"))
  print(f)
  print(summary(fit))
  return(fit)
}

fits <- lapply(c('te', 'protein_len', 'gene_len', 'transcripts', 'mcs_ds'), base_model, te_data2)

### TODO: Compare models with/without wgc

