## Transposable Elements Analysis

library(tidyverse)
library(data.table)

# read data
all_data <- read_csv("./data/alldata_srini.csv")

# rename columns (replace only defined while leaving undefined intact)
col_data <- read_csv("./data/col_names.csv")
original_names <- as.data.frame(names(all_data), stringsAsFactors = FALSE)
names(original_names) <- c('original_name')
revised_names <- left_join(original_names, col_data) %>% 
  select(ifelse(revised_name, revised_name, original_name))
names(all_data) <- revised_names[,1]

# select relevant columns
te_data <- all_data %>% 
  select(symbol, asd, devreg, cns, housekeeping, gene_len, protein_len, transcripts, te, mcs_ds, rvis)

# save a copy for later use
# write_csv(te_data, "./data/gene_te_data.csv")

# add wgc column - using given definition
te_data2 <- 
  te_data %>%
  mutate(asd = ifelse(asd == 'asd-linked', 'Yes', 'No'), 
         devreg = ifelse(devreg == 'Devreg', 'Yes', 'No'), 
         cns = ifelse(cns == 'CNS', 'Yes', 'No'), 
         housekeeping = ifelse(housekeeping == 'housekeeping', 'Yes', 'No')) %>% 
  mutate(wgc = ifelse(asd == 'Yes' | devreg == 'Yes' | cns == 'Yes' | housekeeping == 'Yes', 'No', 'Yes'))

te_data2_1 <- gather(te_data2, key = gene_type, value = gene_linked, c(2:5), -c(symbol))
te_data2_2 <- gather(te_data, key = metric_type, value = metric, c(6:11), -c(symbol)) # use original data with wgc so metric value stays numeric
te_data2_3 <- inner_join(te_data2_2, te_data2_1, by = "symbol") %>%
  filter(!(wgc == 'No' & gene_linked == 'No')) %>% 
  #filter(!metric_type == 'wgc') %>% # not needed anymore since metric does not include wgc
  select(symbol, metric_type, metric, gene_type, wgc )
  
# create summary data
te_data2_3_1 = as.data.table(te_data2_3)
te_data2_sum = te_data2_3_1[, .(mean_val = round(mean(metric, na.rm = TRUE),2)), by = list(wgc, gene_type, metric_type)]

# Graph
p1 = ggplot(data = te_data2_sum, aes(x = gene_type, y = mean_val, group = wgc, colour = wgc) ) + 
  facet_grid(metric_type ~ ., scale = "free_y") +
  geom_line() + geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)); p1  

# Testing
t.test(te ~ asd, data = te_data2)
t.test(te ~ wgc, data = te_data2)

  
#
# simpler definition of wgc to avoid having to drop rows
# e.g. for asd genes, the section above removes all genes not linked to asd, but linked to others
# In this section, treat wgc flag locally for each gene_type

# reshape data
te_data_1 <- gather(te_data, key = gene_type, value = wgc, c(2:5), -c(symbol))
te_data_2 <- gather(te_data, key = metric_type, value = metric, c(6:11), -c(symbol))
te_data_3 <- inner_join(te_data_2, te_data_1, by = "symbol") %>%
  mutate(wgc = ifelse(wgc %in% c('asd-linked', 'Devreg', 'CNS', 'housekeeping'), 'No', 'Yes')) %>% 
  select(symbol, metric_type, metric, gene_type, wgc )

# create summary data
te_data_3_1 = as.data.table(te_data_3)
te_data_sum = te_data_3_1[, .(mean_val = round(mean(metric, na.rm = TRUE),2)), by = list(wgc, gene_type, metric_type)]

# Graph
p2 = ggplot(data = te_data_sum, aes(x = gene_type, y = mean_val, group = wgc, colour = wgc) ) + 
  facet_grid(metric_type ~ ., scale = "free_y") +
  geom_line() + geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)); p2

# Testing
model0 = lm(metric ~ wgc, data = te_data_3) # one factor
summary(model0)

model1 = lm(metric ~ wgc*metric_type, data = te_data_3) # two factors with interaction
summary(model1) 

model2 = lm(metric ~ wgc*metric_type*gene_type, data = te_data_3) # three factors w. interaction
summary(model2)

# use aov to get a summary for each factor 
mod = aov(lm(metric ~ wgc*metric_type*gene_type, data = te_data_3)) # ?
summary(mod)

# compare graphs and tests
