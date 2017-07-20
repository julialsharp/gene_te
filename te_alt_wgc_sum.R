## Transposable Elements Analysis

library(tidyverse)

# read data
all_data <- read_csv("./data/alldata_srini.csv")

# rename columns (replace only defined while leaving undefined intact)
col_data <- read_csv("./data/col_names.csv")
original_names <- as.data.frame(names(all_data), stringsAsFactors = FALSE)
names(original_names) <- c('original_name')
revised_names <- original_names %>% 
  left_join(col_data) %>% 
  select(ifelse(revised_name, revised_name, original_name))
names(all_data) <- revised_names[,1]

# select relevant columns
te_data <- all_data %>% 
  select(symbol, asd, devreg, cns, housekeeping, gene_len, protein_len, transcripts, te, mcs_ds, rvis)

# save a copy for later use
# write_csv(te_data, "./data/gene_te_data.csv")

#
# simpler definition of wgc to avoid having to drop rows
# e.g. for asd genes, the section above removes all genes not linked to asd, but linked to others
# In this section, treat wgc flag locally for each gene_type

# reshape data
te_data_1 <- gather(te_data, key = gene_type, value = wgc, c(2:5), -c(symbol))
te_data_2 <- gather(te_data, key = metric_type, value = metric, c(6:11), -c(symbol))
te_data_3 <- inner_join(te_data_2, te_data_1, by = "symbol") %>%
  mutate(wgc = ifelse(wgc %in% c('asd-linked', 'Devreg', 'CNS', 'housekeeping'), 'No', 'Yes')) %>% 
  select(symbol, gene_type, wgc, metric_type, metric)

# create summary data
te_data_sum <- te_data_3 %>% 
  group_by(gene_type, wgc, metric_type) %>% 
  summarise_if(is.numeric, function(x) { round(mean(x, na.rm = TRUE), 2) })

# Graph
p2 = ggplot(data = te_data_sum, aes(x = gene_type, y = metric, group = wgc, colour = wgc) ) + 
  facet_grid(metric_type ~ ., scale = "free_y") +
  geom_line() + geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)); p2

# Testing
t.test(te ~ asd, data = te_data)

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
