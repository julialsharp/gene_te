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

# add wgc column - using given definition
te_data2 <- te_data %>%
  mutate(asd = ifelse(asd == 'asd-linked', 'Yes', 'No'), 
         devreg = ifelse(devreg == 'Devreg', 'Yes', 'No'), 
         cns = ifelse(cns == 'CNS', 'Yes', 'No'), 
         housekeeping = ifelse(housekeeping == 'housekeeping', 'Yes', 'No')) %>% 
  mutate(wgc = ifelse(asd == 'Yes' | devreg == 'Yes' | cns == 'Yes' | housekeeping == 'Yes', 'No', 'Yes'))

te_data2_1 <- gather(te_data2, key = gene_type, value = gene_linked, c(2:5), -c(symbol))
te_data2_2 <- gather(te_data2, key = metric_type, value = metric, c(6:12), -c(symbol))
te_data2_3 <- inner_join(te_data2_2, te_data2_1, by = "symbol") %>%
  filter(!(wgc == 'No' & gene_linked == 'No')) %>% 
  filter(metric_type != 'wgc') %>% # the column is not numeric due to wgc
  select(symbol, gene_type, wgc, metric_type, metric)
  
# create summary data
te_data2_sum <- te_data2_3 %>% 
  group_by(gene_type, wgc, metric_type) %>% 
  summarise_at(c('metric'), function(x) { round(mean(as.numeric(x), na.rm = TRUE), 2) })

# Graph
p1 = ggplot(data = te_data2_sum, aes(x = gene_type, y = metric, group = wgc, colour = wgc) ) + 
  facet_grid(metric_type ~ ., scale = "free_y") +
  geom_line() + geom_point() + 
  labs(x = 'Gene Type', y = 'Mean Value') +
  ggtitle('Summary Satistics') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)); p1  

ggsave(p1, filename = './output/summary.png')

# Testing
t.test(te ~ asd, data = te_data2)
t.test(te ~ wgc, data = te_data2)

# Testing
model0 = lm(metric ~ wgc, data = te_data2_3) # one factor
summary(model0)

model1 = lm(metric ~ wgc*metric_type, data = te_data2_3) # two factors with interaction
summary(model1) 

model2 = lm(metric ~ wgc*metric_type*gene_type, data = te_data2_3) # three factors w. interaction
summary(model2)

# use aov to get a summary for each factor 
mod = aov(lm(metric ~ wgc*metric_type*gene_type, data = te_data2_3)) # ?
summary(mod)

# compare graphs and tests
