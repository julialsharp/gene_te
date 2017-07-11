# Transposable Elements Analysis

library(tidyverse)
library(data.table)

# read data
all_data <- read_csv("./data/alldata_srini.csv")
spec(genome)

# rename columns (replace only defined while leaving undefined intact)
col_data <- read_csv("./data/col_names.csv")
original_names <- as.data.frame(names(all_data), stringsAsFactors = FALSE)
names(original_names) <- c('original_name')
revised_names <- left_join(original_names, col_data) %>% 
  select(ifelse(revised_name, revised_name, original_name))
names(all_data) <- revised_names[,1]

# select relevant columns
te_data <- all_data %>% 
  select(symbol, asd, devreg, cns, housekeeping, gene_len, protein_len, transcripts, te, rvis)
  

# create a random subset
write_csv(te_data, "./data/gene_te_data.csv")


# add wgc column - given definition
te_data2 <- 
  te_data %>%
  mutate(wgc = ifelse(asd == 'asd-linked' | DevReg == 'Devreg' | cns == 'CNS' | housekeeping == 'housekeeping', 'No', 'Yes'))

# descriptive analysis
te_data_1 <- gather(te_data, key = "gene_type", value = "wgc", c(2:5), -c(symbol))
te_data_2 <- gather(te_data, key = type, value = metric, -c(symbol))
te_data_3 <- inner_join(te_data_2, te_data_1, by = "symbol") %>%
  select(symbol, type, metric, gene_type, wgc )

# create summary data
te_data_3_1 = as.data.table(te_data_3)
te_data_sum = te_data_3_1[, .(means.ave = round(mean(metric),2)), by = list(wgc, gene_type, type)]

# Graph
p2 = ggplot(data = long3.3, aes(x = gene_type, y = means.ave, group = wgc, colour = wgc) ) + 
  facet_grid(type ~ ., scale = "free_y") +
  geom_line() + geom_point() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)); p2

# Testing
model0 = lm(measurement ~ wgc, data = te_data_3) # one factor
summary(model0)

model1 = lm(measurement ~ wgc*type, data = te_data_3)# two factors with interaction
summary(model1) # you can also use update() but not as clear

model2 = lm(measurement ~ wgc *type*gene_type, data = te_data_3) # three factors w. interaction
summary(model2)

# use aov to get a summary for each factor 
mod = aov(lm(measurement ~ wgc*type*gene_type, data = te_data_3)) # not great
summary(mod)

# compare graphics and testing
