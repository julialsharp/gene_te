## Transposable Elements Analysis
source("te_common.R")

IN_DIR <- "./data"
OUT_DIR <- "./out_full"
APPLY_FILTER = FALSE
filtered_symbols <- c('TTN')

outpath <- function(f) {
  return(paste(OUT_DIR, "/", f, sep = '')) 
}

inpath <- function(f) {
  return(paste(IN_DIR, "/", f, sep = ''))
}


# read data
all_data <- read_csv(inpath("alldata_srini.csv"))

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
  filter(!(symbol %in% filtered_symbols & APPLY_FILTER)) %>% 
  select(gene_id, asd, devreg, cns, housekeeping, gene_len, protein_len, transcripts, te, mcs_ds, rvis)

# save a copy for later use
# write_csv(te_data, inpath("gene_te_data.csv"))

# check for duplicates in original data and its subset 
all_data %>% 
  group_by(symbol) %>% 
  filter(n()>1)
# Source: local data frame [4 x 27]
# Groups: symbol [2]

te_data %>% 
  group_by(gene_id) %>% 
  filter(n()>1)

# add wgc column - using the given definition
te_data2 <- te_data %>%
  mutate(asd = ifelse(asd == 'asd-linked', 'Yes', 'No'), 
         devreg = ifelse(devreg == 'Devreg', 'Yes', 'No'), 
         cns = ifelse(cns == 'CNS', 'Yes', 'No'), 
         housekeeping = ifelse(housekeeping == 'housekeeping', 'Yes', 'No')) %>% 
  mutate(wgc = ifelse(asd == 'Yes' | devreg == 'Yes' | cns == 'Yes' | housekeeping == 'Yes', 'No', 'Yes'))

te_data2_1 <- gather(te_data2, key = gene_type, value = gene_linked, c(2:5, 12), -c(gene_id))
te_data2_2 <- gather(te_data2, key = metric_type, value = metric, c(6:11), -c(gene_id))

# keep only data 'TRUE' for each gene type, i.e. filter out all rows where gene_type is not-linked
te_data2_3 <- inner_join(te_data2_2, te_data2_1, by = "gene_id") %>%
  filter(gene_linked != 'No') %>% 
  select(gene_id, gene_type, metric_type, metric, gene_linked)


# Historgrams - Grid
p <- ggplot(data = te_data2_3, aes(metric)) + 
  facet_grid(gene_type ~ metric_type, scale = "free_x") +
  geom_histogram() +
  scale_x_continuous(name = 'Gene Type', breaks = scales::pretty_breaks(n = 3)) + 
  scale_y_continuous(name = 'Frequency', breaks = scales::pretty_breaks(n = 10)) +
  ggtitle('Histogram of Measurements for Gene Types') + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) +
  scale_fill_brewer(palette = "Accent") ; p

ggsave(p, filename = outpath('summary_histogram.png'), width = 45, height = 20, units = "cm")

# Box plot - Grid
p <- ggplot(data = te_data2_3, aes(x = gene_type, y = metric)) + 
  facet_wrap(~ metric_type, scale = "free_y") +
  geom_boxplot() +
  scale_x_discrete(name = "Gene Type") +
  scale_y_continuous(name = "Measurement", breaks = scales::pretty_breaks(n = 15)) +
  ggtitle('Box plot of Measurements') + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 12, family = "Tahoma"),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11, angle = 60, hjust = 1)) +
  scale_fill_brewer(palette = "Accent"); p

ggsave(p, filename = outpath('summary_boxplot.png'), width = 30, height = 20, units = "cm")


p +
  geom_boxplot(outlier.shape = NA) +
  scale_y_continuous(limits = quantile(te_data2_3$metric, c(0.1, 0.9), na.rm = TRUE))

# create summary data with wgc as part the gene category
te_data2_sum <- te_data2_3 %>% 
  group_by(gene_type, metric_type) %>% 
  summarise_at(c('metric'), function(x) { round(mean(as.numeric(x), na.rm = TRUE), 2) })

# Bar chart
p1 = ggplot(data = te_data2_sum, aes(x = gene_type, y = metric) ) + 
  facet_grid(metric_type ~ ., scale = "free_y") +
  geom_bar(aes(fill = gene_type), stat = "identity") + 
  labs(x = 'Gene Type', y = 'Mean Value') +
  ggtitle('Summary Satistics') +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        legend.position="none", 
        plot.title = element_text(hjust = 0.5))

p1  

ggsave(p1, filename = outpath('summary_barchart.png'), width = 10, height = 20, units = "cm")

#
# descriptive stats
te_stats_1 <- te_data2_3 %>% 
  group_by(metric_type) %>% 
  summarize(mean(metric, na.rm = T), sd(metric, na.rm = T), median(metric, na.rm = T), 
            min(metric, na.rm = T), max(metric, na.rm = T))  

# using base R as dplyr rename is not working as expected!
names(te_stats_1)[2:6] <- c("mean_val", "sd_val", "median_val", "min_val", "max_val")
te_stats_1 <- te_stats_1 %>% mutate_if(is.numeric, funs(round(., 2)))
write_csv(te_stats_1, outpath('summary_stat.csv'))

te_stats_2 <- te_data2_3 %>% 
  group_by(gene_type, metric_type) %>% 
  summarize(mean(metric, na.rm = T), sd(metric, na.rm = T), median(metric, na.rm = T), 
            min(metric, na.rm = T), max(metric, na.rm = T))
  
# using base R as dplyr rename is not working as expected for changes using an index!
names(te_stats_2)[3:7] <- c("mean_val", "sd_val", "median_val", "min_val", "max_val")
te_stats_2 <- te_stats_2 %>% mutate_if(is.numeric, funs(round(., 2)))  


#
# Testing - Pairwise t-tests
t.test(te ~ asd, data = te_data)
t.test(protein_len ~ asd, data = te_data)
t.test(gene_len ~ asd, data = te_data)
t.test(transcripts ~ asd, data = te_data)
t.test(mcs_ds ~ asd, data = te_data)
