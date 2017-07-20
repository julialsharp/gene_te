##
## Classification models
##
source("te_common.R")
library(rpart)
library(rpart.plot)

# split data into two
set.seed(1)
train <- sample_frac(te_data2, 0.7)
sid <- as.numeric(rownames(train)) # because rownames() returns character
test <- te_data2[-sid,]

wgc_model <- rpart(wgc ~ gene_len + protein_len + transcripts + te, 
                   data = train, minsplit = 2, minbucket = 1)
rpart.plot(wgc_model)

asd_model <- rpart(asd ~ ., 
                   data = train, control = rpart.control(minsplit = 2))
rpart.plot(asd_model)
