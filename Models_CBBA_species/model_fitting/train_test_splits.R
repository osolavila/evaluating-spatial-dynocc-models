Dir.Base <- "~/evaluating-spatial-dynamic-occupancy-models"
setwd(Dir.Base)
output.path <- file.path(Dir.Base, "colext", "Analysis", "Rdata",
                         "revision", "cross_validation")

set.seed(3)

#Randomly shuffle the data
shuffled_cells_id <- sample(nrow(cells_study))

#Create 5 equally size folds
folds <- cut(seq(1,nrow(cells_study)),breaks=5,labels=FALSE)


test_ids <- list()
train_ids <- list()
#Perform 5 fold cross validation
for(i in 1:5){
  testcut <- which(folds==i)
  test_ids[[i]] <- shuffled_cells_id[testcut]
  train_ids[[i]] <- shuffled_cells_id[-testcut]
}


save(test_ids, train_ids, file=file.path(output.path, "train_test_splits.Rdata"))