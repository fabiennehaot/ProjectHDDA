get_data <- function(seed = 42, training_pct= .7){
  set.seed(seed)
  data("Einecke2010Kidney")
  X <- Einecke2010Kidney[,-1]
  Y <- Einecke2010Kidney$Reject_Status
  
  ### Sample 70 pct random IDs from the rows of X (250 total)
  ### Non balanced dataset -> keep proportions of YN
  nr_training <- floor(as.numeric(table(Y))*.7)
  id_normal <- sample(which(Y==0),nr_training[1])
  id_reject <- sample(which(Y==1),nr_training[2])
  trainID <- sample(nrow(X), floor(training_pct*nrow(X)))
  
  ### Training data
  trainX <- X[trainID, ]
  trainY <- Y[trainID]
  trainX <- scale(trainX, center = TRUE, scale = TRUE)
  scale_factor <- attr(trainX,"scaled:scale")
  scale_translation <- attr(trainX,"scaled:center")
  trainX = as.data.frame(trainX)
  ### Test data: scale accordingly
  testX <- X[-trainID, ]
  testY <- Y[-trainID]
  for (varnr in seq_along(testX)){
    xbar <- scale_translation[varnr]
    xsd <- scale_factor[varnr]
    testX[[varnr]] <- (testX[[varnr]]-xbar )/xsd}
  
  list(tr_x=trainX,tr_y=trainY,test_x=testX,test_y=testY)
  }
