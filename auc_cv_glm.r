# cross validation as a function on glm
#' @param trials subsets of pcs to be evaluated ([1:k])
#' @param data dataset containing a response variable Y
#' @param folds number of folds in the crossvalidation
auc_cv_glm <- function(trials, data, folds){
  full_model <- suppressWarnings(glm(Y ~ ., data = data, family = "binomial"))
  n_trials <- length(trials)
  error_auc <- rep(0,n_trials)
  comp <- 1
  for (k_try in trials){
    full_model_cv <- suppressWarnings(suppressMessages(
      cv.glm(
        data = data[1:(k_try+1)],  
        glmfit = full_model,
        cost = pROC::auc, 
        K = folds # note: specify the auc function (from pROC) without`()`!
      )))
    error_auc[comp] = full_model_cv$delta[1] 
    comp <- comp + 1
  }
  df_auc <- tibble(auc = error_auc,
                   nr_of_pc = trials)
}
