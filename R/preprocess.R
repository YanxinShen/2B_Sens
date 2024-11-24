#' Data Preprocessing Function
#'
#' This function preprocesses the data for treated and control groups by grouping, feature extraction,
#' and machine learning prediction, generating the final preprocessed datasets.
#'
#' @param treated_lst List of indices for the treated group
#' @param control_lst List of indices for the control group
#' @param Z Numeric vector indicating treatment assignment (1 for treated, 0 for control)
#' @param X Data frame or matrix containing covariates
#' @param R Numeric vector containing outcome variables
#' @param gamma Numeric, sensitivity analysis parameter (default is 4)
#' @param predicts Optional, pre-computed predictions to bypass the model training step.
#' @param model Character string specifying the type of machine learning model to use
#'        (e.g., "xgboost", "logistic", "svm", or "rf")
#' @param ... Additional parameters passed to the model
#'
#' @return A list containing two data frames for the treated and control groups, each with various variables
#' @importFrom xgboost xgboost
#' @importFrom e1071 svm
#' @importFrom randomForest randomForest
#' @import dplyr
#'
#' @examples
#' # Extract variables from the dataset
#' treated_lst <- sample_data$treated_lst
#' control_lst <- sample_data$control_lst
#' Z <- sample_data$dataset$Z  # Treatment assignment
#' X <- sample_data$dataset$X  # Covariates
#' R <- sample_data$dataset$R  # Outcome variable
#'
#'
#' # Preprocess the data using the random forest model
#' result <- preprocess(treated_lst, control_lst, Z, X, R, gamma = 4, model = "rf", ntree = 50)
#'
#' # Access preprocessed data for treated and control groups
#' treated_data <- result[[1]]
#' control_data <- result[[2]]
#'
#' # Display the first few rows of the treated group data
#' head(treated_data)
#'
#' # Display the first few rows of the control group data
#' head(control_data)
#'
#' @export

preprocess = function(treated_lst, control_lst, Z, X, R, gamma=4, model='rf', predicts=NULL, ...){
  idx = permutation_split(treated_lst, control_lst)
  group1.index = idx[[1]]
  group2.index = idx[[2]]

  #generate input
  lengths.of.groups1 = sapply(group1.index, length)
  Z.group1 = c()
  index.group1 = c()
  x.group1 <- X[0, ]
  q.group1 = c()
  s.group1 = c()
  upper.group1 = c()
  lower.group1 = c()

  for (i in 1:length(group1.index)) {
    index.group1 = c(index.group1, rep(i, times = lengths.of.groups1[i]))
    upper.group1 = c(upper.group1, rep(gamma/(gamma + (lengths.of.groups1[i] - 1)), times = lengths.of.groups1[i]))
    lower.group1 = c(lower.group1, rep(1/(1 + gamma * (lengths.of.groups1[i] - 1)), times = lengths.of.groups1[i]))
    Z.group1 = c(Z.group1, Z[unlist(group1.index[i])])
    q.group1 = c(q.group1, R[unlist(group1.index[i])])
    if (!is.null(predicts)){
      s.group1 = c(s.group1, predicts[unlist(group1.index[i])])
    }
    for (k in unlist(group1.index[i])) {
      x.group1 = rbind(x.group1, X[k, ])
    }
  }


  lengths.of.groups2 = sapply(group2.index, length)
  Z.group2 = c()
  index.group2 = c()
  x.group2 <- X[0, ]
  q.group2 = c()
  s.group2 = c()
  upper.group2 = c()
  lower.group2 = c()

  for (i in 1:length(group2.index)) {
    index.group2 = c(index.group2, rep(i, times = lengths.of.groups2[i]))
    upper.group2 = c(upper.group2, rep(gamma/(gamma + (lengths.of.groups2[i] - 1)), times = lengths.of.groups2[i]))
    lower.group2 = c(lower.group2, rep(1/(1 + gamma * (lengths.of.groups2[i] - 1)), times = lengths.of.groups2[i]))
    Z.group2 = c(Z.group2, Z[unlist(group2.index[i])])
    q.group2 = c(q.group2, R[unlist(group2.index[i])])
    if (!is.null(predicts)){
      s.group2 = c(s.group2, predicts[unlist(group2.index[i])])
    }
    for (k in unlist(group2.index[i])) {
      x.group2 = rbind(x.group2, X[k, ])
    }
  }

  x.group1 <- as.data.frame(x.group1)
  x.group2 <- as.data.frame(x.group2)


  #machine learning
  if (is.null(predicts)){
    if(model == 'rf'){
      x.group1$Z <- as.factor(Z.group1)
      x.group2$Z <- as.factor(Z.group2)
      default_params <- list(
        ntree = 100,
        nodesize = 5,
        maxnodes = NULL
      )
      user_params <- list(...)
      params <- modifyList(default_params, user_params)
      rf_model.group2 <- randomForest(Z ~ ., data = x.group2, ntree = params$ntree)
      s.group1 <- predict(rf_model.group2, newdata = x.group1, type = "prob")
      s.group1 <- as.vector(s.group1[, 2])

      rf_model.group1 <- randomForest(Z ~ ., data = x.group1, ntree = params$ntree)
      s.group2 <- predict(rf_model.group1, newdata = x.group2, type = "prob")
      s.group2 <- as.vector(s.group2[, 2])
    }else{

      s.group2 = get_predict(Z.group1, x.group1, x.group2, model=model, ...)
      s.group1 = get_predict(Z.group2, x.group2, x.group1, model=model, ...)
    }
  }


  # generate phat
  s.group1[s.group1 == 1] = 0.9999
  s.group2[s.group2 == 1] = 0.9999
  s.group1[s.group1 == 0] = 0.0001
  s.group2[s.group2 == 0] = 0.0001
  phat.group2 = c()
  phat.group1 = c()
  # group2
  for (i in (1:length(group2.index))) {
    e.sum = sum(s.group2[index.group2 == i])
    match_type = sum(Z.group2[index.group2 == i])
    n_i = length(group2.index[[i]])
    phat.prod = 0
    if (match_type > 1) {
      for (j in (1:n_i)) {
        phat.prod =
          phat.prod + prod(s.group2[index.group2 == i]) / (s.group2[index.group2 == i][j]) * (1 - s.group2[index.group2 == i][j])
      }
    } else{
      for (j in (1:n_i)) {
        phat.prod =
          phat.prod + prod(1 - s.group2[index.group2 == i]) / (1 - s.group2[index.group2 == i][j]) * s.group2[index.group2 == i][j]
      }
    }


    if (match_type > 1) {
      for (j in (1:n_i)) {
        phat.group2 =
          c(phat.group2,
            prod(s.group2[index.group2 == i]) / (s.group2[index.group2 == i][j]) * (1 - s.group2[index.group2 == i][j]) / phat.prod)
      }
    } else{
      for (j in (1:n_i)) {
        phat.group2 =
          c(phat.group2,
            prod(1 - s.group2[index.group2 == i]) / (1 - s.group2[index.group2 == i][j]) * s.group2[index.group2 == i][j] / phat.prod)
      }
    }
  }

  #group1
  for (i in (1:length(group1.index))) {
    e.sum = sum(s.group1[index.group1 == i])
    match_type = sum(Z.group1[index.group1 == i])
    n_i = length(group1.index[[i]])
    phat.prod = 0
    if (match_type > 1) {
      for (j in (1:n_i)) {
        phat.prod =
          phat.prod + prod(s.group1[index.group1 == i]) / (s.group1[index.group1 == i][j]) * (1 - s.group1[index.group1 == i][j])
      }
    } else{
      for (j in (1:n_i)) {
        phat.prod =
          phat.prod + prod(1 - s.group1[index.group1 == i]) / (1 - s.group1[index.group1 == i][j]) * s.group1[index.group1 == i][j]
      }
    }


    if (match_type > 1) {
      for (j in (1:n_i)) {
        phat.group1 =
          c(phat.group1,
            prod(s.group1[index.group1 == i]) / (s.group1[index.group1 == i][j]) * (1 - s.group1[index.group1 == i][j]) / phat.prod)
      }
    } else{
      for (j in (1:n_i)) {
        phat.group1 =
          c(phat.group1,
            prod(1 - s.group1[index.group1 == i]) / (1 - s.group1[index.group1 == i][j]) * s.group1[index.group1 == i][j] / phat.prod)
      }
    }
  }
  data.group1 = data.frame(Z.group1, q.group1, phat.group1, s.group1, index.group1, upper.group1, lower.group1)
  colnames(data.group1) <- c("Z", "q", "phat", "s", "index", "upper", "lower")
  data.group2 = data.frame(Z.group2, q.group2, phat.group2, s.group2, index.group2, upper.group2, lower.group2)
  colnames(data.group2) <- c("Z", "q", "phat", "s", "index", "upper", "lower")

  return(list(data.group1, data.group2))
}
