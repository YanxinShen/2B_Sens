# utils.R
# =============================================================================
# Utility functions for data preprocessing, machine learning, and sensitivity analysis.
# =============================================================================

# =============================================================================
# Data Processing Functions
# =============================================================================

#' Permutation-Based Group Splitting Function
#'
#' This function performs a permutation-based split of treated and control groups, shuffling and dividing
#' them into two balanced groups for further analysis. It ensures that single-element indices are grouped first.
#'
#' @param treated_lst List of indices for the treated group
#' @param control_lst List of indices for the control group
#' @return A list containing two grouped indices, each as a combination of treated and control indices
#' @export

permutation_split = function(treated_lst, control_lst){
  random_rows = sample(1:length(treated_lst), length(control_lst))
  treated.index.shuffle = treated_lst[random_rows]
  control.index.shuffle = control_lst[random_rows]
  half_nrow = length(treated_lst) %/% 2
  treated.index.group1 = treated.index.shuffle[1:half_nrow]
  control.index.group1 = control.index.shuffle[1:half_nrow]
  treated.index.group2 = treated.index.shuffle[(half_nrow + 1):length(treated_lst)]
  control.index.group2 = control.index.shuffle[(half_nrow + 1):length(treated_lst)]

  single_element_indices = which(sapply(treated.index.group1, function(x) is.numeric(x) && length(x) == 1))
  new_order = c(single_element_indices, setdiff(seq_along(treated.index.group1), single_element_indices))
  treated.index.group1 = treated.index.group1[new_order]
  control.index.group1 = control.index.group1[new_order]

  single_element_indices = which(sapply(treated.index.group2, function(x) is.numeric(x) && length(x) == 1))
  new_order = c(single_element_indices, setdiff(seq_along(treated.index.group2), single_element_indices))
  treated.index.group2 = treated.index.group2[new_order]
  control.index.group2 = control.index.group2[new_order]

  group1.index = mapply(c, treated.index.group1, control.index.group1, SIMPLIFY = FALSE)
  group2.index = mapply(c, treated.index.group2, control.index.group2, SIMPLIFY = FALSE)

  return(list(group1.index, group2.index))
}


# =============================================================================
# Machine Learning Prediction Functions
# =============================================================================

#' XGBoost Prediction Function
#'
#' @param Y Numeric vector of response variable for the training data (1 for treated, 0 for control)
#' @param x Data frame or matrix of predictor variables for training
#' @param x.2 Data frame or matrix of predictor variables for prediction
#' @param ... Additional parameters to pass to the XGBoost model
#' @return A numeric vector of predicted values (probabilities)
#' @export
m.xgboost <- function(Y, x, x.2, ...) {
  default_params <- list(
    max.depth = 6,
    eta = 0.1,
    nthread = 5,
    objective = "binary:logistic",
    eval_metric = "error"
  )
  user_params <- list(...)
  params <- modifyList(default_params, user_params)
  nrounds = ifelse(is.null(user_params$nrounds), 100, user_params$nrounds)
  bst = xgboost(data = as.matrix(x), label = Y, params = params, nrounds = nrounds)
  y.hat = predict(bst, as.matrix(x.2))
  return(y.hat)
}

#' Support Vector Machine (SVM) Prediction Function
#'
#' @param Y Response variable for training
#' @param x Predictors for training
#' @param x.2 Predictors for prediction
#' @param ... Additional SVM parameters
#' @return A vector of predicted values
#' @export
m.svm <- function(Y, x, x.2, ...) {
  default_params <- list(
    kernel = "radial"
  )
  user_params <- list(...)
  params <- modifyList(default_params, user_params)
  x$Z <- Y
  svm_model <- do.call(svm, c(list(formula = Z ~ ., data = x), params))
  y.hat <- RBSA:::predict(svm_model, newdata = x.2)

  return(as.vector(y.hat))
}

#' Random Forest Prediction Function
#'
#' @param Y Response variable for training
#' @param x Predictors for training
#' @param x.2 Predictors for prediction
#' @param ... Additional RF parameters
#' @return A vector of predicted probabilities
#' @export
m.rf <- function(Y, x, x.2, ...) {

  x$Z <- as.factor(Y)
  rf_model <- randomForest(Z ~ ., data = x, ntree = 100)
  y.hat <- RBSA:::predict(rf_model, newdata = x.2, type = "prob")

  return(as.vector(y.hat[, 2]))
}

#' Logistic Regression Prediction Function
#'
#' @param Y Response variable for training
#' @param x Predictors for training
#' @param x.2 Predictors for prediction
#' @param ... Additional GLM parameters
#' @return A vector of predicted probabilities
#' @export
m.logistic <- function(Y, x, x.2, ...) {
  default_params <- list(
    family = binomial
  )
  user_params <- list(...)
  params <- modifyList(default_params, user_params)
  x$Z <- Y
  model <- do.call(glm, c(list(formula = Z ~ ., data = x), params))
  y.hat <- RBSA:::predict(model, newdata = x.2, type = "response")

  return(as.vector(y.hat))
}

# =============================================================================
# Matrix and Vector Calculation
# =============================================================================

#' Generate Quadratic Coefficient Matrix for Quadratic Programming
#' @param vec Numeric vector
#' @param Z Treatment assignment vector
#' @param index Grouping index
#' @param c Chi-squared critical value
#' @return Symmetric matrix Q
#' @export
get_Q = function(vec, Z, index, c){
  N = length(index)
  ns = max(index)
  grouped_Z = tapply(Z, index, function(x) sum(x == 1) == 1)
  num_single_ones = sum(grouped_Z)
  ns1 = sum(index %in% 1:num_single_ones)
  if (num_single_ones == 0) {
    ns1 = 0
  }

  ns2 = N - ns1


  Q = outer(vec, vec, FUN = "*")
  for (i in (1:ns)) {
    iindex = which(index == i)
    Q[iindex, iindex] = (1 + c) * Q[iindex, iindex]
  }
  if (ns1 > 0 && ns2 > 0) {
    Q[1:ns1, (ns1 + 1):N] = -1 * Q[1:ns1, (ns1 + 1):N]
    Q[(ns1 + 1):N, 1:ns1] = -1 * Q[(ns1 + 1):N, 1:ns1]
  }
  return(Q)
}

#' Calculate Linear Coefficients for Quadratic Programming
#'
#' @param vec Numeric vector
#' @param c Chi-squared critical value
#' @param tstar Test statistic
#' @param ns1 Groups for one-to-one matching
#' @param ns2 Groups for many-to-one matching
#' @return Linear coefficients w
#' @export
get_w = function(vec, c, tstar, ns1, ns2){
  N = length(vec)
  I1 = vec[1:ns1]
  I2 = vec[(ns1 + 1):N]
  sum.I2 = sum(I2)
  if (is.na(sum.I2)) {
    sum.I2 = 0
  }
  w.1 = (2 * sum.I2 - 2 * tstar) * I1 - c * I1 ^ 2
  w.2 = (-2 * sum.I2 + 2 * tstar) * I2 - c * I2 ^ 2
  if (is.na(sum(w.2))) {
    w.2 = 0
  }
  w = c(w.1, w.2)
  if (ns1 == 0) {
    w = w.2
  } else if (ns2 == 0) {
    w = w.1
  }
  return(w)
}

#' Generate Constraint Matrix for Quadratic Programming
#' @param index Grouping index
#' @param gamma Sensitivity parameter
#' @return Constraint matrix A
#' @export

constr = function(index, gamma=4){

  N = length(index)
  ns = max(index)

  final_nrow = ns + sum(sapply(unique(index), function(x) sum(index == x) * (sum(index == x) - 1)))
  A = matrix(0, nrow = final_nrow, ncol = N)
  row_counter = ns + 1
  start = 1

  for (i in 1:ns) {
    iindex = which(index == i)
    A[i, iindex] = 1
    num = length(iindex)
    if (num > 1) {
      combinations = expand.grid(j=1:num, k=1:num)
      combinations = combinations[combinations$j != combinations$k, ]
      for (row in 1:nrow(combinations)) {
        j = combinations$j[row]
        k = combinations$k[row]
        A[row_counter, iindex[j]] = 1
        A[row_counter, iindex[k]] = -gamma
        row_counter = row_counter + 1
      }
    }
  }
  return(A)
}


# =============================================================================
# Prediction Wrapper
# =============================================================================

#' Prediction Wrapper Function
#' @param Y Response variable
#' @param x Predictors for training
#' @param x.2 Predictors for prediction
#' @param model Model type
#' @param ... Additional parameters
#' @return Predicted values
#' @export
get_predict = function(Y, x, x.2, model='rf', ...){
  func_name = paste0("m.", model)
  func <- tryCatch({
    match.fun(RBSA:::func_name)
  }, error = function(e) {
    stop(paste("Unsupported Model:", model, "- ensure that the function", func_name, "exists."))
  })
  Y.hat = func(Y, x, x.2, ...)
  return(Y.hat)
}

