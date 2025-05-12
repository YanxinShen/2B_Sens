#' Main Hypothesis Testing Function
#'
#' This function performs hypothesis testing between treated and matched control groups.
#' It includes data preprocessing, sensitivity analysis, and testing
#' based on specified parameters and model types.
#'
#' @param treated_lst List of indices for treated subjects
#' @param control_lst List of indices for matched control subjects
#' @param Z Numeric vector indicating treatment assignment (1 for treated, 0 for control)
#' @param X Data frame containing covariates
#' @param R Numeric vector containing outcome variables
#' @param test_beta Numeric value representing the treatment effect to test against (default is 0)
#' @param rank Logical, whether to apply a rank test on `q_vec` (default is FALSE)
#' @param alpha Numeric, significance level for the main test (default is 0.05)
#' @param alpha1 Numeric, significance level for the first stage gamma test (default is NA)
#' @param alpha2 Numeric, significance level for the second stage sensitivity analysis (default is 0.025)
#' @param gamma Numeric, sensitivity analysis parameter (default is 4)
#' @param conv_test Logical, whether to perform conventional sensitivity analysis (default is FALSE)
#' @param model Character string specifying the machine learning model used for propensity score estimation. Options include:
#'   \itemize{
#'     \item \code{"rf"} - Random forest (default).
#'     \item \code{"xgboost"} - Extreme Gradient Boosting.
#'     \item \code{"logistic"} - Logistic regression.
#'     \item \code{"svm"} - Support vector machine.
#'   }
#' @param predicts Optional, pre-computed predictions to bypass the model training step.
#' @param seed Numeric, optional random seed for reproducibility.
#' @param gurobi_params1 List of additional parameters for the Gurobi optimizer during the gamma test.
#' @param gurobi_params2 List of additional parameters for the Gurobi optimizer during the main hypothesis test.
#' @param ... Additional parameters passed to the machine learning model.
#'
#' @return A data frame with test results and sensitivity analysis metrics, including:
#' \itemize{
#'   \item \code{obj1} - Objective value for group 1
#'   \item \code{obj2} - Objective value for group 2
#'   \item \code{group1} - Hypothesis result for group 1 (accept(0) or reject(1))
#'   \item \code{group2} - Hypothesis result for group 2 (accept(0) or reject(1))
#'   \item \code{F0.group1} - Stage 1 result for group 1 (accept(0) or reject(1))
#'   \item \code{F0.group2} - Stage 1 result for group 2 (accept(0) or reject(1))
#'   \item \code{bigH} - Overall hypothesis decision (0 for acceptance, 1 for rejection)
#' }
#'
#' @importFrom gurobi gurobi
#' @importFrom Matrix Matrix
#' @importFrom matrixStats rowSums2
#'
#' @details
#' The function executes the following steps:
#' \itemize{
#'   \item Preprocess data using the \code{\link{preprocess}} function, generating datasets for tests.
#'   \item Perform a gamma test using the \code{\link{gamma_test}} function to evaluate sensitivity.
#'   \item If the gamma test passes, proceed to the main hypothesis test using the \code{\link{Htest}} function.
#'   \item Combine results from both tests to generate a final decision.
#' }
#' The significance levels (\code{alpha1} and \code{alpha2}) are automatically adjusted if not provided.
#'
#' @examples
#' # Example usage of the maintest function
#'
#' # Load a preprocessed dataset (replace 'sample_data' with your actual data object)
#' data(sample_data)
#'
#' # Extract variables from the dataset
#' treated_lst <- sample_data$treated_lst
#' control_lst <- sample_data$control_lst
#' Z <- sample_data$dataset$Z  # Treatment assignment
#' X <- sample_data$dataset$X  # Covariates
#' R <- sample_data$dataset$R  # Outcome variable
#'
#' # Perform the main hypothesis test
#' result <- maintest(
#'   treated_lst = treated_lst,
#'   control_lst = control_lst,
#'   Z = Z,
#'   X = X,
#'   R = R,
#'   test_beta = 0,
#'   rank = FALSE,
#'   alpha = 0.05,
#'   gamma = 4,
#'   model = "rf",  # Random forest
#'   seed = 123,
#'   gurobi_params1 = list(TimeLimit = 1000),
#'   gurobi_params2 = list(TimeLimit = 1000)
#' )
#'
#' # View the result
#' print(result)
#'
#' # Interpretation:
#' # - `obj1` and `obj2`: Objective values for groups 1 and 2.
#' # - `group1` and `group2`: Hypothesis test results for groups 1 and 2.
#' # - `bigH`: Overall hypothesis decision (0 = accept null hypothesis, 1 = reject null hypothesis).
#'
#' @export

# switch conv_test to TRUE if you want to implement conventional sensitivity analysis

maintest =
  function(treated_lst, control_lst, Z, X, R, test_beta = 0, rank = FALSE, alpha = 0.05, alpha1 = NA, alpha2 = 0.025,
           gamma = 4, conv_test = FALSE, model='rf', predicts=NULL, seed=NULL, gurobi_params1=list(), gurobi_params2=list(), ...) {

    set.seed(seed)
    alpha = alpha / 2
    alpha1 = alpha1 / 2
    alpha2 = alpha2 / 2
    if (!(is.na(alpha2))) {
      alpha1 <- (alpha - alpha2) / (1 - alpha2)
    } else if (!(is.na(alpha1))) {
      alpha2 <- (alpha - alpha1) / (1 - alpha1)
    } else{
      alpha2 <- alpha / 2
      alpha1 <- (alpha - alpha2) / (1 - alpha2)
    }
    if (alpha1 < 0 || alpha2 < 0) {
      stop("alpha1 or alpha2 is smaller than zero")
    }


    dts = preprocess(treated_lst, control_lst, Z, X, R, gamma=gamma, model=model, predicts=predicts, ...)
    d.group1 = dts[[1]]
    d.group2 = dts[[2]]

    #### test gamma ####
    gamma.results = gamma_test(d.group1, d.group2, alpha2=alpha2, rank=rank,
                               test_beta=test_beta, gamma=gamma, params=gurobi_params1)
    bigF=gamma.results[[1]]
    F0.group1 = gamma.results[[2]]
    F0.group2 = gamma.results[[3]]




    if (F < 2){
      main.results = Htest(d.group1, d.group2, F0.1=F0.group1, F0.2=F0.group2, alpha1=alpha1, alpha2=alpha2, rank=rank,
                              test_beta=test_beta, gamma=gamma, conv_test=conv_test, params=gurobi_params2)
    }
    bigH=main.results[[1]]
    H0.group1 = main.results[[2]]
    H0.group2 = main.results[[3]]
    obj.group1 = main.results[[4]]
    obj.group2 = main.results[[5]]


    ##############group1 and group2 complete#############


    if(is.na(H0.group1) && is.na(H0.group2)){
      bigH = NA
      print("Unable to reject or accpet the null hypothesis due to failure in sensitivity analysis or obtaining the test results.")
    }else{
      if(bigH > 0){
        bigH = 1
        print("Reject")
      }else{
        print("Accept")
      }
    }

    return(
      data.frame(
        obj1 = obj.group1,
        obj2 = obj.group2,
        group1 = H0.group1,
        group2 = H0.group2,
        F0.group1,
        F0.group2,
        bigH
      )
    )
  }



