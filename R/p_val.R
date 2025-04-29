#' Calculate P-Value with Gamma Test and Sensitivity Analysis
#'
#' This function calculates the p-value for a given treated and control dataset
#' using a combination of a gamma test and sensitivity analysis. It processes the data,
#' applies optimization techniques, and determines whether to reject the null hypothesis.
#'
#' @param treated_lst A list containing data for the treated group.
#' @param control_lst A list containing data for the control group.
#' @param Z A binary vector indicating treatment assignment.
#' @param X A matrix or data frame of covariates.
#' @param R A vector of outcomes.
#' @param test_beta A numeric value used to adjust \code{q_vec} when \code{Z == 1}. Default is 0.
#' @param rank Logical. If \code{TRUE}, ranks the \code{q_vec} values. Default is \code{FALSE}.
#' @param alpha A numeric value specifying the overall significance level. Default is 0.05.
#' @param alpha1 A numeric value specifying the first stage significance level. Default is \code{NA}.
#' @param alpha2 A numeric value specifying the second stage significance level. Default is 0.025.
#' @param gamma A numeric value specifying the gamma parameter for constraints. Default is 4.
#' @param conv_test Logical. If \code{TRUE}, performs a convergence test. Default is \code{FALSE}.
#' @param model A character string specifying the prediction model to use. Options include \code{"rf"}, \code{"svm"}, etc. Default is \code{"rf"}.
#' @param predicts Optional pre-computed predictions to use for the analysis. Default is \code{NULL}.
#' @param seed An integer specifying the random seed for reproducibility. Default is \code{NULL}.
#' @param gurobi_params1 A list of parameters to pass to the Gurobi solver for the gamma test. Default is an empty list.
#' @param gurobi_params2 A list of parameters to pass to the Gurobi solver for the sensitivity analysis. Default is an empty list.
#' @param tol A numeric value specifying the tolerance for optimization. Default is \code{1e-5}.
#' @param verbose Logical. If \code{TRUE}, prints additional details during execution. Default is \code{TRUE}.
#' @param ... Additional arguments passed to the preprocessing or prediction functions.
#' @param nbreak An integer specifying the maximum number of iterations allowed before terminating the process in case of errors. Default is 10.
#'
#' @return A data frame with the following columns:
#'   - \code{obj1}: Objective value for group 1.
#'   - \code{obj2}: Objective value for group 2.
#'   - \code{group1}: Result of the null hypothesis test for group 1.
#'   - \code{group2}: Result of the null hypothesis test for group 2.
#'   - \code{F0.group1}: Gamma test result for group 1.
#'   - \code{F0.group2}: Gamma test result for group 2.
#'   - \code{bigH}: Combined result of the sensitivity analysis.
#'   - \code{pval}: Combined p-value from the analysis.
#'   - \code{pval.group1}: P-value for group 1.
#'   - \code{pval.group2}: P-value for group 2.
#'
#' @details
#' This function first preprocesses the treated and control datasets using the \code{preprocess} function.
#' It then performs a gamma test using \code{gamma_test} to calculate quadratic programming objectives
#' and constraints. Finally, it conducts a sensitivity analysis using \code{get_p} to compute p-values
#' and determine whether to reject or accept the null hypothesis.
#'
#' The significance levels (\code{alpha}, \code{alpha1}, \code{alpha2}) are adjusted dynamically based on the input.
#' The function assumes that required helper functions such as \code{preprocess}, \code{gamma_test},
#' and \code{get_p} are defined elsewhere and available in the environment.
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
#' result <- pvalue(
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
#'   gurobi_params2 = list(TimeLimit = 1000),
#'   tol=1e-5,
#'   verbose=TRUE,
#'   nbreak=10
#' )
#'
#' # View the result
#' print(result)
#'
#' @export


pvalue =
  function(treated_lst, control_lst, Z, X, R, test_beta = 0, rank = FALSE, alpha = 0.05, alpha1 = NA, alpha2 = 0.025,
           gamma = 4, conv_test = FALSE, model='rf', predicts=NULL, seed=NULL, gurobi_params1=list(), gurobi_params2=list(), tol=1e-5, verbose=TRUE, nbreak=10, ...) {

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
      main.results = get_p(d.group1, d.group2, F0.1=F0.group1, F0.2=F0.group2, alpha1=alpha1, alpha2=alpha2, rank=rank,
                           test_beta=test_beta, gamma=gamma, conv_test=conv_test, params=gurobi_params2, verbose=verbose, tol=tol, nbreak=nbreak)
    }

    pval.group1 = main.results[[1]]
    pval.group2 = main.results[[2]]

    pval = 2 * min(pval.group1, pval.group2, na.rm = TRUE)

    cat("P value=", pval, "\n")


    return(
      data.frame(
        pval=pval,
        pval.group1=pval.group1,
        pval.group2=pval.group2
      )
    )
  }
