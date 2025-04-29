#' Perform Sensitivity Analysis and Compute P-Values
#'
#' This function performs sensitivity analysis and computes p-values for two datasets, using
#' quadratic programming optimization and hypothesis testing. The results include objectives,
#' p-values, and hypothesis test results for each dataset.
#'
#' @param d1 A list or data frame representing the first dataset. It should contain:
#'   - \code{q}: A vector of quantile values.
#'   - \code{phat}: A vector of predicted probabilities.
#'   - \code{Z}: A binary vector indicating observed responses.
#'   - \code{index}: A grouping index for the data.
#'   - \code{upper}: Upper bounds for constraints.
#'   - \code{lower}: Lower bounds for constraints.
#' @param d2 A list or data frame representing the second dataset, with the same structure as \code{d1}.
#' @param F0.1 A numeric value representing the initial objective for dataset 1. Default is 0.
#' @param F0.2 A numeric value representing the initial objective for dataset 2. Default is 0.
#' @param alpha1 A numeric value specifying the first-stage significance level. Default is 0.0125.
#' @param alpha2 A numeric value specifying the second-stage significance level. Default is \code{NA}.
#' @param rank Logical. If \code{TRUE}, ranks the \code{q_vec} values. Default is \code{FALSE}.
#' @param test_beta A numeric value used to adjust \code{q_vec} when \code{Z == 1}. Default is 0.
#' @param gamma A numeric value specifying the gamma parameter for constraints.
#' @param conv_test Logical. If \code{TRUE}, performs a convergence test. Default is \code{TRUE}.
#' @param params A list of parameters to pass to the Gurobi solver.
#' @param tol A numeric value specifying the tolerance for optimization convergence. Default is \code{1e-5}.
#' @param verbose Logical. If \code{TRUE}, prints detailed progress information. Default is \code{TRUE}.
#' @param nbreak An integer specifying the maximum number of iterations allowed before terminating the process in case of errors. Default is 10.
#'
#'
#' @return A list containing:
#'   - \code{bigH}: Combined result of the hypothesis test.
#'   - \code{H0.1}: Hypothesis test result for dataset 1.
#'   - \code{H0.2}: Hypothesis test result for dataset 2.
#'   - \code{obj.1}: Objective value for dataset 1.
#'   - \code{obj.2}: Objective value for dataset 2.
#'   - \code{pvalue.1}: P-value for dataset 1.
#'   - \code{pvalue.2}: P-value for dataset 2.
#'
#' @details
#' The function iteratively adjusts the chi-squared critical value \code{c1} until convergence, using
#' quadratic programming to compute objectives and constraints. Sensitivity analysis is performed to test
#' the null hypothesis for each dataset, and p-values are computed based on the resulting chi-squared statistics.
#'
#' The function assumes that the required helper functions such as \code{get_Q}, \code{get_w}, \code{constr},
#' and \code{cal_gurobi} are defined elsewhere and available in the environment.
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
#' # Preprocess
#' dts = preprocess(treated_lst, control_lst, Z, X, R)
#' d1 = dts[[1]]
#' d2 = dts[[2]]
#'
#' # Set parameters for hypothesis testing
#' F0.1 <- 0  # Gamma test result for the first group
#' F0.2 <- 0  # Gamma test result for the second group
#' alpha1 <- 0.0125  # Significance level for the first stage
#' alpha2 <- 0.025   # Significance level for the second stage
#' rank <- FALSE     # Whether to use rank transformation
#' test_beta <- 0    # Adjust observed treatment effect
#' gamma <- 4        # Sensitivity parameter
#' conv_test <- FALSE # Use main test
#' params <- list(TimeLimit = 1000)  # Gurobi solver parameters
#'
#' # Perform hypothesis testing
#' result <- get_p(
#'   d1 = d1,
#'   d2 = d2,
#'   F0.1 = F0.1,
#'   F0.2 = F0.2,
#'   alpha1 = alpha1,
#'   alpha2 = alpha2,
#'   rank = rank,
#'   test_beta = test_beta,
#'   gamma = gamma,
#'   conv_test = conv_test,
#'   params = params
#' )
#'
#' # View the result
#' print(result)
#'
#' @export


get_p = function(d1, d2, F0.1=0, F0.2=0, alpha1=0.0125, alpha2=NA, rank=FALSE, test_beta=0, gamma=gamma, conv_test=TRUE, params=list(), tol=1e-5, verbose=TRUE, nbreak = 10){
  bigH=0
  H0.1 = NA
  H0.2 = NA
  obj.1 = NA
  obj.2 = NA
  pvalue.1 = NA
  pvalue.2 = NA
  ERR=1
  for (num in c(1, 2)) {

    data <- get(paste0("d", num))
    # preprocess
    F0 = get(paste0("F0.", num))
    if (F0 == 1){
      next
    }
    true_p = NA
    cold = -Inf
    dt = get(paste0("d", num))
    q_vec = dt$q
    s_vec = dt$phat
    Z = dt$Z
    index = dt$index
    upper = dt$upper
    lower = dt$lower

    c1 = qchisq((1 - alpha1), df = 1)
    c2 = qchisq((1 - alpha2), df = 1)
    N = length(index)
    ns = max(index)
    grouped_Z = tapply(Z, index, function(x) sum(x == 1) == 1)
    num_single_ones = sum(grouped_Z)
    ns1 = sum(index %in% 1:num_single_ones)
    if (num_single_ones == 0) {
      ns1 = 0
    }

    ns2 = N - ns1
    t = sum(Z * q_vec)
    tstar = sum(Z * s_vec)

    # rank or add test beta
    if (rank == TRUE) {
      q_vec = rank(q_vec) # Rank test
    }
    q_vec[Z == 1] = q_vec[Z == 1] - test_beta


    while (abs(cold - c1) > tol) {

      if (verbose==TRUE){
        cat('--------Group ', num, ': Try chi^2= ', c1,"--------\n")
      }
      #### objective ####
      # quadratic
      Q = get_Q(q_vec, Z, index, c1)
      # linear
      w = get_w(q_vec, c1, t, ns1, ns2)


      #### quadratic constraint ####
      # quadratic
      Q1 = get_Q(s_vec, Z, index, c2)
      # linear
      w1 = get_w(s_vec, c2, tstar, ns1, ns2)

      #### constraint ####
      A = constr(index, gamma=gamma)
      quadcon = list()

      if (alpha2 != 0 & conv_test==FALSE){
        if (ns2 > 0) {
          quadcon[[1]] = list(
            Qc = Q1,
            q = w1,
            rhs = -sum(s_vec[(ns1 + 1):N]) ^ 2 - (tstar) ^ 2 + 2 * tstar * sum(s_vec[(ns1 +
                                                                                        1):N]),
            sense = "<="
          )
        } else {
          quadcon[[1]] = list(
            Qc = Q1,
            q = w1,
            rhs = -(tstar) ^ 2,
            sense = "<="
          )
        }
      }
      params$OutputFlag = 0
      #### gurobi ####
      res = cal_gurobi(Q, w, A, q_vec, N, ns, ns1, ns2, upper, lower, t, quadcon=quadcon, test=TRUE, params=params)
      H0 = res[[1]]
      if (!is.na(H0)){
        p = res[[3]]
        pp = p
        if (ns2 > 0){
          p[(ns1 + 1):N] = 1 - p[(ns1 + 1):N]
        }

        exp = p * q_vec
        ET <- tapply(exp, index, sum)
        var1 <- tapply(pp * q_vec^2, index, sum)
        var2 <- tapply(pp * q_vec, index, sum)^2
        VART = var1 - var2
        cold = c1
        c1 = (t - sum(ET))^2 / sum(VART)

      }else{
        ERR = ERR + 1
        if (ERR > nbreak){
          ERR=NA
          break
        }
      }
    }

    if (is.na(ERR)){
      true_p = NA
    }else{
      true_p = 2 * (1 - alpha2) * (1 - pnorm(sqrt(c1))) + alpha2
    }
    var_name <- paste0("pvalue.", num)
    assign(var_name, true_p)

  }

  return(list(pvalue.1, pvalue.2))

}
