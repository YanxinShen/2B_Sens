#' Gamma Sensitivity Analysis Test
#'
#' This function performs a gamma sensitivity analysis test to check whether the feasible set is empty.
#' The function computes test statistics and evaluates constraints using quadratic programming (e.g., Gurobi).
#'
#' @param d1 Data frame for the first group. It must contain the following variables:
#' \itemize{
#'   \item \code{q}: Numeric vector, representing the outcome variable for the first group.
#'   \item \code{phat}: Numeric vector, predicted probabilities of Z for the first group.
#'   \item \code{Z}: Numeric vector (0 or 1), indicating treatment assignment (1 = treated, 0 = control).
#'   \item \code{index}: Numeric vector, grouping index for observations in the first group.
#'   \item \code{upper}: Numeric vector, upper bounds of optimization constraints for decision variables.
#'   \item \code{lower}: Numeric vector, lower bounds of optimization constraints for decision variables.
#' }
#' @param d2 A list or data frame representing the second dataset with the same structure as \code{d1}.
#' @param alpha2 Numeric, significance level for the sensitivity analysis (default is 0.0125).
#' @param rank Logical, whether to apply a rank test on \code{q} (default is \code{FALSE}).
#' @param test_beta Numeric, a value to adjust the observed effect size (default is 0).
#' @param gamma Numeric, sensitivity parameter representing the degree of departure from random assignment (default is \code{gamma}).
#' @param params List, additional parameters for the Gurobi optimization solver.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{bigF} - Combined test result for both groups.
#'   \item \code{F0.1} - Test result for the first group.
#'   \item \code{F0.2} - Test result for the second group.
#' }
#'
#' @details
#' The data frames \code{d1} and \code{d2} are assumed to follow a specific structure for the \code{index} column:
#' \itemize{
#'   \item Each \code{index} corresponds to a group of matched individuals. Groups are ordered such that all observations
#'         within a group appear together in the data frame.
#'   \item The \strong{one-to-many} or \strong{one-to-one} groups are listed first.
#'   \item The \strong{many-to-one} groups are listed afterward.
#'   \item Within each group, observations with \code{Z == 1} (treatment) always appear before those with \code{Z == 0} (control).
#'     }
#' These conventions ensure consistent processing for sensitivity testing and constraint generation.
#'
#' @examples
#' # Example usage of the gamma_test function
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
#' # Set parameters for the gamma sensitivity analysis
#' alpha2 <- 0.0125  # Significance level for sensitivity analysis
#' rank <- FALSE     # Whether to apply rank transformation to q
#' test_beta <- 0    # Adjust observed treatment effect
#' gamma <- 4        # Sensitivity parameter
#' params <- list(TimeLimit = 1000)  # Gurobi solver parameters
#'
#' # Perform gamma sensitivity analysis
#' result <- gamma_test(
#'   d1 = d1,
#'   d2 = d2,
#'   alpha2 = alpha2,
#'   rank = rank,
#'   test_beta = test_beta,
#'   gamma = gamma,
#'   params = params
#' )
#'
#' # View the result
#' print(result)
#'
#' # Interpretation:
#' # - `bigF`: Combined gamma sensitivity test result for both groups.
#' # - `F0.1` and `F0.2`: Gamma sensitivity test results for the first and second groups, respectively.
#' #   - 0 indicates the null hypothesis is accepted.
#' #   - 1 indicates the null hypothesis is rejected.
#'
#' @export



gamma_test = function(d1, d2, alpha2=0.0125, rank=FALSE, test_beta=0, gamma=4, params=list()){
  bigF = 0
  for (num in c(1, 2)) {
    data <- get(paste0("d", num))
    # preprocess for test
    q_vec = data$q
    s_vec = data$phat
    Z = data$Z
    index = data$index
    upper = data$upper
    lower = data$lower

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
    tstar = sum(Z * s_vec)

    # rank or add test beta
    if (rank == TRUE) {
      q_vec = rank(q_vec) # Rank test
    }
    q_vec[Z == 1] = q_vec[Z == 1] - test_beta


    #### objective ####
    # quadratic
    Q1 = get_Q(s_vec, Z, index, c2)
    # linear
    w1 = get_w(s_vec, c2, tstar, ns1, ns2)

    #### constraint ####
    A = constr(index, gamma=gamma)

    #### gurobi ####

    F0 = cal_gurobi(Q1, w1, A, s_vec, N, ns, ns1, ns2, upper, lower, tstar, params=params)[[1]]

    if(is.na(F0)){
      F0 = 1
    }

    var_name = paste0("F0.", num)
    assign(var_name, F0)
    bigF = bigF + F0

    cat("Gamma Test H0: ", F0, "\n")
  }
  return(list(bigF, F0.1, F0.2))
}
