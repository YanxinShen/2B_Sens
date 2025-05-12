#' Perform Hypothesis Testing
#'
#' This function conducts hypothesis testing for sensitivity analysis in observational studies.
#' It evaluates the null hypothesis (\code{H0}) for both treated and control groups by solving
#' quadratic programming problems with constraints.
#'
#' @param d1 Data frame for the first group. It must contain the following variables:
#' \itemize{
#'   \item \code{q}: Numeric vector, outcome variable for the group.
#'   \item \code{phat}: Numeric vector, predicted probabilities or propensity scores for the group.
#'   \item \code{Z}: Numeric vector (0 or 1), treatment assignment (1 = treated, 0 = control).
#'   \item \code{index}: Numeric vector, grouping index for matched individuals.
#'   \item \code{upper}: Numeric vector, upper bounds for decision variables.
#'   \item \code{lower}: Numeric vector, lower bounds for decision variables.
#' }
#' @param d2 Data frame for the second group. It has the same structure as \code{d1}.
#' @param F0.1 Integer (default is 0), gamma test result for the first group.
#' Set to 1 if the first group hypothesis was already rejected in a prior step.
#' @param F0.2 Integer (default is 0), gamma test result for the second group.
#' Set to 1 if the second group hypothesis was already rejected in a prior step.
#' @param alpha1 Numeric (default is 0.0125), significance level for the main hypothesis testing step.
#' @param alpha2 Numeric (default is NA), significance level for sensitivity analysis.
#' If \code{conv_test = TRUE}, \code{alpha2} is not used.
#' @param rank Logical (default is \code{FALSE}), whether to apply a rank transformation to \code{q}.
#' @param test_beta Numeric (default is 0), value to adjust the observed treatment effect.
#' @param gamma Numeric, sensitivity parameter representing the degree of departure from random assignment.
#' @param conv_test Logical (default is \code{TRUE}), whether to perform a conventional sensitivity analysis.
#' If \code{FALSE}, quadratic constraints are added for advanced sensitivity testing.
#' @param params List of additional parameters for the Gurobi optimization solver.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{bigH} - Combined hypothesis result for both groups (0 = accept null hypothesis, 1 = reject null hypothesis).
#'   \item \code{H0.1} - Hypothesis result for the first group.
#'   \item \code{H0.2} - Hypothesis result for the second group.
#'   \item \code{obj.1} - Objective value for the first group.
#'   \item \code{obj.2} - Objective value for the second group.
#' }
#'
#' @details
#' The function performs hypothesis testing for two groups (\code{d1} and \code{d2}).
#' It constructs quadratic coefficient matrices and linear terms using helper functions
#' (\code{get_Q} and \code{get_w}). Constraints are built using \code{constr}, and
#' optimization is performed with the Gurobi solver. The function can optionally apply
#' advanced sensitivity testing by incorporating quadratic constraints (\code{quadcon}).
#'
#' The testing procedure involves:
#' \itemize{
#'   \item Constructing and solving a quadratic programming problem for each group.
#'   \item Evaluating the hypothesis result (\code{H0}) based on the solver's output.
#'   \item Combining the results (\code{H0.1}, \code{H0.2}) into a final decision (\code{bigH}).
#' }
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
#' result <- Htest(
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
#' # Interpretation:
#' # - `bigH`: Combined hypothesis result for both groups (0 = accept, 1 = reject).
#' # - `H0.1` and `H0.2`: Hypothesis results for the first and second groups.
#' # - `obj.1` and `obj.2`: Objective values for the first and second groups.
#'
#' @export


Htest = function(d1, d2, F0.1=0, F0.2=0, alpha1=0.0125, alpha2=NA, rank=FALSE, test_beta=0, gamma=gamma, conv_test=TRUE, params=list()){
  bigH=0
  H0.1 = NA
  H0.2 = NA
  obj.1 = NA
  obj.2 = NA
  for (num in c(1, 2)) {
    data <- get(paste0("d", num))
    # preprocess
    F0 = get(paste0("F0.", num))
    if (F0 == 1){
      next
    }
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

    #### gurobi ####
    res = cal_gurobi(Q, w, A, q_vec, N, ns, ns1, ns2, upper, lower, t, quadcon=quadcon, test=TRUE, params=params)
    H0 = res[[1]]


    var_name = paste0("H0.", num)
    assign(var_name, H0)

    if (!is.na(H0)){
      bigH = bigH + H0
    }


    obj = res[[2]]
    var_name = paste0("obj.", num)
    assign(var_name, obj)

    cat("H0: ", H0, "\n")

  }

  return(list(bigH, H0.1, H0.2, obj.1, obj.2))

}
