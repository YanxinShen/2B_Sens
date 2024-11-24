#' Solve Quadratic Programming Problem Using Gurobi
#'
#' This function solves a quadratic programming problem using the Gurobi optimizer.
#' It is designed to evaluate hypotheses for sensitivity analysis in observational studies.
#'
#' @param Q Matrix, the quadratic coefficient matrix for the objective function.
#' @param w Numeric vector, the linear coefficients for the objective function.
#' @param A Matrix, the constraint matrix defining equality and inequality constraints.
#' @param vec Numeric vector, typically representing predicted probabilities or propensity scores.
#' @param N Integer, the total number of observations.
#' @param ns Integer, the total number of groups in the matching structure.
#' @param ns1 Integer, the total number of groups for \strong{one-to-one} and \strong{one-to-many} matching.
#' @param ns2 Integer, the total number of groups for \strong{many-to-one} matching.
#' @param upper Numeric vector, the upper bounds for the decision variables.
#' @param lower Numeric vector, the lower bounds for the decision variables.
#' @param tstar Numeric, the test statistic for the current group.
#' @param params List, additional parameters for the Gurobi optimization solver.
#' @param quadcon List, optional quadratic constraints for the optimization problem (default is \code{NULL}).
#' @param test Logical (default is \code{FALSE}), whether to perform iterative testing with increasing \code{NumericFocus} parameter.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{H0} - Hypothesis result (0 = accept null hypothesis, 1 = reject null hypothesis, \code{NA} if no solution is found).
#'   \item \code{calculation} - Final objective value, representing the test statistic.
#' }
#'
#' @details
#' The function constructs a quadratic programming model and solves it using the Gurobi optimizer.
#' The optimization problem includes:
#' \itemize{
#'   \item An objective function with quadratic (\code{Q}) and linear (\code{w}) terms.
#'   \item Constraints defined by the matrix \code{A} and the bounds \code{upper} and \code{lower}.
#'   \item Optional quadratic constraints (\code{quadcon}) for sensitivity analysis.
#' }
#'
#' The function supports iterative testing with increasing \code{NumericFocus} values to improve solver robustness when \code{test = TRUE}.
#'
#' The hypothesis result (\code{H0}) is determined based on:
#' \itemize{
#'   \item If \code{ns2 == 0}: \code{calculation = result$objval + tstar^2}.
#'   \item If \code{ns2 > 0}: \code{calculation = result$objval + sum(vec[(ns1 + 1):N])^2 + tstar^2 - 2 * tstar * sum(vec[(ns1 + 1):N])}.
#' }
#' \code{H0} is set to 0 if the calculation is less than or equal to 0, and 1 otherwise.
#'
#' @export






cal_gurobi = function(Q, w, A, vec, N, ns, ns1, ns2, upper, lower, tstar, params = list(), quadcon=NULL, test=FALSE){

  #### input ####
  set.seed(123)
  model = list()
  model$obj = w
  model$Q = Q
  model$vtype = c(rep("C", N))
  model$modelsense = "min"
  model$A = A
  model$lb = lower
  model$ub = upper
  model$rhs = c(rep(1, ns), rep(0, nrow(A) - ns))
  model$sense = c(rep('=', ns), rep("<=", nrow(A)- ns))

  if (!is.null(quadcon)){
    model$quadcon = quadcon
  }

  if (test==TRUE){
    i = 0
    while (i <= 3) {
      params$NumericFocus = i
      result = gurobi(model, params)
      if (!is.null(result$objval)){
        break
      }
      i = i + 1
    }
  }else{
    result = gurobi(model, params)
  }

  #### calculation ####
  result$objval
  if (result$status == "OPTIMAL" ||
      result$status == "SUBOPTIMAL") {
    # For the condition when ns2 is 0
    if (ns2 == 0) {
      calculation = result$objval + tstar ^ 2
      if (calculation <= 0) {
        H0 = 0
      } else if (calculation > 0) {
        H0 = 1
      }
      # For the condition when ns2 is greater than 0
    } else if (ns2 > 0) {
      calculation = result$objval + sum(vec[(ns1 + 1):N]) ^ 2 + (tstar) ^ 2 - 2 * tstar * sum(vec[(ns1 + 1):N])
      if (calculation < 0) {
        H0 = 0
      } else if (calculation > 0) {
        H0 = 1
      }
    }
  } else {
    print("No solution found or an error occurred.")
    calculation=NA
    H0 = NA
  }
  cat("Objective = ", calculation, "\n")
  return(list(H0, calculation))
}
