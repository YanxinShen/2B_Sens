#' Preprocessed Matched Dataset
#'
#' This dataset contains information related to the balance of covariates, index of treated and matched control subjects,
#' and the main dataset used in a treatment effect analysis. It is designed for use in evaluating the treatment and control group balance before and after matching, as well as analyzing treatment effects.
#'
#' @format A list containing the following objects:
#' \describe{
#'   \item{treated_lst}{A list of indices representing the subjects in the treated group.
#'   Each element of the list corresponds to a group of treated subjects and stores their respective indices,
#'   facilitating subsequent analyses.}
#'
#'   \item{control_lst}{A list of indices representing the matched control subjects for each treated group.
#'   This list helps to identify the control group subjects that have been matched to each treated subject, based on the chosen matching algorithm.}
#'
#'   \item{dataset}{A data frame containing the data used for treatment effect analysis, including covariates, treatment assignment, and outcome variables:
#'     \describe{
#'       \item{X}{Covariates: A set of variables representing baseline characteristics of the subjects. These are used to adjust for confounding effects and ensure the comparability between the treated and control groups.}
#'       \item{Z}{Treatment assignment: A numeric vector indicating whether a subject received the treatment (1 for treated, 0 for control).}
#'       \item{R}{Outcome variable: A numeric vector representing the outcome of interest, which could be a continuous or binary variable depending on the analysis.}
#'       \item{P}{Propensity score: A numeric vector representing the estimated probability that a subject received the treatment given their baseline covariates. This variable is included in the dataset but is not required for the final analysis as it primarily serves as an intermediate result in the matching process.}
#'       \item{R_t}{Treatment outcome: The observed outcome for subjects who received the treatment. This variable is included in the dataset but is not needed for the final analysis.}
#'       \item{R_c}{Control outcome: The observed outcome for control subjects. Similar to `R_t`, this variable is included in the dataset but is not necessary for the final analysis.}
#'     }
#'   }
#' }
#'
#' @details This dataset is used in the context of a matching analysis to estimate the causal effect of a treatment.
#' The `balance` table provides crucial information on how well the covariates are balanced between treated and control groups before and after matching.
#' The indices (`treated.subject.index` and `matched.control.subject.index`) facilitate matching operations,
#' while the `dataset` contains all necessary variables for estimating treatment effects.
#' While `P`, `R_t`, and `R_c` are included in the dataset for completeness, they are not used in the final analysis.
#' The primary variables of interest are `X`, `Z`, and `R`.
#'
#' @source Generated for illustrative purposes
"sample_data"
