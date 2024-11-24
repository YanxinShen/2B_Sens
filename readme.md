# RBSA: Reconciling Overt Bias and Hidden Bias in Sensitivity Analysis

`RBSA` is an R package designed for robust sensitivity analysis and hypothesis testing in observational studies. It incorporates advanced statistical methods, including gamma sensitivity tests and quadratic programming, to evaluate treatment effects under potential biases.

---

## Features

- **Sensitivity Analysis**: Perform gamma tests to assess the robustness of results to unobserved confounding.
- **Hypothesis Testing**: Evaluate treatment effects using quadratic programming and machine learning-assisted propensity score modeling.
- **Flexible Model Options**: Support for various models including Random Forest, XGBoost, Logistic Regression, and SVM.
- **Customizable Analysis**: Fine-tune parameters such as significance levels, gamma sensitivity, and machine learning model configurations.
- **Automated Workflow**: Combines preprocessing, sensitivity analysis, and hypothesis testing in a seamless pipeline.

---

## Installation

To install the package, use:

```R
# Install directly from the tar.gz file
devtools::install_github("YanxinShen/RBSA")
