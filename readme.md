# 2B_Sens: Reconciling Overt Bias and Hidden Bias in Sensitivity Analysis for Matched Observational Studies

`2B_Sens` is an R package designed for conducting a Rosenbaum-type sensitivity analysis informed by post-matching overt bias information.

---

## Features

- **Sensitivity Analysis**: Calculate the worst-case p-values under the Rosenbaum bounds constraint. 
- **Hypothesis Testing**: Conduct hypothesis testing under the Rosenbaum bounds constraint.
- **Flexible Model Options**: Support for various models including Random Forest, XGBoost, Logistic Regression, and SVM.
- **Customizable Analysis**: Parameters such as significance levels, sensitivity parameter Gamma, and machine learning model configurations.
- **Automated Workflow**: Combines preprocessing, sensitivity analysis, and hypothesis testing in a seamless pipeline.

---

## Installation

To install the package, use:

```R
# Install directly from the tar.gz file
devtools::install_github("YanxinShen/RBSA")
