# 2B_Sens: Reconciling Overt Bias and Hidden Bias in Sensitivity Analysis for Matched Observational Studies

`2B_Sens` is an R package designed for conducting a Rosenbaum-type sensitivity analysis informed by post-matching overt bias information.

---

## Authors

- **Siyu Heng**  
  Department of Biostatistics, New York University, New York, U.S.A.  
  Email: siyuheng@nyu.edu

- **Yanxin Shen** (**Maintainer**)  
  Institute of Finance, School of Economics, Nankai University, Tianjin, China.  
  Email: shenyanxin5@163.com

- **Pengyun Wang**  
  Data Science Institute, The University of Chicago, Chicago, U.S.A.  
  Email: pengyunwang@uchicago.edu

---

## Features

- **Sensitivity Analysis**: Calculate the worst-case p-values under the Rosenbaum bounds constraint. 
- **Hypothesis Testing**: Conduct hypothesis testing under the Rosenbaum bounds constraint.
- **Flexible Model Options**: Support for various models including Random Forest, XGBoost, Logistic Regression, and SVM.
- **Customizable Analysis**: Parameters such as significance levels, sensitivity parameter Gamma, and machine learning model configurations.
- **Automated Workflow**: Combines preprocessing, sensitivity analysis, and hypothesis testing in a seamless pipeline.

---

## File Structure and Functions

- **maintest**: Complete test workflow combining preprocessing, sensitivity analysis, and hypothesis testing.
- **preprocess**: Preprocess the data by splitting it into two parts and formatting it for subsequent analysis.
- **gammatest**: Check whether the feasible set is empty.
- **testH**: If the feasible set is not empty, determine whether to accept or reject the null hypothesis.
- **p_val**: Full workflow from data preprocessing to p-value computation.
- **get_p**: Detailed process of calculating p-values using machine learning models.
- **gurobi**: Solve the optimization problems involved in the analysis.
- **utils**: Contains helper functions, including machine learning model settings and coefficient matrix calculations used during optimization.

---

## Installation

To install the package, use:

```R
# Install directly from the tar.gz file
devtools::install_github("YanxinShen/2B_Sens")