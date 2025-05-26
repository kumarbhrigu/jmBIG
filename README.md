# 📦 jmBIG

> **Joint Longitudinal and Survival Models for Big Data**

[![CRAN Status](https://www.r-pkg.org/badges/version/jmBIG)](https://CRAN.R-project.org/package=jmBIG)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)

**`jmBIG`** is an R package providing efficient and flexible tools for joint modeling of **longitudinal** and **time-to-event** data, tailored for **large-scale ("BIG") datasets**.

It is designed especially for applications in **medical research**, enabling the analysis of complex relationships between repeated biomarker measurements and clinical events, such as survival or disease progression.

---

## ✨ Features

- 📈 **Joint modeling** of longitudinal and survival outcomes
- ⚡ Optimized for **Big Data**: handles large sample sizes efficiently
- 🔁 Compatible with:
  - `JMbayes2` (Bayesian)
  - `FastJM` (Frequentist)
  - `rstanarm` (Bayesian)
  - `joineRML` (Frequentist)
- 📊 Prediction of:
  - Survival probabilities
  - Longitudinal trajectories
- 🛠 Tools for:
  - Posterior prediction
  - Visualization
  - Bootstrapped confidence intervals

---

## 🛠 Installation

```r
# Install from CRAN (when available)
install.packages("jmBIG")

# Or install the development version from GitHub:
# install.packages("remotes")
remotes::install_github("kumarbhrigu/jmBIG")
