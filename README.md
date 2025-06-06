# jmBIG: Joint Longitudinal and Survival Models for Big Data
<img src="man/figures/jmBIGlogo.png" align="right" alt="jmBIG logo" width="200">

<!-- badges: start -->

[![CRAN Status](https://www.r-pkg.org/badges/version/jmBIG)](https://CRAN.R-project.org/package=jmBIG)
[![License: GPL-3](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/jmBIG)](https://cran.r-project.org/package=jmBIG)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

<br/>

> **`jmBIG`** is an R package that provides a flexible, modular framework for **joint modeling of longitudinal and time-to-event data**, specifically optimized for **large datasets** commonly found in real-world clinical and biomedical studies.

---

## ✨ Features

- 🧠 **Joint modeling** of longitudinal biomarkers and survival outcomes
- ⚡ **Big data readiness**: designed for efficiency on large sample sizes
- 🔁 Compatible with:
  - `JMbayes2` (Bayesian joint modeling)
  - `FastJM` (Frequentist joint modeling)
  - `rstanarm` (Bayesian GLMMs)
  - `joineRML` (MLE-based methods)
- 🧾 Posterior prediction of:
  - Individual survival probabilities
  - Longitudinal trajectories
- 📊 Tools for:
  - Dynamic predictions
  - Bootstrapped confidence intervals
  - Visual summaries

---

## 📦 Installation

### From GitHub (development version)

```r
# If not already installed
install.packages("remotes")

# Install jmBIG from GitHub
remotes::install_github("kumarbhrigu/jmBIG")
```
