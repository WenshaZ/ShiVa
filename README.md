# ShiVa: Detection of Evolutionary Shifts in both Optimal Value and Variance

**Version:** 1.0.1  
**Authors:**  
- Wensha Zhang <wn209685@dal.ca>  
- Lam Si Tung Ho <Lam.Ho@dal.ca>  
- Toby Kenney <tb432381@dal.ca>  

## Overview

**ShiVa** is an R package for detecting evolutionary shifts in both optimal trait values (mean) and diffusion variance across phylogenetic trees, based on the Ornstein-Uhlenbeck (OU) model. While many existing methods detect shifts in optimal values only, ShiVa simultaneously identifies both types of shifts using a penalized likelihood framework with L1 regularization.

The method is designed for researchers in evolutionary biology and comparative methods, particularly those interested in identifying abrupt regime changes along a phylogeny.

See our preprint for methodological details: [arXiv:2312.17480](https://arxiv.org/abs/2312.17480)

## Features

- Automatively estimates shifts in both mean (θ) and variance (σ²) under the OU model
- Supports model selection using BIC, mBIC, or pBIC
- Includes backward correction for improved parsimony
- Can optionally account for measurement error
- Fast optimization using coordinate descent and soft-thresholding

## Installation

To install the package from GitHub:

```r
# install.packages("devtools") # if not already installed
devtools::install_github("WenshaZ/ShiVa")
```
