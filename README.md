# vaxstrat 

## Overview

**vaxstrat** implements principal stratification methods for estimating causal effects of vaccination on post-infection outcomes.

The package includes multiple estimators for evaluating vaccine effects in the following principal strata:

- the **Doomed**, individuals who would be infected regardless of vaccination),
- the **Naturally Infected**, individuals who would be infected only in the absence of a vaccine.

The package also includes estimators of nonparametric bounds for Naturally Infected effects, as well as efficient estimators of marginal effects.

## Installation

```r
# Install the development version from GitHub
# install.packages("devtools")
devtools::install_github("allicodi/vaxstrat")
```

## Key Features

- **Multiple estimands**: Naturally Infected, Doomed, and population-level (marginal) effects
- **Multiple estimators**: G-computation, IPW, AIPW, TMLE, nonparametric bounds, covariate-adjusted bounds, and sensitivity analyses
- **Flexible nuisance estimation**: Standard GLMs or SuperLearner ensembles 
- **Principled inference**: Closed-form standard errors for AIPW/TMLE; bootstrap for all other methods; permutation tests for bounds
- **Identification flexibility**: Support for the exclusion restriction, partial principal ignorability assumptions, both, or neither (bounds only)

## Conceptual Background

`vaxstrat` considers estimating effects of vaccines on post-infection outcomes using principal stratification. Under the assumption that the vaccine cannot cause an increased risk of infection, individuals can be partitioned based on their potential infection status under each treatment arm ($S(z)$ for $z = 0,1$).

| Stratum | $S(1)$ | $S(0)$ | Description |
|---|---|---|---|
| **Immune** | 0 | 0 | Never infected |
| **Protected** | 0 | 1 | Infected only in absence of vaccine |
| **Doomed** | 1 | 1 | Infected regardless of vaccine |

The **Naturally Infected** stratum is consists of the Protected and Doomed individuals, i.e., individuals who would be infected in absence of vaccination. 

Identification of effects in the Naturally Infected requires either an **exclusion restriction** (vaccination cannot affect the post-infection outcome in absence of an infection), **partial principal ignorability** (that Protected and Immune individuals are exchangeable conditional on a set of covariates $X$), or both. Without either assumption, nonparametric bounds are available.

## Example: Rotavirus Vaccine Trial Reanalysis
 
The package includes `provide`, a simulated dataset inspired by the PROVIDE study, which examined rotavirus vaccination effects on antibiotic use in infants. The variables map directly to the main function arguments:
 
| Variable | Role | Description |
|---|---|---|
| `rotaarm` | Treatment (`Z`) | Rotavirus vaccine assignment (0/1) |
| `rotaepi` | Infection (`S`) | Rotavirus infection episode (0/1) |
| `any_abx_wk52` | Outcome (`Y`) | Any antibiotic use by 52 weeks (0/1) |
| `wk10_haz` | Covariate (`X`) | Height-for-age Z-score at 10 weeks |
| `gender` | Covariate (`X`) | Infant gender |
| `num_hh_sleep` | Covariate (`X`) | Household size (sleeping members) |
 
```r
library(vaxstrat)
 
# Load the included PROVIDE-inspired dataset
data(provide)
 
fit <- vaxstrat(
  data     = provide,
  Y_name   = "any_abx_wk52",                              # antibiotic use by 52 weeks
  Z_name   = "rotaarm",                                   # rotavirus vaccine assignment
  S_name   = "rotaepi",                                   # rotavirus infection episode
  X_name   = c("wk10_haz", "gender", "num_hh_sleep"),     # baseline covariates
  estimand = c("nat_inf", "doomed", "pop"),               # naturally infected strata, doomed strata, population (marginal)
  method   = c("aipw"),                                   # aipw estimators for all 
  seed     = 12345
)
 
print(fit)
```

## Estimation Methods

Estimation of effects requires estimation of certain regression quantities. These include: 
- the conditional mean of the postinfection outcome ($Y$) given infection status ($S$), vaccine status ($Z$) and covariates ($X$),
- the conditional probability of infection ($S$) given vaccine ($Z$) and covariates ($X$), and
- the conditional probability of vaccination ($Z$) given covariates ($X$).

Based on (a subset of) these estimates, the following estimation techniques are available. 

| Method | Description | SE source |
|---|---|---|
| `gcomp` | G-computation (outcome regression) | Bootstrap |
| `ipw` | Inverse probability weighting (propensity score) | Bootstrap |
| `aipw` | Augmented IPW (multiply robust) | Closed-form or bootstrap |
| `tmle` | Targeted maximum likelihood estimation | Closed-form or bootstrap |
| `bound` | Nonparametric bounds (no cross-world assumptions) | Bootstrap / permutation |
| `cov_adj_bound` | Covariate-adjusted bounds | Bootstrap |
| `sens` | Sensitivity analysis over range of `epsilon` values | Closed-form or bootstrap |
 
AIPW and TMLE estimators are **multiply robust**: they remain consistent if certain combinations of models are correctly specified. 

When `ml = TRUE`, nuisance parameters are estimated using **SuperLearner** ensemble learning.
 
## Sensitivity Analysis
 
When `method = "sens"`, the package sweeps over a grid of sensitivity parameter values `epsilon` (deviations from the partial principal ignorability assumption) and reports point estimates, standard errors, and hypothesis test results at each value:
 
```r
fit_sens <- vaxstrat(
  data     = provide,
  Y_name   = "any_abx_wk52",
  Z_name   = "rotaarm",
  S_name   = "rotaepi",
  X_name   = c("wk10_haz", "gender", "num_hh_sleep"),
  estimand = "nat_inf",
  method   = "sens",
  epsilon  = exp(seq(log(0.5), log(2), length = 49))
)

plot.sens(fit_sens$nat_inf$sens)
```
 
## Output Structure
 
`vaxstrat()` returns an object of class `"vaxstrat"` with results organized by estimand and method:
 
```
fit
├── nat_inf
│   ├── aipw_ER_CW
│   │   ├── pt_est        # point estimates (additive and log-multiplicative)
│   │   ├── boot_se       # bootstrap SEs and CIs (if applicable)
│   │   ├── test_stat     # test statistics
│   │   ├── p_val         # p-values
│   │   └── reject        # rejection indicators at alpha_level
│   ├── tmle
│   ├── bound
│   │   └── permutation   # permutation test results (if permutation = TRUE)
│   └── sens
│       └── reject        # data frame of results across epsilon grid
├── doomed
├── pop
└── models                # fitted nuisance models (if return_models = TRUE)
```
 
Use `print(fit)` for a formatted table of results.
 

## Citation
 
If you use `vaxstrat` in your research, please cite our forthcoming manuscript.
