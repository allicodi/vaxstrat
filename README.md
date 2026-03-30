# vaxstrat 

## Overview

**vaxstrat** implements principal stratification methods for estimating causal effects of vaccination on post-infection outcomes.

Existing approaches focus on the **Doomed** stratum (individuals who would be infected regardless of vaccination), which may understate vaccine benefit by excluding individuals whose post-infection outcomes would have been improved had vaccination prevented infection entirely. `vaxstrat` additionally targets estimands for the **Naturally Infected** — individuals who would be infected only in the absence of a vaccine — and provides a suite of estimation strategies with accompanying inferential tools.

## Installation

```r
# Install the development version from GitHub
# install.packages("devtools")
devtools::install_github("allicodi/vaxstrat")
```

## Key Features

- **Multiple estimands**: Naturally Infected, Doomed stratum, and population-level (marginal) effects
- **Multiple estimators**: G-computation, IPW, AIPW, TMLE, nonparametric bounds, covariate-adjusted bounds, and sensitivity analyses
- **Flexible nuisance estimation**: Standard GLMs or SuperLearner ensembles (`ml = TRUE`)
- **Principled inference**: Closed-form standard errors for AIPW/TMLE; bootstrap for all other methods; permutation tests for bounds
- **Identification flexibility**: Support for the exclusion restriction, cross-world (partial principal ignorability) assumptions, or neither (bounds only)

## Conceptual Background

The key challenge is that post-infection outcomes (e.g., disease severity, antibiotic use, growth faltering) are only observed in infected individuals, and infection status is itself affected by treatment. Failure to account for this may introduce selection bias. 

`vaxstrat` handles this via principal stratification, partitioning individuals by their potential infection status under each treatment arm:

| Stratum | $S(1)$ | $S(0)$ | Description |
|---|---|---|---|
| **Immune** | 0 | 0 | Never infected |
| **Naturally Infected** | 0 | 1 | Infected in absence of vaccine |
| **Doomed** | 1 | 1 | Infected regardless of vaccine |

The **Naturally Infected** stratum is of particular interest: it captures individuals for whom vaccination prevented infection, meaning any post-infection outcomes observed under vaccination in this group reflect vaccine-modified disease. 

## Quick Start

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

## Estimands

### Naturally Infected (`nat_inf`)

The average treatment effect on the outcome among individuals who would be infected in the absence of vaccination ($S(0) = 1$). Requires either the **exclusion restriction** (vaccination does not affect the outcome except through preventing infection) or **partial principal ignorability** (cross-world assumption) for point identification. Without these, nonparametric bounds are available.

### Doomed (`doomed`)

The average treatment effect among individuals who would be infected regardless of vaccination ($S(0) = S(1) = 1$). This is the classical principal stratification estimand for post-infection outcomes.

### Population (`pop`)

The marginal average treatment effect on the outcome across all individuals.

## Estimation Methods
 
| Method | Description | SE source |
|---|---|---|
| `gcomp` | G-computation (outcome regression) | Bootstrap |
| `ipw` | Inverse probability weighting | Bootstrap |
| `aipw` | Augmented IPW (doubly robust) | Closed-form or bootstrap |
| `tmle` | Targeted maximum likelihood estimation | Closed-form or bootstrap |
| `bound` | Nonparametric bounds (no cross-world assumptions) | Bootstrap / permutation |
| `cov_adj_bound` | Covariate-adjusted bounds | Bootstrap |
| `sens` | Sensitivity analysis over range of `epsilon` values | Closed-form or bootstrap |
 
AIPW and TMLE estimators are **doubly robust**: they remain consistent if either the outcome model or the propensity/infection model is correctly specified. When `ml = TRUE`, nuisance parameters are estimated using **SuperLearner** ensemble learning.
 
## Identifying Assumptions
 
The function exposes two key assumptions, both of which can be toggled:
 
- **`exclusion_restriction`**: Vaccination affects the outcome only by preventing infection. Under this assumption, the Naturally Infected estimand simplifies because vaccinated individuals in this stratum are never observed to be infected.
- **`cross_world`** (partial principal ignorability): Conditional on covariates, potential infection status is independent of the potential outcome under the opposite treatment arm. This enables point identification without the exclusion restriction.
 
When neither assumption is imposed, the package returns **nonparametric bounds** on the estimand.
 
## Sensitivity Analysis
 
When `method = "sens"`, the package sweeps over a grid of sensitivity parameter values `epsilon` (deviations from the cross-world assumption) and reports point estimates, standard errors, and hypothesis test results at each value:
 
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
 
## Citation
 
If you use `vaxstrat` in your research, please cite:
 
> [Authors]. (Year). *[Paper title]*. [Journal]. [DOI]
