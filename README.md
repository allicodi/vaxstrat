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
  method   = c("aipw", "gcomp", "bound"),                 # estimators
  return_se = TRUE,                                       # closed form standard error for AIPW
  n_boot = 100,                                           # small number of bootstrap replicates gcomp standard error and bounds for quick example
  seed     = 12345,
  family = "binomial"
)
 
print(fit)

                         Growth Effect Estimation Results: Additive
------------------------------------------------------------------------------------------ 
Estimand                 Method         Point Est.     95% CI: Lower  95% CI: Upper  
------------------------------------------------------------------------------------------ 
Naturally Infected - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation (Cross-World only)         -0.1589        -0.2023        -0.1148        
                         Lower Bound                              -0.5621         -0.6434        -0.4890        
                         Upper Bound                              0.0339         -0.0047        0.0767         
                         AIPW (Cross-World only)                  -0.1685        -0.2256        -0.1113        
Naturally Infected (ER) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation (ER + Cross-World)         -0.1605        -0.2072        -0.1174        
                         G-Computation (ER only)                  -0.1606        -0.2951        -0.0443        
                         AIPW (ER + Cross-World)                  -0.1687        -0.2232        -0.1141        
                         AIPW (ER only)                           -0.1606        -0.2900        -0.0313        
Doomed - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation  -0.0307        -0.0853        0.0229         
                         Lower Bound    -0.1000        -0.1655        -0.0516        
                         Upper Bound    0.0829         -0.0125        0.2009         
                         AIPW           -0.0267        -0.0879        0.0346         
Population - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation  -0.0647        -0.1171        -0.0174        
                         AIPW           -0.0647        -0.1178        -0.0117    
                         
print(fit, "mult")
                         Growth Effect Estimation Results: Multiplicative
------------------------------------------------------------------------------------------ 
Estimand                 Method         Point Est.     95% CI: Lower  95% CI: Upper  
------------------------------------------------------------------------------------------ 
Naturally Infected - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation (Cross-World only)0.8275         0.7816         0.8726         
                         Lower Bound    0.3928         0.3087         0.4652         
                         Upper Bound    1.0366         0.9951         1.0861         
                         AIPW (Cross-World only)0.8184         0.7638         0.8769         
Naturally Infected (ER) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation (ER + Cross-World)0.8259         0.7828         0.8688         
                         G-Computation (ER only)0.8257         0.6840         0.9522         
                         AIPW (ER + Cross-World)0.8182         0.7666         0.8732         
                         AIPW (ER only) 0.8268         0.6995         0.9774         
Doomed - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation  0.9674         0.9085         1.0247         
                         Lower Bound    0.9000         0.8345         0.9484         
                         Upper Bound    1.1015         0.9860         1.2821         
                         AIPW           0.9714         0.9080         1.0392         
Population - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                         G-Computation  0.9147         0.8521         0.9765         
                         AIPW           0.9147         0.8497         0.9845         
                         
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
 
`vaxstrat()` returns an object of class `"vaxstrat"` with results organized by estimand and method. This example shows results from `fit` above, but the structure follows for all estimand and estimator combinations.
 
```
fit
├── nat_inf               # object of class "nat_inf" containing all estimates for naturally infected estimand
│   ├── gcomp_ER          # object of class "gcomp_ER" (g-computation with exclusion restriction assumption only)
│   │   ├── pt_est        # point estimates (named numeric with additive effect, log multiplicative effect, psi_1, psi_0 )
│   │   ├── boot_se       # bootstrap SEs and CIs (data.frame with bootstrap standard error and 95% confidence interval for the point estimates above)
│   │   ├── test_stat     # test statistics (additive and multiplicative scales)
│   │   |   ├── additive
│   │   |   ├── multiplicative 
│   │   ├── p_val         # p-values (additive and multiplicative scales)
│   │   |   ├── additive
│   │   |   ├── multiplicative 
│   │   └── reject        # rejection indicators at alpha_level (additive and multiplicative scales)
│   │   |   ├── additive
│   │   |   ├── multiplicative 
│   ├── gcomp_CW          # same results structure as gcomp_ER, but for g-computation with cross-world assumption only
│   │   ├── ...
│   ├── gcomp_ER_CW       # same results structure as gcomp_ER and gcomp_CW, but for g-computation with cross-world and exclusion restriction assumptions
│   │   ├── ...
│   ├── aipw_ER           # object of class "aipw_ER" (augmented inverse probability weighted estimator with exclusion restriction only)
│   |   ├── pt_est        # named numeric vector with same arguments as gcomp. if closed form standard error specified (return_SE = TRUE), will also include additive_se, log_multiplicative_se, se_psi_1, and se_psi_0 in addition to point estimates.
│   |   ├── boot_se       # will only be returned if return_SE = FALSE
│   │   ├── ...           # all other arguments same as gcomp  
│   ├── aipw_CW
│   │   ├── ...           # same format as aipw_ER
│   ├── aipw_ER_CW
│   │   ├── ...           # same format as aipw_ER
│   ├── bound
│   │   ├── pt_est        # named numeric vector with E_Y0__S0_1, E_Y1__S0_1_lower, E_Y1__S0_1_upper, additive_effect_lower, additive_effect_upper, mult_effect_lower, mult_effect_upper 
│   |   ├── boot_se       # data.frame with bootstrap standard error and confidence intervals for all of the above
│   │   ├── ...           # all other arguments same as gcomp  
├── doomed                # object of class "doomed" containing all estimates for doomed estimand
│   ├── gcomp             # object of class "gcomp" 
│   │   ├── ...           # same nested structure as nat_inf estimator
│   ├── aipw              
│   │   ├── ...           # same nested structure as nat_inf estimator
│   ├── bound        
│   │   ├── pt_est        # named numeric vector with E_Y1__S0_1, E_Y0__S0_1_lower, E_Y0__S0_1_upper, additive_effect_lower, additive_effect_upper, mult_effect_lower, mult_effect_upper
│   │   ├── ...           # remainder of nested structure as nat_inf estimator
├── pop                   # object of class "pop" containing all estimates for population level estimand
│   ├── gcomp             
│   │   ├── ...           # same nested structure as nat_inf estimator
│   ├── aipw              
│   │   ├── ...           # same nested structure as nat_inf estimator
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
