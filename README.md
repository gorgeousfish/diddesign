# diddesign

**Double Difference-in-Differences for Stata**

[![Stata 16+](https://img.shields.io/badge/Stata-16%2B-blue.svg)](https://www.stata.com/)
[![License: AGPL-3.0](https://img.shields.io/badge/License-AGPL--3.0-blue.svg)](LICENSE)
[![Version: 0.1.0](https://img.shields.io/badge/Version-0.1.0-green.svg)]()

![diddesign](image/image.png)

## Overview

`diddesign` implements the **Double Difference-in-Differences** (Double DID) method proposed by Egami and Yamauchi (2023, *Political Analysis*) for Stata. The method uses the generalized method of moments (GMM) to combine two estimators:

- **Standard DID**: Consistent under the parallel trends assumption
- **Sequential DID (s-DID)**: Consistent under the weaker parallel trends-in-trends assumption

When multiple pre-treatment periods are available, researchers can exploit three uses:

1. **Assessing the parallel trends assumption**: Pre-treatment periods enable placebo tests by applying the DID estimator to pre-treatment outcomes. The package also reports equivalence confidence intervals, which provide positive evidence for approximate parallel trends when the interval is narrow.

2. **Improving estimation accuracy**: Under the extended parallel trends assumption, multiple pre-treatment periods provide additional information to estimate counterfactual outcomes, reducing variance.

3. **Allowing for a more flexible parallel trends assumption**: The sequential DID estimator requires only the parallel trends-in-trends assumption, which permits linear time-varying confounding. When the extended parallel trends assumption is violated but parallel trends-in-trends holds, the Double DID converges to the sequential DID.

The Double DID combines these uses within a unified GMM framework. In panel data with multiple pre-treatment periods, the two-way fixed effects (TWFE) estimator is numerically equivalent to a weighted average of DIDs with equal weights. The Double DID generalizes this by selecting optimal weights via GMM.

**Features**:

- Double DID estimation via GMM
- Staggered adoption (SA) design support
- Parallel trends diagnostics with equivalence confidence intervals
- Visualization for trends, placebo tests, and treatment effects

## Key Concepts

### Identifying Assumptions

The Double DID framework relies on different identifying assumptions:

```
Extended Parallel Trends
        ↓ implies both
   ┌────┴────┐
   ↓         ↓
Standard    Parallel
Parallel    Trends-in-
Trends      Trends
```

Standard parallel trends and parallel trends-in-trends are distinct assumptions that do not imply each other.

| Assumption | Definition | Bias Over Time |
|------------|------------|----------------|
| **Extended Parallel Trends** | Parallel trends holds for all adjacent time periods | Constant |
| **Standard Parallel Trends** | Outcome trends are equal from last pre-treatment to post-treatment | Constant |
| **Parallel Trends-in-Trends** | The *change* in outcome trends is equal across groups | Linear change permitted |

Note: "Bias" refers to the difference in mean potential outcomes between treatment and control groups (selection bias). Under standard/extended parallel trends, this bias is constant over time. Under parallel trends-in-trends, bias can change at a constant rate.

### Estimator-Assumption Mapping

| Estimator | Required Assumption | Notes |
|-----------|---------------------|-------|
| **Standard DID** | Standard Parallel Trends | Uses one pre-treatment and one post-treatment period |
| **Extended DID** | Extended Parallel Trends | Uses all pre-treatment periods; equivalent to TWFE in panel data |
| **Sequential DID (s-DID)** | Parallel Trends-in-Trends | Requires at least two pre-treatment periods |
| **Double DID** | Extended (for efficiency) or Trends-in-Trends (for consistency) | GMM-optimal weighted combination of DID and s-DID |

### When to Use Each Approach

| Diagnostic Result | Recommendation |
|-------------------|----------------|
| Placebo close to zero, narrow equivalence CI | Use **Double-DID** |
| Placebo non-zero but trends have same direction | Consider **s-DID** |
| Trends have opposite directions | Caution: neither assumption may hold |

### Equivalence Confidence Intervals

Traditional hypothesis testing interprets non-rejection as evidence for parallel trends, but this conflates "no evidence against" with "evidence for." The equivalence approach reverses the null hypothesis: rejection of the null (that trends differ by more than some threshold) provides positive evidence for approximate parallel trends.

The 95% standardized equivalence CI is reported in units of baseline control group standard deviations. Narrower intervals provide stronger evidence for approximate parallel trends. Researchers should interpret the width of the equivalence CI in the context of their specific application, combining it with domain knowledge about what constitutes a substantively meaningful deviation for the outcome of interest.

## Requirements

- Stata 16.0 or later
- No additional dependencies (Mata included in Stata)

## Installation

### From GitHub

```stata
net install diddesign, from("https://raw.githubusercontent.com/gorgeousfish/diddesign/main")
```

### From Local Directory

```stata
net install diddesign, from("/path/to/diddesign-stata/")
```

### Get Example Files

```stata
net get diddesign
```

### Loading Example Data

The `diddesign_data` command provides a convenient way to load bundled example datasets:

```stata
* List all available datasets
diddesign_data, list

* Load a dataset by name
diddesign_data malesky2014, clear

* Describe a dataset without loading
diddesign_data anzia2012, describe
```

Available datasets: `anzia2012`, `malesky2014`, `paglayan2019`.

### Verify Installation

```stata
which diddesign
help diddesign
```

## Recommended Workflow

The Double DID follows a two-step workflow:

### Step 1: Assess Assumptions with `diddesign_check`

```stata
* Panel data
diddesign_check outcome, treatment(treat) id(unit_id) time(time_var) lag(1 2 3)

* Repeated cross-section data
diddesign_check outcome, treatment(treat) time(time_var) post(post_var) rcs lag(1 2)
```

**Output interpretation:**

| Component | What to Look For |
|-----------|------------------|
| **Placebo estimates** | Values close to zero support extended parallel trends |
| **Equivalence CI** | Narrower intervals = stronger evidence for parallel trends |

**Decision guide:**

| Scenario | Recommendation |
|----------|----------------|
| Placebo ≈ 0, narrow equivalence CI | Use **Double-DID** |
| Placebo ≈ 0, wide equivalence CI | Use **Double-DID** with caution |
| Significant placebo, trends same direction | Consider **s-DID** |
| Significant placebo, trends opposite | Neither assumption may hold |

Visualize trends with `diddesign_plot, type(both)` to complement statistical tests.

### Step 2: Estimate Treatment Effects with `diddesign`

```stata
* Panel data
diddesign outcome, treatment(treat) id(unit_id) time(time_var) nboot(200)

* Repeated cross-section data
diddesign outcome, treatment(treat) time(time_var) post(post_var) rcs nboot(200)
```

The command reports three estimators:

| Estimator | When to Use |
|-----------|-------------|
| **Double-DID** | Diagnostics support extended parallel trends |
| **DID** | Baseline comparison |
| **s-DID** | Diagnostics suggest linear divergence |

**Interpreting GMM weights**: `w_DID ≈ 1` favors standard DID; `w_sDID ≈ 1` favors sequential DID.

**Comparing estimates**: If Double-DID ≈ DID ≈ s-DID, results are robust across assumptions. If DID ≠ s-DID, consider potential bias in standard DID.

## Quick Start

### Example 1: Repeated Cross-Section Data (Malesky et al. 2014)

This example replicates Malesky, Nguyen, and Tran (2014) on the abolition of elected councils in Vietnam (Table 2 and Figure 4 in Egami and Yamauchi 2023). Treatment was assigned in 2009, with pre-treatment periods in 2006 and 2008. This is repeated cross-sectional (RCS) data where different communes are observed at each time period.

```stata
* Load example data (Vietnamese communes, treatment in 2009)
diddesign_data malesky2014, clear

* Drop observations with missing covariates (matches R replication code)
drop if missing(lnarea) | missing(lnpopden) | missing(city)

* Create region fixed effects dummies for covariates
tabulate reg8, gen(reg8_)

* Set random seed for reproducibility
set seed 1234

* Step 1: Assess parallel trends with covariates
* Covariates: area size, population density, city indicator, regional FE
* (matches the original study and R replication code)
diddesign_check pro4 lnarea lnpopden city ///
    reg8_2 reg8_3 reg8_4 reg8_5 reg8_6 reg8_7, ///
    treatment(treatment) time(year) ///
    post(post_treat) rcs ///
    cluster(id_district) ///
    lag(1) nboot(2000)

* Visualize diagnostics: trends plot + equivalence CI
diddesign_plot, type(both)

* Step 2: Estimate treatment effects with covariates
* level(90) matches the 90% CI reported in the original paper
diddesign pro4 lnarea lnpopden city ///
    reg8_2 reg8_3 reg8_4 reg8_5 reg8_6 reg8_7, ///
    treatment(treatment) time(year) ///
    post(post_treat) rcs ///
    cluster(id_district) ///
    nboot(2000) level(90)

* Plot treatment effect estimates with confidence intervals
diddesign_plot, type(estimates)
```

**Output — Parallel Trends Assessment (pro4):**

```
Parallel Trends Assessment              Number of obs         =      6,265
Design:        Standard DID             Number of clusters    =      1,080
Data type:     Repeated Cross-Section   Bootstrap reps        =      2,000
Standardized:  Yes                      Bootstrap valid       =      2,000
Equiv. CI:     95%

Clustering:    id_district

------------------------------------------------------------------------------
    Lag|  Estimate  Std. err.   Est.(raw)   SE(raw)       [95% eq. CI]
--------+---------------------------------------------------------------------
      1|   -0.0073    0.0975     -0.0032    0.0421   [-0.1683,  0.1683]
------------------------------------------------------------------------------
```

**Output — Treatment Effect Estimation (pro4):**

```
Double Difference-in-Differences        Number of obs         =      6,265
Design:        Standard DID
Data type:     Repeated Cross-Section   Number of periods     =          3
                                        Bootstrap reps        =      2,000

Clustering:    id_district

------------------------------------------------------------------------------
         pro4 |  Coef.   Std. err.     z   P>|z|   [90% conf. interval]  Weight
--------------+---------------------------------------------------------------
   Double DID |   0.0817    0.0516   1.58  0.113    -0.0032    0.1666        .
          DID |   0.0837    0.0575   1.46  0.145    -0.0109    0.1783   1.6120
         sDID |   0.0869    0.0845   1.03  0.304    -0.0520    0.2259  -0.6120
------------------------------------------------------------------------------

Note: Double DID combines DID and sDID via optimal GMM weights.
```

**Interpreting the results:**

The results replicate Table 2 and Figure 4 of Egami and Yamauchi (2023):

*Parallel trends assessment (Table 2):* The placebo estimate for "Education and Cultural Program" (pro4) is -0.007 (SE = 0.097, p = 0.940), indicating no significant deviation from parallel trends. The 95% standardized equivalence CI of [-0.168, 0.168] is narrow, providing positive evidence that pre-treatment trend differences are bounded within approximately 0.17 SD of the baseline control group mean. This supports the extended parallel trends (EPT) assumption.

*Treatment effect estimation (Figure 4):* The standard DID estimate of the ATT is 0.084 (90% CI [-0.011, 0.178]), while the Double-DID estimate is 0.082 (90% CI [-0.003, 0.167]). The Double-DID achieves a roughly 10% reduction in standard errors by optimally combining DID and sDID within the GMM framework. Since EPT is supported, the Double-DID is the preferred estimator. The positive point estimate suggests that abolishing elected councils modestly increased education and cultural programs, though the effect is marginally significant.

For the other two outcomes analyzed in the paper:
- *Tap water:* The placebo estimate is 0.166 (SE = 0.080, p = 0.038), suggesting EPT is questionable. The sDID estimate (-0.119) differs substantially from the standard DID (-0.078), indicating the standard DID underestimates the negative effect by about 50% due to non-parallel pre-treatment trends.
- *Agricultural center:* The placebo estimate is 0.198 (SE = 0.084, p = 0.019) with opposite-sign pre-treatment trends, meaning neither EPT nor parallel trends-in-trends is plausible. No credible causal estimate can be obtained without stronger assumptions.

### Example 2: Staggered Adoption Design (Paglayan 2019)

This example replicates Paglayan (2019) on collective bargaining rights for teachers (Appendix H.2 in Egami and Yamauchi 2023). States adopted policies at different times (staggered adoption). This is panel data where the same 49 states are observed from 1959 through 2000, with 14 unique treatment timings.

```stata
* Load example data (US states, staggered policy adoption 1959-2000)
diddesign_data paglayan2019, clear

* Note: WI and DC are already excluded in the Stata dataset,
* matching the R replication code: filter(!(state %in% c("WI", "DC")))

* Create log-transformed outcome variable
gen log_expenditure = log(pupil_expenditure + 1)

* Recompute log_salary to match R code: log(teacher_salary + 1)
drop log_salary
gen log_salary = log(teacher_salary + 1)

* Set random seed
set seed 1234

* Step 1: Assess parallel trends across multiple lags
* design(sa) activates staggered adoption mode
* thres(1) includes cohorts with at least 1 treated unit
diddesign_check log_expenditure, ///
    treatment(treatment) time(year) ///
    id(stateid) design(sa) ///
    lag(1 2 3 4 5) thres(1) nboot(2000)

* Visualize treatment timing pattern and placebo results
diddesign_plot, type(both)

* Step 2: Estimate dynamic treatment effects (lead 0-9)
* SA-Double-DID aggregates effects across treatment timing cohorts
diddesign log_expenditure, ///
    treatment(treatment) time(year) ///
    id(stateid) design(sa) ///
    thres(1) lead(0 1 2 3 4 5 6 7 8 9) nboot(2000)

* Plot treatment effect estimates
diddesign_plot, type(estimates)
```

**Output — Parallel Trends Assessment (log_expenditure):**

```
Parallel Trends Assessment              Number of obs         =      2,058
Design:        Staggered Adoption       Number of clusters    =         49
Data type:     Panel                    Bootstrap reps        =      2,000
Standardized:  Yes                      Bootstrap valid       =      2,000
Equiv. CI:     95%

------------------------------------------------------------------------------
    Lag|  Estimate  Std. err.   Est.(raw)   SE(raw)       [95% eq. CI]
--------+---------------------------------------------------------------------
      1|   -0.0197    0.0541     -0.0027    0.0087   [-0.1090,  0.1090]
      2|   -0.0691    0.0581     -0.0124    0.0091   [-0.1650,  0.1650]
      3|    0.0163    0.0639      0.0023    0.0111   [-0.1210,  0.1210]
      4|   -0.0447    0.0675     -0.0076    0.0120   [-0.1560,  0.1560]
      5|   -0.0488    0.0473     -0.0107    0.0094   [-0.1270,  0.1270]
------------------------------------------------------------------------------
```

**Output — Treatment Effect Estimation (log_expenditure, lead 0):**

```
Staggered Adoption Double DID           Number of obs         =      2,058
Design:        Staggered Adoption       Number of units       =         49
Data type:     Panel                    Number of periods     =         42
                                        Bootstrap reps        =      2,000

Clustering:    stateid

------------------------------------------------------------------------------
log_expendit~e |  Coef.   Std. err.     z   P>|z|   [95% conf. interval]  Weight
--------------+---------------------------------------------------------------
SA-Double DID |   0.0113    0.0127   0.89  0.373    -0.0130    0.0369        .
       SA-DID |   0.0110    0.0127   0.87  0.386    -0.0135    0.0363   0.8860
      SA-sDID |   0.0137    0.0148   0.93  0.354    -0.0159    0.0426   0.1140
------------------------------------------------------------------------------

Note: Double DID combines DID and sDID via optimal GMM weights.
```

**Interpreting the results:**

The results replicate Appendix H.2 (Figures H.3–H.5) of Egami and Yamauchi (2023):

*Parallel trends assessment (Figure H.4):* All five placebo estimates are close to zero, and the 95% standardized equivalence CIs are narrow for all lags (the widest is [-0.165, 0.165] at lag 2). This provides strong evidence supporting the extended parallel trends assumption for the SA design, consistent with the paper's finding.

*Treatment effect estimation (Figure H.5):* The SA-Double-DID estimates for per-pupil expenditure are small and not statistically significant for most lead periods (e.g., lead 0: 0.011, 95% CI [-0.013, 0.037]). The GMM weights favor SA-DID (w = 0.886 at lead 0), reflecting that the extended parallel trends assumption is well-supported. Dynamic effects across leads 0–9 show no systematic pattern of increasing or decreasing treatment effects, consistent with the original finding of Paglayan (2019) that granting collective bargaining rights did not significantly increase education resources.

For log teacher salary, the results are qualitatively similar: placebo tests support parallel trends (all equivalence CIs are narrow), and treatment effects are not statistically significant.

### Example 3: Panel Data (Standard DID)

For balanced panel data where the same units are observed across all time periods:

```stata
* Panel data structure:
* - unit_id: unit identifier (numeric or string; strings are auto-encoded)
* - time: time period identifier
* - treat: binary treatment indicator (0/1)
* - outcome: outcome variable

* Step 1: Check parallel trends
diddesign_check outcome, ///
    treatment(treat) id(unit_id) time(time) ///
    lag(1 2 3) nboot(200)

* Visualize trends and placebo results
diddesign_plot, type(both)

* Step 2: Estimate Double DID
diddesign outcome, ///
    treatment(treat) id(unit_id) time(time) ///
    nboot(200)

* With covariates (two equivalent approaches):
* Approach 1: list after outcome in varlist
diddesign outcome covariate1 covariate2, ///
    treatment(treat) id(unit_id) time(time) ///
    nboot(200)

* Approach 2: use covariates() option (supports factor variables)
diddesign outcome, ///
    treatment(treat) id(unit_id) time(time) ///
    covariates(i.region c.gdp) nboot(200)
```

**Note:** For publication-quality results, use `nboot(2000)`. The Staggered Adoption (SA) design only supports panel data, not RCS data.

## Commands

| Command             | Description                          |
| ------------------- | ------------------------------------ |
| `diddesign`       | Main estimation command              |
| `diddesign_check` | Parallel trends diagnostic           |
| `diddesign_plot`  | Visualization                        |
| `diddesign_data`  | Load bundled example datasets        |

## Options

### diddesign Options

| Option                    | Description                                   | Default            |
| ------------------------- | --------------------------------------------- | ------------------ |
| `treatment(varname)`    | Treatment indicator variable                  | Required           |
| `id(varname)`           | Unit identifier variable (for panel data)     | Required for panel |
| `time(varname)`         | Time identifier variable                      | Required           |
| `post(varname)`         | Post-treatment indicator (for RCS data)       | Required for RCS   |
| `nboot(#)`              | Number of bootstrap replications              | 30                 |
| `lead(numlist)`         | Lead periods for dynamic effects              | 0                  |
| `design(string)`        | Design type:`did` or `sa`                 | did                |
| `thres(#)`              | Minimum treated units per period (SA only)    | 2                  |
| `cluster(varname)`      | Cluster variable for bootstrap                | -                  |
| `covariates(fvvarlist)` | Control variables (supports factor variables) | -                  |
| `level(#)`              | Confidence level (%)                          | 95                 |
| `seed(#)`               | Random seed for reproducibility               | -                  |
| `seboot`                | Use bootstrap percentile confidence intervals | -                  |
| `panel`                 | Panel data format                             | -                  |
| `rcs`                   | Repeated cross-section data format            | -                  |

**Data Type Options:**

- For **panel data**: use `id(varname)` to specify unit identifiers
- For **RCS data**: use `post(varname)` with `rcs` option to specify post-treatment periods

**Factor Variable Support:**

The `covariates()` option supports Stata's factor variable notation:

```stata
* Categorical covariates with automatic dummy expansion
diddesign y, treatment(treat) id(id) time(time) covariates(i.region)

* Interaction terms
diddesign y, treatment(treat) id(id) time(time) covariates(c.x1#c.x2)

* Mixed factor and continuous
diddesign y, treatment(treat) id(id) time(time) covariates(i.region c.gdp)
```

### diddesign_check Options

| Option                 | Description                                | Default            |
| ---------------------- | ------------------------------------------ | ------------------ |
| `treatment(varname)` | Treatment indicator variable               | Required           |
| `time(varname)`      | Time identifier variable                   | Required           |
| `id(varname)`        | Unit identifier (for panel data)           | Required for panel |
| `post(varname)`      | Post-treatment indicator (for RCS data)    | Required for RCS   |
| `lag(numlist)`       | Lag periods for placebo test               | 1                  |
| `design(string)`     | Design type:`did` or `sa`              | did                |
| `thres(#)`           | Minimum treated units per period (SA only) | 2                  |
| `nboot(#)`           | Number of bootstrap replications           | 30                 |
| `cluster(varname)`   | Cluster variable for bootstrap             | -                  |
| `seed(#)`            | Random seed for reproducibility            | -                  |
| `panel`              | Panel data format                          | -                  |
| `rcs`                | Repeated cross-section data format         | -                  |

### diddesign_plot Options

| Option               | Description                                                             | Default        |
| -------------------- | ----------------------------------------------------------------------- | -------------- |
| `type(string)`     | Plot type:`trends`, `placebo`, `pattern`, `both`, `estimates` | estimates/both |
| `saving(filename)` | Save graph to file                                                      | -              |
| `replace`          | Overwrite existing file                                                 | -              |
| `scheme(string)`   | Graph scheme name                                                       | -              |
| `title(string)`    | Graph title                                                             | -              |
| `xtitle(string)`   | X-axis title                                                            | -              |
| `ytitle(string)`   | Y-axis title                                                            | -              |
| `ci`               | Display confidence bands on trends plot                                 | -              |
| `band`             | Display CI as ribbon on estimates plot                                  | -              |
| `name(string)`     | Graph name in memory                                                    | -              |

**Plot Types:**

- `estimates`: Double-DID estimates across lead values with 90% CI (default after `diddesign`)
- `trends`: Outcome trajectories for treated and control groups (standard DID only)
- `placebo`: Pre-treatment placebo test with equivalence CI
- `pattern`: Treatment timing heatmap (SA design only)
- `both`: Combined diagnostic plot (default after `diddesign_check`)

## Stored Results

### diddesign

**Scalars:**

| Result           | Description                       |
| ---------------- | --------------------------------- |
| `e(N)`         | Number of observations            |
| `e(n_units)`   | Number of units (panel data only) |
| `e(n_periods)` | Number of time periods            |
| `e(n_boot)`    | Number of bootstrap replications  |
| `e(n_lead)`    | Number of lead periods            |
| `e(level)`     | Confidence level                  |
| `e(is_panel)`  | Panel data indicator              |
| `e(seboot)`    | Bootstrap CI indicator            |

**Macros:**

| Result            | Description                               |
| ----------------- | ----------------------------------------- |
| `e(cmd)`        | Command name:`diddesign`                |
| `e(cmdline)`    | Full command line                         |
| `e(design)`     | Design type:`did` or `sa`             |
| `e(depvar)`     | Dependent variable name                   |
| `e(treatment)`  | Treatment variable name                   |
| `e(covariates)` | Covariate list                            |
| `e(id)`         | Unit identifier variable                  |
| `e(time)`       | Time variable                             |
| `e(clustvar)`   | Cluster variable                          |
| `e(lead)`       | Lead values                               |
| `e(ci_method)`  | CI method:`asymptotic` or `bootstrap` |

**Matrices:**

| Result           | Description                                   |
| ---------------- | --------------------------------------------- |
| `e(b)`         | Coefficient vector [Double-DID, DID, sDID]    |
| `e(V)`         | Variance-covariance matrix                    |
| `e(estimates)` | Results table: lead, estimate, SE, CI, weight |
| `e(weights)`   | GMM weights matrix: w_did, w_sdid             |
| `e(W)`         | GMM optimal weight matrix (2×2)              |
| `e(vcov_gmm)`  | GMM variance-covariance matrix (2×2)         |

### diddesign_check

**Scalars:**

| Result              | Description                          |
| ------------------- | ------------------------------------ |
| `e(N)`            | Number of observations               |
| `e(n_lags)`       | Number of lag periods                |
| `e(n_boot)`       | Number of bootstrap replications     |
| `e(n_boot_valid)` | Number of valid bootstrap iterations |
| `e(n_clusters)`   | Number of clusters                   |
| `e(is_panel)`     | Panel data indicator                 |

**Macros:**

| Result           | Description                      |
| ---------------- | -------------------------------- |
| `e(cmd)`       | Command name:`diddesign_check` |
| `e(cmdline)`   | Full command line                |
| `e(design)`    | Design type:`did` or `sa`    |
| `e(depvar)`    | Dependent variable name          |
| `e(treatment)` | Treatment variable name          |
| `e(clustvar)`  | Cluster variable                 |
| `e(datatype)`  | Data type:`panel` or `rcs`   |

**Matrices:**

| Result         | Description                                   |
| -------------- | --------------------------------------------- |
| `e(placebo)` | Placebo test results: lag, estimate, SE, EqCI |
| `e(trends)`  | Trends data: time, group, mean, SD, N         |
| `e(Gmat)`    | Treatment timing matrix (SA only)             |

## Methodology

### GMM Framework

The Double DID is formulated as a GMM problem:

```
τ̂_d-DID = argmin  [τ - τ̂_DID  ]ᵀ     [τ - τ̂_DID  ]
             τ     [τ - τ̂_sDID ]  W   [τ - τ̂_sDID ]
```

where **W** is a 2×2 weight matrix. Setting **W** = Σ⁻¹ (inverse of variance-covariance) yields the optimal Double DID. Other choices of **W** recover standard DID, sequential DID, or extended DID as special cases. In panel data, the extended DID with equal weights is numerically equivalent to the TWFE estimator.

### Optimal Weights

The optimal weights minimize variance:

```
w₁ = [Var(τ̂_sDID) - Cov(τ̂_DID, τ̂_sDID)] / [Var(τ̂_DID) + Var(τ̂_sDID) - 2·Cov(τ̂_DID, τ̂_sDID)]
w₂ = 1 - w₁
```

The Double DID estimator is: `τ̂_d-DID = w₁ · τ̂_DID + w₂ · τ̂_sDID`

### Theoretical Properties

Under the extended parallel trends assumption:

```
Var(τ̂_d-DID) ≤ min{ Var(τ̂_DID), Var(τ̂_sDID), Var(τ̂_e-DID) }
```

That is, the Double DID achieves variance no larger than standard DID, sequential DID, or extended DID when the extended parallel trends assumption holds. When only the parallel trends-in-trends assumption holds, the Double DID reduces to the sequential DID estimator.

### Staggered Adoption (SA) Extension

The SA design handles settings where units receive treatment at different times.

| Feature | Description |
|---------|-------------|
| **Control group** | Only not-yet-treated units serve as controls |
| **Cohort-specific effects** | Estimates SA-ATT(t) for each treatment timing t |
| **Aggregation** | Time-average SA-ATT = Σ π_t × SA-ATT(t), where π_t is cohort proportion |

The SA extension applies the Double DID to each treatment cohort separately, then aggregates.

## References

Egami, N., & Yamauchi, S. (2023). Using Multiple Pretreatment Periods to Improve Difference-in-Differences and Staggered Adoption Designs. *Political Analysis*, 31(2), 195-212. https://doi.org/10.1017/pan.2022.8

## Authors

**Stata Implementation:**

- **Xuanyu Cai**, City University of Macau
  Email: [xuanyuCAI@outlook.com](mailto:xuanyuCAI@outlook.com)
- **Wenli Xu**, City University of Macau
  Email: [wlxu@cityu.edu.mo](mailto:wlxu@cityu.edu.mo)

**Methodology:**

- **Naoki Egami**, Columbia University
- **Soichiro Yamauchi**, Harvard University

## License

AGPL-3.0. See [LICENSE](LICENSE) for details.

## Citation

If you use this package in your research, please cite both the methodology paper and the Stata implementation:

**APA Format:**

> Cai, X., & Xu, W. (2025). *diddesign: Stata module for Double Difference-in-Differences estimation* (Version 0.1.0) [Computer software]. GitHub. https://github.com/gorgeousfish/diddesign
>
> Egami, N., & Yamauchi, S. (2023). Using Multiple Pretreatment Periods to Improve Difference-in-Differences and Staggered Adoption Designs. *Political Analysis*, 31(2), 195-212. https://doi.org/10.1017/pan.2022.8

**BibTeX:**

```bibtex
@software{diddesign2025stata,
  title={diddesign: Stata module for Double Difference-in-Differences estimation},
  author={Xuanyu Cai and Wenli Xu},
  year={2025},
  version={0.1.0},
  url={https://github.com/gorgeousfish/diddesign}
}

@article{egami2023using,
  title={Using Multiple Pretreatment Periods to Improve Difference-in-Differences and Staggered Adoption Designs},
  author={Egami, Naoki and Yamauchi, Soichiro},
  journal={Political Analysis},
  volume={31},
  number={2},
  pages={195--212},
  year={2023},
  doi={10.1017/pan.2022.8}
}
```

## See Also

- Original R package by Egami and Yamauchi: https://github.com/naoki-egami/DIDdesign
- Paper: Egami, N., & Yamauchi, S. (2023). https://doi.org/10.1017/pan.2022.8
