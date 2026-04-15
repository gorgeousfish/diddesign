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

The package now supports the **generalized K-DID** extension from Appendix E
of the paper. When more than two pre-treatment periods are available, specify
`kmax(K)` to combine up to K component estimators via GMM, allowing for
(K-1)-th degree polynomial time-varying confounding. The default `kmax(2)`
preserves backward compatibility with the basic 2-moment Double DID.

When multiple pre-treatment periods are available, researchers can exploit three uses:

1. **Assessing the parallel trends assumption**: Pre-treatment periods enable placebo tests by applying the DID estimator to pre-treatment outcomes. The package also reports equivalence confidence intervals, which provide positive evidence for approximate parallel trends when the interval is narrow.

2. **Improving estimation accuracy**: Under the extended parallel trends assumption, multiple pre-treatment periods provide additional information to estimate counterfactual outcomes, reducing variance.

3. **Allowing for a more flexible parallel trends assumption**: The sequential DID estimator requires only the parallel trends-in-trends assumption, which permits linear time-varying confounding. In the paper's theoretical moment-selection result, when only the parallel trends-in-trends assumption is plausible, the identifying moment reduces to the sequential DID estimator. The current `diddesign` command still reports the GMM-weighted `Double-DID` row together with `DID` and `s-DID`, so users should interpret those estimates alongside `diddesign_check` diagnostics rather than expecting the displayed `Double-DID` row to automatically collapse to `s-DID`.

The Double DID combines these uses within a unified GMM framework. In panel
data with multiple pre-treatment periods, the two-way fixed effects (TWFE)
estimator is numerically equivalent to a weighted average of DIDs with equal
weights.

**Features**:

- Double DID estimation via GMM (K=2 default)
- **Generalized K-DID** for K>2 pre-treatment periods (`kmax()` option)
- **J-test moment selection** for adaptive moment choice (`jtest()` option)
- Staggered adoption (SA) design support (including K-DID extension)
- **Parallel bootstrap** for faster computation (`parallel` option, requires SSC `parallel` package)
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
| Placebo close to zero, narrow equivalence CI | Supports the extended parallel trends interpretation, with final judgment still tied to substantive domain knowledge |
| Placebo non-zero but trends have same direction | May motivate a parallel-trends-in-trends interpretation only with substantive justification |
| Trends have opposite directions | May indicate neither assumption is credible without stronger assumptions |

### Equivalence Confidence Intervals

Traditional hypothesis testing interprets non-rejection as evidence for parallel trends, but this conflates "no evidence against" with "evidence for." The equivalence approach reverses the null hypothesis: rejection of the null (that trends differ by more than some threshold) provides positive evidence for approximate parallel trends.

The 95% standardized equivalence CI is reported in units of baseline control group standard deviations. It summarizes the smallest symmetric equivalence range supported by the observed data at the 5% significance level. As in the paper and the reference R package, there is no universal EqCI cutoff: researchers should use substantive domain knowledge to decide whether the reported interval is narrow enough for their application.

## Requirements

- Stata 16.0 or later
- No additional dependencies (Mata included in Stata)

## Installation

### Install the package

```stata
net install diddesign, from("https://raw.githubusercontent.com/gorgeousfish/diddesign/main")
```

This automatically installs:
- All commands and help files
- Example data files (`malesky2014.dta`, `paglayan2019.dta`)

### Verify Installation

```stata
which diddesign
help diddesign
sysuse dir
```

## Quick Start with Examples

After installation, the example datasets are automatically available:

```stata
* Load Vietnam communes data
sysuse malesky2014, clear

* Load US states panel data  
sysuse paglayan2019, clear

* Download optional example scripts to the current directory
net get diddesign

* Run complete example analyses
do example_malesky.do
do example_paglayan.do
```

**Note:** Bundled datasets are installed to Stata's example-data search path, so use `sysuse`, not `use`.

## Recommended Workflow

The Double DID follows a two-step workflow:

### Step 1: Assess Assumptions with `diddesign_check`

```stata
diddesign_check outcome, treatment(treatment_var) id(unit_id) time(time_var) lag(1 2 3)
```

**Output interpretation:**

| Component | What to Look For |
|-----------|------------------|
| **Placebo estimates** | Values close to zero support extended parallel trends |
| **Std. Err.** | Smaller standard errors make placebo deviations easier to interpret |
| **Equivalence CI** | Narrower intervals = stronger evidence that pre-trends are substantively small |

**Decision guide:**

| Scenario | Recommendation |
|----------|----------------|
| Placebo ≈ 0, narrow equivalence CI | Supports the extended parallel trends interpretation, subject to substantive domain knowledge |
| Placebo ≈ 0, wide equivalence CI | Evidence for extended parallel trends is weak, so avoid a mechanical estimator choice |
| Placebo clearly away from zero, trends same direction | A parallel-trends-in-trends reading may be plausible only with substantive justification |
| Placebo clearly away from zero, trends opposite | Neither assumption is credible without stronger assumptions |

Visualize trends with `diddesign_plot, type(both)` to complement statistical tests.

### Step 2: Estimate Treatment Effects with `diddesign`

```stata
diddesign outcome, treatment(treatment_var) id(unit_id) time(time_var) nboot(200)
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

### Example 1: Repeated Cross-Section (Malesky et al. 2014)

This example uses the bundled `malesky2014.dta` data, which is a repeated cross-section (RCS) design with pre-treatment periods in 2006 and 2008 and a post-treatment period in 2010.

```stata
* Load bundled example data
sysuse malesky2014, clear

* Set random seed for reproducibility
set seed 1234

* Step 1: Assess parallel trends using pre-treatment periods
* lag(1) compares 2006-2008 as a placebo DID
* Cluster at the treatment-assignment level for RCS inference
diddesign_check pro4, ///
    treatment(treatment) ///
    time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    lag(1)

* Visualize diagnostics: trends plot + equivalence CI
diddesign_plot, type(both)

* Step 2: Estimate treatment effects
* Output shows Double-DID, DID, and sDID estimates
diddesign pro4, ///
    treatment(treatment) ///
    time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

* Plot treatment effect estimates with confidence intervals
diddesign_plot, type(estimates)
```

### Example 2: Staggered Adoption Design

This example replicates Paglayan (2019) on collective bargaining rights for teachers. States adopted policies at different times (staggered adoption).

```stata
* Load bundled example data
sysuse paglayan2019, clear
gen log_expenditure = log(pupil_expenditure + 1)
encode state, gen(id_subject)

* Set random seed
set seed 1234

* Step 1: Assess parallel trends across multiple lags
* thres(1) requires at least 1 treated unit per cohort
diddesign_check log_expenditure, ///
    treatment(treatment) ///
    id(id_subject) time(year) design(sa) lag(1 2 3 4 5) thres(1)

* Visualize treatment timing pattern (which states treated when)
diddesign_plot, type(pattern)

* Step 2: Estimate time-weighted average treatment effect
* SA-Double-DID aggregates effects across treatment timing cohorts
diddesign log_expenditure, ///
    treatment(treatment) ///
    id(id_subject) time(year) design(sa) thres(1) nboot(200)
```

### Example 3: Repeated Cross-Section (RCS) Data

For data where different units are observed at each time period (no panel tracking):

```stata
* RCS data structure:
* - treat_group: binary (1 = treatment group, 0 = control group)
* - post: binary (1 = post-treatment period, 0 = pre-treatment)
*   values within 1e-6 of 0 or 1 are canonicalized to exact binary values
* - time: time period identifier
* - cluster_id: treatment-assignment cluster for bootstrap inference

* Step 1: Check parallel trends
diddesign_check outcome, treatment(treat_group) time(time) post(post) rcs ///
    cluster(cluster_id) lag(1 2)

* Visualize
diddesign_plot

* Step 2: Estimate Double DID
diddesign outcome, treatment(treat_group) time(time) post(post) rcs ///
    cluster(cluster_id) nboot(200)
```

**Note:** The Staggered Adoption (SA) design only supports panel data, not RCS data. For RCS inference, `cluster()` is required and should point to the treatment-assignment bootstrap block; `diddesign` and `diddesign_check` now reject RCS runs that omit `cluster()`.

### Example 4: Generalized K-DID with Multiple Pre-treatment Periods

When 3+ pre-treatment periods are available, the generalized K-DID can allow for higher-order polynomial confounding:

```stata
* Panel with 5+ time periods (≥3 pre-treatment)
* K=3 allows for quadratic time-varying confounding
diddesign outcome, treatment(treat) id(unit_id) time(time) ///
    nboot(200) kmax(3)

* With J-test moment selection (adaptive)
diddesign outcome, treatment(treat) id(unit_id) time(time) ///
    nboot(200) kmax(3) jtest(on)

* SA design with K=3
diddesign outcome, treatment(treat) id(unit_id) time(time) ///
    design(sa) nboot(200) kmax(3)
```

**Interpreting K-DID output:**
- **K-DID (K=n)**: GMM-optimal combination of n component estimators
- **k=1 (DID)**: Standard parallel trends assumption
- **k=2 (sDID)**: Parallel trends-in-trends (linear confounding)
- **k=3**: Quadratic time-varying confounding
- **Weights**: GMM-optimal weights for each component

## Commands

| Command             | Description                |
| ------------------- | -------------------------- |
| `diddesign`       | Main estimation command    |
| `diddesign_check` | Parallel trends diagnostic |
| `diddesign_plot`  | Visualization              |

## Options

### diddesign Options

| Option                    | Description                                   | Default            |
| ------------------------- | --------------------------------------------- | ------------------ |
| `treatment(varname)`    | Treatment indicator variable                  | Required           |
| `id(varname)`           | Unit identifier variable (for panel data)     | Required for panel |
| `time(varname)`         | Time identifier variable                      | Required           |
| `post(varname)`         | Post-treatment indicator for RCS data; values within 1e-6 of 0/1 are canonicalized before validation and `post()` triggers RCS auto-detection | Required for RCS   |
| `nboot(#)`              | Number of bootstrap replications (`# >= 2`)  | 30                 |
| `lead(numlist)`         | Lead periods for dynamic effects              | 0                  |
| `design(string)`        | Design type:`did` or `sa`                 | did                |
| `thres(#)`              | Minimum treated units per period (SA only)    | 2                  |
| `cluster(varname)`      | Cluster variable for bootstrap; required for RCS and should identify the treatment-assignment block | Optional for panel |
| `covariates(fvvarlist)` | Control variables (supports factor variables) | -                  |
| `level(#)`              | Confidence level (%)                          | 95                 |
| `seed(#)`               | Random seed for reproducibility               | -                  |
| `kmax(#)`               | Max K-DID components; requires ≥K pre-treatment periods | 2              |
| `jtest(on\|off)`        | J-test moment selection (Hansen overidentification test) | off            |
| `seboot`                | Use bootstrap percentile confidence intervals | -                  |
| `parallel`              | Enable parallel bootstrap (requires SSC `parallel` package) | -               |
| `panel`                 | Panel data format                             | -                  |
| `rcs`                   | Repeated cross-section data format; optional when `post()` already implies RCS | -                  |

**Data Type Options:**

- For **panel data**: use `id(varname)` to specify unit identifiers
- For **RCS data**: use `post(varname)` to specify post-treatment periods; values within 1e-6 of 0/1 are canonicalized to exact binary values before validation, and `rcs` is optional and only makes the data type explicit

When `parallel` is specified, bootstrap replications are distributed across multiple Stata instances using the SSC `parallel` package (`ssc install parallel`). If the package is not available, computation falls back to sequential bootstrap with a warning.

**SA Stored Results Note:**

When `design(sa)` is used with multiple requested `lead()` values, `e(time_weights)` summarizes the union of adoption periods that enter the posted common-support SA targets, while `e(time_weights_by_lead)` contains the lead-specific weight surface with one column per requested lead and column-wise re-normalization only over periods where both SA-DID and SA-sDID are jointly identified for that lead.

**Factor Variable Support:**

The `covariates()` option supports Stata's factor variable notation. String categorical variables are accepted directly inside factor notation and will be encoded automatically before expansion:

```stata
* Categorical covariates with automatic dummy expansion
diddesign y, treatment(treat) id(id) time(time) covariates(i.region)

* Interaction terms
diddesign y, treatment(treat) id(id) time(time) covariates(c.x1#c.x2)

* Mixed factor and continuous
diddesign y, treatment(treat) id(id) time(time) covariates(i.region c.gdp)
```

**Time Variable Requirement:**

`time()` may be numeric or string. String time labels are automatically encoded to numeric period indices only when that encoding preserves the observed order of distinct labels in the estimation sample. If labels such as `t1`, `t2`, and `t10` would be silently reordered lexicographically or by their numeric suffixes during encoding, `diddesign` stops with an error and asks users to recode `time()` to a numeric or lexically ordered variable first.

### diddesign_check Options

| Option                 | Description                                | Default            |
| ---------------------- | ------------------------------------------ | ------------------ |
| `treatment(varname)` | Treatment indicator variable               | Required           |
| `time(varname)`      | Time identifier variable                   | Required           |
| `id(varname)`        | Unit identifier (for panel data)           | Required for panel |
| `post(varname)`      | Post-treatment indicator for RCS data; values within 1e-6 of 0/1 are canonicalized before validation and `post()` triggers RCS auto-detection | Required for RCS   |
| `lag(numlist)`       | Lag periods for placebo test               | 1                  |
| `design(string)`     | Design type:`did` or `sa`              | did                |
| `thres(#)`           | Minimum treated units per period (SA only) | 2                  |
| `nboot(#)`           | Number of bootstrap replications (`# >= 2`) | 30                 |
| `cluster(varname)`   | Cluster variable for bootstrap; required for RCS and should identify the treatment-assignment block | Optional for panel |
| `seed(#)`            | Random seed for reproducibility            | -                  |
| `quiet`              | Suppress bootstrap progress display | - |
| `parallel`           | Enable parallel bootstrap (requires SSC `parallel` package) | -           |
| `panel`              | Panel data format                          | -                  |
| `rcs`                | Repeated cross-section data format; optional when `post()` already implies RCS | -                  |

`lag(0)` is allowed in `diddesign_check`. It compares the treatment period to the immediately preceding period and therefore is reported with a warning instead of being treated as a true placebo test.

For the main `diddesign` command, `quiet` suppresses the bootstrap progress display only; non-fatal warnings may still appear when the command needs to explain filtered leads or bootstrap support issues.

When `quiet` is specified, `diddesign_check` suppresses the bootstrap progress display and the display-stage placebo diagnostic notes that explicitly follow the `quiet` contract, including the interpretive `lag(0)` note and warnings about dropped or raw-only lags.

`quiet` does not suppress input-validation or preprocessing notices that occur before the diagnostic display stage, such as duplicate covariates being removed.

`thres()` is only valid for `design(sa)`. Standard DID workflows reject it rather than silently ignoring it.

`parallel` enables parallel bootstrap in `diddesign_check` using the SSC `parallel` package; falls back to sequential if unavailable.

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
| `ci`               | Display 90% SD-based bands around group-period mean outcomes on trends plot | -              |
| `band`             | Display CI as ribbon on estimates plot                                  | -              |
| `name(string)`     | Graph name in memory                                                    | -              |
| `use_check(name)`  | Overlay stored `diddesign_check` placebo results on an estimates plot   | -              |

**Plot Types:**

- `estimates`: Double-DID estimates across lead values with 90% CI (default after `diddesign`)
- `trends`: Outcome trajectories for treated and control groups (standard DID only)
- `placebo`: Pre-treatment placebo test with equivalence CI
- `pattern`: Treatment timing heatmap (SA design only)
- `both`: Combined diagnostic plot (default after `diddesign_check`)

When `use_check(name)` is used, the stored `diddesign_check` result must match the
current `diddesign` result on design, dependent variable, treatment variable,
data type, cluster variable, covariates, and `if/in` sample. For panel results,
the two objects must also use the same `id()` and `time()` definitions. For RCS results,
the two objects must also use the same `post()` definition.

## Stored Results

### diddesign

**Scalars:**

| Result           | Description                       |
| ---------------- | --------------------------------- |
| `e(N)`         | Number of observations            |
| `e(n_units)`   | Number of units (panel data only) |
| `e(n_periods)` | Number of time periods            |
| `e(n_boot)`    | Number of bootstrap replications  |
| `e(n_boot_success)` | Number of bootstrap iterations with at least one finite DID/sDID component estimate |
| `e(n_lead)`    | Number of retained lead periods with at least one posted coefficient |
| `e(n_lead_requested)` | Number of leads requested in `lead()` |
| `e(n_lead_filtered)` | Number of requested leads filtered out before estimation |
| `e(n_lead_identified)` | Number of retained leads with at least one posted coefficient in `e(b)` |
| `e(n_clusters)` | Number of clusters in the posted sample   |
| `e(level)`     | Confidence level                  |
| `e(is_panel)`  | Panel data indicator              |
| `e(seboot)`    | Bootstrap CI indicator            |
| `e(kmax)`      | Maximum K-DID components requested |
| `e(jtest_on)`  | J-test moment selection enabled (1/0) |
| `e(parallel)`  | Parallel bootstrap was used (1/0)  |
| `e(n_periods_valid)` | Number of adoption periods that actually enter SA estimation |
| `e(thres)`     | Minimum treated units threshold used |

**Macros:**

| Result            | Description                               |
| ----------------- | ----------------------------------------- |
| `e(cmd)`        | Command name:`diddesign`                |
| `e(cmdline)`    | Full command line                         |
| `e(design)`     | Design type:`did` or `sa`             |
| `e(depvar)`     | Dependent variable name                   |
| `e(treatment)`  | Treatment variable name                   |
| `e(covariates)` | Covariate list                            |
| `e(covars)`    | Alias of `e(covariates)` for cross-command compatibility |
| `e(id)`         | Unit identifier variable                  |
| `e(time)`       | Time variable                             |
| `e(post)`       | Post-treatment indicator (RCS only)       |
| `e(datatype)`   | Data type: `panel` or `rcs`               |
| `e(clustvar)`   | Cluster variable                          |
| `e(sample_ifin)` | Stored `if/in` sample restriction; used to verify `use_check()` sample compatibility |
| `e(lead)`       | Lead values with at least one coefficient posted in `e(b)` and `e(V)` |
| `e(requested_lead)` | Lead values originally requested in `lead()` |
| `e(filtered_lead)` | Requested lead values filtered out before estimation |
| `e(identified_lead)` | Retained lead values with at least one coefficient posted in `e(b)` and `e(V)` |
| `e(unidentified_lead)` | Requested lead values with no coefficient posted in `e(b)` and `e(V)` |
| `e(time_weight_labels)` | Pipe-delimited adoption period labels for `e(time_weights)` (SA only) |
| `e(properties)` | `b V` |
| `e(ci_method)`  | CI method:`asymptotic` or `bootstrap` |

**Matrices:**

| Result           | Description                                   |
| ---------------- | --------------------------------------------- |
| `e(b)`         | Coefficient vector [Double-DID, DID, sDID]    |
| `e(V)`         | Variance-covariance matrix of the posted estimator vector. In mixed-support multi-lead SA runs, this matrix is rebuilt from jointly observed posted bootstrap draws for postestimation. |
| `e(estimates)` | Results table: lead, estimate, std_error, ci_lo, ci_hi, weight. In mixed-support multi-lead SA runs, component DID/sDID `std_error` and `ci_*` use lead-specific marginal bootstrap draws, so they can differ from `diag(e(V))`. |
| `e(lead_values)` | posted lead values vector aligned with `e(b)` and `e(V)` |
| `e(bootstrap_support)` | Per-lead bootstrap support counts with columns `n_valid_did`, `n_valid_sdid`, and `n_joint_valid`; rows align with `e(weights)` and retained `lead()` values |
| `e(weights)`   | GMM weights matrix: w_did, w_sdid             |
| `e(W)`         | GMM optimal weight matrix; single lead is `2×2`, multiple leads are one row per lead with a flattened `2×2` matrix |
| `e(vcov_gmm)`  | GMM variance-covariance matrix; single lead is `2×2`, multiple leads are one row per lead with a flattened `2×2` matrix |

**Matrices (K-DID with `kmax()>2` only):**

| Result           | Description                                   |
| ---------------- | --------------------------------------------- |
| `e(k_summary)` | K summary: K_init, K_sel, K_final per lead    |
| `e(moment_selected)` | Moment selection mask (0/1) per lead/component |
| `e(moment_dropped_jtest)` | Moments dropped by J-test (0/1)       |
| `e(moment_dropped_numerical)` | Moments dropped for numerical reasons (0/1) |
| `e(jtest_stats)` | J-test statistics: J_stat, J_df, J_pval per lead |

**Matrices (SA design only):**

| Result           | Description                                   |
| ---------------- | --------------------------------------------- |
| `e(time_weights)` | Union-level SA time weights over adoption periods that enter the posted common-support SA targets under the requested `lead()` |
| `e(time_weight_periods)` | Adoption periods corresponding to each row of `e(time_weights)` |
| `e(time_weights_by_lead)` | Lead-specific SA time weights with rows aligned to `e(time_weight_periods)` and one re-normalized column per requested lead over periods where both SA-DID and SA-sDID are jointly identified for that lead |

**Functions:**

| Result         | Description            |
| -------------- | ---------------------- |
| `e(sample)`  | Marks the posted estimation sample, that is, the union of retained lead-specific support rows that enter at least one DID/sDID component after listwise deletion |

### diddesign_check

**Scalars:**

| Result              | Description                          |
| ------------------- | ------------------------------------ |
| `e(N)`            | Number of observations in the posted diagnostic sample represented by `e(placebo)` and `e(sample)` |
| `e(n_lags)`       | Number of identified lag periods     |
| `e(n_lags_posted)` | Number of lag values posted in `e(b)` and `e(V)` |
| `e(n_boot)`       | Number of bootstrap replications     |
| `e(n_boot_valid)` | Number of valid bootstrap iterations |
| `e(n_clusters)`   | Number of clusters in the posted union diagnostic sample |
| `e(level)`        | Confidence level used for placebo confidence intervals (90) |
| `e(is_panel)`     | Panel data indicator                 |

**Macros:**

| Result           | Description                      |
| ---------------- | -------------------------------- |
| `e(cmd)`       | Command name:`diddesign_check` |
| `e(cmdline)`   | Full command line                |
| `e(design)`    | Design type:`did` or `sa`    |
| `e(depvar)`    | Dependent variable name          |
| `e(treatment)` | Treatment variable name          |
| `e(id)`        | Unit identifier variable for panel results |
| `e(time)`      | Time variable used to index diagnostic periods |
| `e(post)`      | Post-treatment indicator for RCS results |
| `e(clustvar)`  | Cluster variable                 |
| `e(covariates)` | Covariate variables in original model syntax |
| `e(covars)`    | Alias of `e(covariates)` for cross-command compatibility |
| `e(identified_lags)` | Lag values retained in `e(placebo)` after support checks |
| `e(posted_lags)` | Lag values with standardized placebo estimates posted in `e(b)` and `e(V)` |
| `e(raw_only_lags)` | Lag values retained with raw placebo results only |
| `e(unidentified_lags)` | Requested lag values dropped because placebo inference was not identifiable |
| `e(datatype)`  | Data type:`panel` or `rcs`   |
| `e(sample_ifin)` | Stored `if/in` sample restriction for the diagnostic object; used to match compatible estimation samples |
| `e(properties)` | `b V` when standardized placebo results are posted; otherwise empty |

**Matrices:**

| Result         | Description                                   |
| -------------- | --------------------------------------------- |
| `e(b)`       | Posted standardized placebo coefficient vector when at least one lag is posted |
| `e(V)`       | Variance-covariance matrix of posted standardized placebo estimates |
| `e(placebo)` | Placebo test results: `lag`, `estimate`, `std_error`, `estimate_orig`, `std_error_orig`, `EqCI95_LB`, `EqCI95_UB` |
| `e(trends)`  | Trends data: time, group, mean, SD, N         |
| `e(n_boot_valid_lag)` | Lag-specific bootstrap support counts for standardized and raw placebo inference |
| `e(n_clusters_lag)` | Lag-specific cluster support matrix with one `n_clusters` column (standard DID only) |
| `e(Gmat)`    | Treatment timing matrix (SA only)             |

**Functions:**

| Result         | Description                                   |
| -------------- | --------------------------------------------- |
| `e(sample)`  | Posted diagnostic sample indicator for the union of retained placebo support rows, not the full raw command sample |

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

That is, the Double DID achieves variance no larger than standard DID, sequential DID, or extended DID when the extended parallel trends assumption holds. In the paper's theoretical moment-selection result, when only the parallel trends-in-trends assumption holds, the Double DID reduces to the sequential DID estimator. The current `diddesign` command still reports the GMM-weighted `Double-DID` row together with `DID` and `s-DID`; users should interpret those estimates alongside `diddesign_check` diagnostics rather than expecting the reported `Double-DID` row to automatically collapse to `s-DID`.

When `lead()>0` is requested, `diddesign` follows the dynamic extension in Appendix E.2.3 of the paper rather than reusing the static contemporaneous contrast. In that case, the DID and sequential DID components are evaluated for the requested post-treatment lead and remain anchored at the last pre-treatment period (`-1` on the standardized event-time scale), matching the lead-aware implementation used in the reference R code.

### Staggered Adoption (SA) Extension

The SA design handles settings where units receive treatment at different times.

The current `design(sa)` implementation requires complete and unique `id()` x `time()` cells, that is, a balanced panel with exactly one observation per unit-period. It also requires at least three distinct time periods because the SA DID / sDID estimators use the `{t-2, t-1, t}` window.

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
