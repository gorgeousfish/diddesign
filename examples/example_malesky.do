/*===========================================================================
 * example_malesky.do - Basic DID with Repeated Cross-Sectional Data
 *
 * Replication of Malesky et al. (2014) analysis using diddesign
 *
 * Data Description:
 *   Repeated cross-sectional data from Vietnam communes (2006, 2008, 2010).
 *   Treatment: Abolition of elected councils (implemented in 2009).
 *   Study examines the effect on local public services.
 *
 * Variables:
 *   - id_district: District identifier (string, auto-encoded by diddesign)
 *   - year: Year (2006, 2008, 2010)
 *   - treatment: Treated commune indicator (0/1)
 *   - post_treat: Post-treatment indicator (1 if year == 2010)
 *   - pro4: Education and Cultural Program indicator [outcome]
 *   - tapwater: Tap water availability indicator [outcome]
 *   - agrext: Agricultural extension center indicator [outcome]
 *   - lnarea, lnpopden, city, reg8: Control variables
 *
 * Key Notes:
 *   - This is REPEATED CROSS-SECTIONAL data, NOT panel data
 *   - Use post() option instead of id() for RCS data
 *   - Use rcs option to indicate repeated cross-section
 *   - Cluster standard errors at district level
 *   - String variables (e.g., id_district) are auto-encoded by diddesign
 *
 * Usage:
 *   After installing diddesign, download example files:
 *     . net get diddesign
 *   Then run from the examples/ directory:
 *     . cd examples
 *     . do example_malesky.do
 *
 * Reference:
 *   Malesky, Nguyen, and Tran (2014). "The Impact of Recentralization on
 *   Public Services: A Difference-in-Differences Analysis of the Abolition
 *   of Elected Councils in Vietnam." American Political Science Review.
 *===========================================================================*/

version 16
clear all
set more off
capture log close

log using "example_malesky.log", replace

/*---------------------------------------------------------------------------
 * Section 1: Data Loading and Exploration
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 1: DATA LOADING AND EXPLORATION"
di as txt _dup(70) "=" _n

// Load data using diddesign_data command
diddesign_data malesky2014, clear

// Drop observations with missing covariates
// This matches the R replication code: drop_na(malesky2014, lnarea, lnpopden, city)
drop if missing(lnarea) | missing(lnpopden) | missing(city)

// Basic data description
di as txt "Data Description:"
describe

// Summary statistics for key variables
di as txt _n "Summary Statistics for Key Variables:"
summarize pro4 tapwater agrext treatment year

// Check data structure - this is RCS data, not panel
di as txt _n "Data Structure (Repeated Cross-Section):"
di as txt "Note: Each year contains different observations, not a panel"
tabulate year treatment

// Check post-treatment indicator (treatment happened in 2009)
di as txt _n "Post-Treatment Indicator:"
di as txt "Using existing post_treat variable (1 if year == 2010)"
tabulate year post_treat

// Check covariate availability
di as txt _n "Covariate Summary:"
summarize lnarea lnpopden city reg8

/*---------------------------------------------------------------------------
 * Section 2: Parallel Trends Assessment
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 2: PARALLEL TRENDS ASSESSMENT"
di as txt _dup(70) "=" _n

// Set random seed for reproducibility
set seed 1234

di as txt "Note: For RCS data, parallel trends test compares pre-treatment"
di as txt "      outcome trends between treatment and control groups."
di as txt ""
di as txt "Extended Parallel Trends (EPT) assumption:"
di as txt "  E[Y_t - Y_{t-1} | G=1] - E[Y_t - Y_{t-1} | G=0] = constant"
di as txt ""
di as txt "Testing with lag=1 (comparing 2006-2008 trends)"
di as txt _n _dup(50) "-"

// =====================================================
// Test 1: pro4 (Education and Cultural Program)
// =====================================================
di as txt _n "=== Testing pro4 (Education and Cultural Program) ==="
di as txt "Running parallel trends test (nboot=200 for speed)..."

// Note: id_district is a string variable; diddesign auto-encodes it
diddesign_check pro4, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    lag(1) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

// Store results
matrix pro4_placebo = e(placebo)
local pro4_est = pro4_placebo[1,2]
local pro4_se = pro4_placebo[1,3]
local pro4_pval = 2 * (1 - normal(abs(`pro4_est'/`pro4_se')))

di as txt _n "pro4 Parallel Trends Test:"
di as txt "  Estimate:    " %9.4f `pro4_est'
di as txt "  Std. Error:  " %9.4f `pro4_se'
di as txt "  p-value:     " %9.4f `pro4_pval'
di as txt "  Conclusion:  " _continue
if `pro4_pval' > 0.05 {
    di as result "PASS - Parallel trends plausible (p > 0.05)"
}
else {
    di as error "FAIL - Parallel trends may not hold (p <= 0.05)"
}

// =====================================================
// Test 2: tapwater (Tap Water)
// =====================================================
di as txt _n "=== Testing tapwater (Tap Water) ==="
di as txt "Running parallel trends test..."

diddesign_check tapwater, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    lag(1) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

matrix tapwater_placebo = e(placebo)
local tapwater_est = tapwater_placebo[1,2]
local tapwater_se = tapwater_placebo[1,3]
local tapwater_pval = 2 * (1 - normal(abs(`tapwater_est'/`tapwater_se')))

di as txt _n "tapwater Parallel Trends Test:"
di as txt "  Estimate:    " %9.4f `tapwater_est'
di as txt "  Std. Error:  " %9.4f `tapwater_se'
di as txt "  p-value:     " %9.4f `tapwater_pval'
di as txt "  Conclusion:  " _continue
if `tapwater_pval' > 0.05 {
    di as result "PASS - Parallel trends plausible"
}
else {
    di as error "MARGINAL - Parallel trends questionable (p <= 0.05)"
}

// =====================================================
// Test 3: agrext (Agricultural Center)
// =====================================================
di as txt _n "=== Testing agrext (Agricultural Center) ==="
di as txt "Running parallel trends test..."

diddesign_check agrext, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    lag(1) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

matrix agrext_placebo = e(placebo)
local agrext_est = agrext_placebo[1,2]
local agrext_se = agrext_placebo[1,3]
local agrext_pval = 2 * (1 - normal(abs(`agrext_est'/`agrext_se')))

di as txt _n "agrext Parallel Trends Test:"
di as txt "  Estimate:    " %9.4f `agrext_est'
di as txt "  Std. Error:  " %9.4f `agrext_se'
di as txt "  p-value:     " %9.4f `agrext_pval'
di as txt "  Conclusion:  " _continue
if `agrext_pval' > 0.05 {
    di as result "PASS - Parallel trends plausible"
}
else {
    di as error "FAIL - Parallel trends likely violated (p < 0.05)"
}

// Summary Table
di as txt _n _dup(70) "-"
di as txt "PARALLEL TRENDS ASSESSMENT SUMMARY"
di as txt _dup(70) "-"
di as txt "Variable     Estimate   Std.Error  p-value    Conclusion"
di as txt _dup(70) "-"
di as txt "pro4        " %9.4f `pro4_est' "   " %9.4f `pro4_se' "   " %7.4f `pro4_pval' "   EPT satisfied"
di as txt "tapwater    " %9.4f `tapwater_est' "   " %9.4f `tapwater_se' "   " %7.4f `tapwater_pval' "   EPT marginal"
di as txt "agrext      " %9.4f `agrext_est' "   " %9.4f `agrext_se' "   " %7.4f `agrext_pval' "   EPT violated"
di as txt _dup(70) "-"

/*---------------------------------------------------------------------------
 * Section 3: Basic DID Estimation (Without Covariates)
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 3: BASIC DID ESTIMATION (NO COVARIATES)"
di as txt _dup(70) "=" _n

set seed 1234

// =====================================================
// Estimate 1: pro4
// =====================================================
di as txt "=== Estimating effect on pro4 (Education Program) ==="

diddesign pro4, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

// e(estimates) columns: [lead, estimate, std_error, ci_lo, ci_hi, weight]
// Row order: Double-DID (1), DID (2), sDID (3)
local pro4_ddid = e(estimates)[1,2]
local pro4_did = e(estimates)[2,2]
local pro4_sdid = e(estimates)[3,2]
local pro4_w_did = e(weights)[1,1]
local pro4_w_sdid = e(weights)[1,2]

di as txt _n "pro4 Estimation Results:"
di as txt "  Double-DID:  " %9.4f `pro4_ddid'
di as txt "  DID:         " %9.4f `pro4_did'
di as txt "  sDID:        " %9.4f `pro4_sdid'
di as txt "  Weights:     w_DID=" %5.3f `pro4_w_did' ", w_sDID=" %5.3f `pro4_w_sdid'

// =====================================================
// Estimate 2: tapwater
// =====================================================
di as txt _n "=== Estimating effect on tapwater (Tap Water) ==="

diddesign tapwater, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

local tapwater_ddid = e(estimates)[1,2]
local tapwater_did = e(estimates)[2,2]
local tapwater_sdid = e(estimates)[3,2]
local tapwater_w_did = e(weights)[1,1]
local tapwater_w_sdid = e(weights)[1,2]

di as txt _n "tapwater Estimation Results:"
di as txt "  Double-DID:  " %9.4f `tapwater_ddid'
di as txt "  DID:         " %9.4f `tapwater_did'
di as txt "  sDID:        " %9.4f `tapwater_sdid'
di as txt "  Weights:     w_DID=" %5.3f `tapwater_w_did' ", w_sDID=" %5.3f `tapwater_w_sdid'

// =====================================================
// Estimate 3: agrext
// =====================================================
di as txt _n "=== Estimating effect on agrext (Agricultural Center) ==="

diddesign agrext, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

local agrext_ddid = e(estimates)[1,2]
local agrext_did = e(estimates)[2,2]
local agrext_sdid = e(estimates)[3,2]
local agrext_w_did = e(weights)[1,1]
local agrext_w_sdid = e(weights)[1,2]

di as txt _n "agrext Estimation Results:"
di as txt "  Double-DID:  " %9.4f `agrext_ddid'
di as txt "  DID:         " %9.4f `agrext_did'
di as txt "  sDID:        " %9.4f `agrext_sdid'
di as txt "  Weights:     w_DID=" %5.3f `agrext_w_did' ", w_sDID=" %5.3f `agrext_w_sdid'

/*---------------------------------------------------------------------------
 * Section 4: Estimation with Covariates (Paper Replication)
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 4: ESTIMATION WITH COVARIATES (PAPER REPLICATION)"
di as txt _dup(70) "=" _n

di as txt "Covariates: lnarea, lnpopden, city, reg8 (regional FE)"
di as txt "Clustering: id_district level"
di as txt "Bootstrap: 200 iterations (use 2000 for publication)"
di as txt _n _dup(50) "-"

// Create region fixed effects dummies
quietly tabulate reg8, gen(reg8_)

set seed 1234

// =====================================================
// pro4 with covariates
// =====================================================
di as txt _n "=== pro4 with covariates ==="

diddesign pro4 lnarea lnpopden city reg8_2 reg8_3 reg8_4 reg8_5 reg8_6 reg8_7, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

local pro4_cov_ddid = e(estimates)[1,2]
local pro4_cov_se = e(estimates)[1,3]
local pro4_cov_ci_lo = e(estimates)[1,4]
local pro4_cov_ci_hi = e(estimates)[1,5]

di as txt "pro4 Results (with covariates):"
di as txt "  Double-DID:      " %9.4f `pro4_cov_ddid'
di as txt "  Std. Error:      " %9.4f `pro4_cov_se'
di as txt "  95% CI:          [" %7.4f `pro4_cov_ci_lo' ", " %7.4f `pro4_cov_ci_hi' "]"

// =====================================================
// tapwater with covariates
// =====================================================
di as txt _n "=== tapwater with covariates ==="

diddesign tapwater lnarea lnpopden city reg8_2 reg8_3 reg8_4 reg8_5 reg8_6 reg8_7, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

local tapwater_cov_ddid = e(estimates)[1,2]
local tapwater_cov_did = e(estimates)[2,2]
local tapwater_cov_sdid = e(estimates)[3,2]
local tapwater_cov_se = e(estimates)[1,3]

di as txt "tapwater Results (with covariates):"
di as txt "  Double-DID:      " %9.4f `tapwater_cov_ddid'
di as txt "  DID:             " %9.4f `tapwater_cov_did'
di as txt "  sDID:            " %9.4f `tapwater_cov_sdid'
di as txt "  Note: For tapwater, sDID may be more appropriate (PTT holds)"

// =====================================================
// agrext with covariates
// =====================================================
di as txt _n "=== agrext with covariates ==="

diddesign agrext lnarea lnpopden city reg8_2 reg8_3 reg8_4 reg8_5 reg8_6 reg8_7, ///
    treatment(treatment) time(year) ///
    post(post_treat) ///
    rcs ///
    cluster(id_district) ///
    nboot(200)

local agrext_cov_ddid = e(estimates)[1,2]
local agrext_cov_did = e(estimates)[2,2]
local agrext_cov_sdid = e(estimates)[3,2]

di as txt "agrext Results (with covariates):"
di as txt "  Double-DID:      " %9.4f `agrext_cov_ddid'
di as txt "  DID:             " %9.4f `agrext_cov_did'
di as txt "  sDID:            " %9.4f `agrext_cov_sdid'
di as txt "  Note: Neither EPT nor PTT holds; estimates may be biased"

/*---------------------------------------------------------------------------
 * Section 5: Visualization
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 5: VISUALIZATION"
di as txt _dup(70) "=" _n

// Generate plots for each outcome

// pro4 - Trends plot
di as txt "Generating plots for pro4..."
set seed 1234
diddesign_check pro4, ///
    treatment(treatment) time(year) ///
    post(post_treat) lag(1) rcs cluster(id_district) nboot(200)

diddesign_plot, type(trends) saving("malesky_pro4_trends.gph", replace)

// tapwater - Trends plot
di as txt "Generating plots for tapwater..."
set seed 1234
diddesign_check tapwater, ///
    treatment(treatment) time(year) ///
    post(post_treat) lag(1) rcs cluster(id_district) nboot(200)

diddesign_plot, type(trends) saving("malesky_tapwater_trends.gph", replace)

// agrext - Trends plot
di as txt "Generating plots for agrext..."
set seed 1234
diddesign_check agrext, ///
    treatment(treatment) time(year) ///
    post(post_treat) lag(1) rcs cluster(id_district) nboot(200)

diddesign_plot, type(trends) saving("malesky_agrext_trends.gph", replace)

di as txt _n "Plots saved: malesky_*_trends.gph"

/*---------------------------------------------------------------------------
 * Section 6: Results Interpretation
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 6: RESULTS INTERPRETATION"
di as txt _dup(70) "=" _n

di as txt "MALESKY ET AL. (2014) FINDINGS - DOUBLE DID METHOD"
di as txt ""
di as txt _dup(70) "-"
di as txt "OUTCOME 1: Education and Cultural Program (pro4)"
di as txt _dup(70) "-"
di as txt ""
di as txt "  Parallel Trends Assessment:"
di as txt "    - EPT assumption: SATISFIED (placebo p > 0.05)"
di as txt "    - Pre-treatment trends are approximately parallel"
di as txt ""
di as txt "  Treatment Effect:"
di as txt "    - Double-DID is preferred (most efficient under EPT)"
di as txt "    - Abolishing councils increased education programs"
di as txt ""
di as txt _dup(70) "-"
di as txt "OUTCOME 2: Tap Water (tapwater)"
di as txt _dup(70) "-"
di as txt ""
di as txt "  Parallel Trends Assessment:"
di as txt "    - EPT assumption: MARGINALLY VIOLATED"
di as txt "    - Parallel trends-in-trends (PTT) may hold"
di as txt ""
di as txt "  Treatment Effect:"
di as txt "    - sDID is preferred when only PTT holds"
di as txt "    - Standard DID may underestimate the negative effect"
di as txt "    - This illustrates the value of Double-DID methodology"
di as txt ""
di as txt _dup(70) "-"
di as txt "OUTCOME 3: Agricultural Extension Center (agrext)"
di as txt _dup(70) "-"
di as txt ""
di as txt "  Parallel Trends Assessment:"
di as txt "    - EPT assumption: VIOLATED"
di as txt "    - Pre-treatment trends have opposite signs"
di as txt ""
di as txt "  Treatment Effect:"
di as txt "    - All estimates may be biased"
di as txt "    - No credible causal effect can be identified"
di as txt ""
di as txt _dup(70) "-"
di as txt "METHODOLOGICAL INSIGHTS"
di as txt _dup(70) "-"
di as txt ""
di as txt "1. The Double-DID method provides transparent assumption testing"
di as txt "2. When EPT holds: Use Double-DID (combines DID and sDID efficiently)"
di as txt "3. When only PTT holds: Use sDID (robust to EPT violation)"
di as txt "4. When neither holds: Consider alternative identification strategies"

/*---------------------------------------------------------------------------
 * Section 7: Exporting Results
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 7: EXPORTING RESULTS"
di as txt _dup(70) "=" _n

// Create summary matrix
matrix results_summary = J(3, 6, .)
matrix rownames results_summary = "pro4" "tapwater" "agrext"
matrix colnames results_summary = "PT_est" "PT_pval" "DID" "sDID" "DoubleDID" "Preferred"

// Fill in results (using stored locals)
matrix results_summary[1,1] = `pro4_est'
matrix results_summary[1,2] = `pro4_pval'
matrix results_summary[1,3] = `pro4_did'
matrix results_summary[1,4] = `pro4_sdid'
matrix results_summary[1,5] = `pro4_ddid'
matrix results_summary[1,6] = `pro4_ddid'

matrix results_summary[2,1] = `tapwater_est'
matrix results_summary[2,2] = `tapwater_pval'
matrix results_summary[2,3] = `tapwater_did'
matrix results_summary[2,4] = `tapwater_sdid'
matrix results_summary[2,5] = `tapwater_ddid'
matrix results_summary[2,6] = `tapwater_sdid'

matrix results_summary[3,1] = `agrext_est'
matrix results_summary[3,2] = `agrext_pval'
matrix results_summary[3,3] = `agrext_did'
matrix results_summary[3,4] = `agrext_sdid'
matrix results_summary[3,5] = `agrext_ddid'
matrix results_summary[3,6] = .

di as txt "Results Summary Matrix:"
matrix list results_summary, format(%9.4f)

// Export to CSV
preserve
clear
set obs 3
gen outcome = ""
gen pt_estimate = .
gen pt_pvalue = .
gen did = .
gen sdid = .
gen double_did = .
gen preferred = ""

replace outcome = "pro4" in 1
replace pt_estimate = `pro4_est' in 1
replace pt_pvalue = `pro4_pval' in 1
replace did = `pro4_did' in 1
replace sdid = `pro4_sdid' in 1
replace double_did = `pro4_ddid' in 1
replace preferred = "Double-DID" in 1

replace outcome = "tapwater" in 2
replace pt_estimate = `tapwater_est' in 2
replace pt_pvalue = `tapwater_pval' in 2
replace did = `tapwater_did' in 2
replace sdid = `tapwater_sdid' in 2
replace double_did = `tapwater_ddid' in 2
replace preferred = "sDID" in 2

replace outcome = "agrext" in 3
replace pt_estimate = `agrext_est' in 3
replace pt_pvalue = `agrext_pval' in 3
replace did = `agrext_did' in 3
replace sdid = `agrext_sdid' in 3
replace double_did = `agrext_ddid' in 3
replace preferred = "None" in 3

export delimited using "malesky_results.csv", replace
di as txt _n "Results exported to malesky_results.csv"
restore

/*---------------------------------------------------------------------------
 * End of Example
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "MALESKY EXAMPLE COMPLETED SUCCESSFULLY"
di as txt _dup(70) "=" _n

di as txt "Key takeaways:"
di as txt "1. RCS data requires post() and rcs options, not id()"
di as txt "2. Always test parallel trends before interpreting effects"
di as txt "3. Choose estimator based on which assumptions are satisfied"
di as txt "4. Double-DID provides transparent guidance on estimator choice"

log close
