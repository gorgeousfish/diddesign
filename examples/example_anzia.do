/*===========================================================================
 * example_anzia.do - Standard DID Example with Panel Data
 *
 * Replication of Anzia (2012) analysis using diddesign
 *
 * Data Description:
 *   Panel data of Texas school districts, examining the effect of moving
 *   school board elections from off-cycle to on-cycle timing on teacher
 *   salaries.
 *
 * Variables:
 *   - district: School district identifier
 *   - year: Year (2003-2009)
 *   - lnavgsalary_cpi: Log average teacher salary (CPI-adjusted) [outcome]
 *   - oncycle: On-cycle election indicator (0/1) [treatment]
 *   - teachers_avg_yrs_exper: Average years of teacher experience [covariate]
 *   - ami_pc, asian_pc, black_pc, hisp_pc: Demographics [covariates]
 *
 * Key Notes:
 *   - This is PANEL data with the same districts observed across years
 *   - Use id() option to specify unit identifier
 *   - Treatment is staggered: districts switch to on-cycle at different times
 *
 * Usage:
 *   After installing diddesign, download example files:
 *     . net get diddesign
 *   Then run from the examples/ directory:
 *     . cd examples
 *     . do example_anzia.do
 *
 * Reference:
 *   Anzia, S. F. (2012). "The Election Timing Effect: Evidence from a
 *   Policy Intervention in Texas." Quarterly Journal of Political Science.
 *===========================================================================*/

version 16
clear all
set more off
capture log close

log using "example_anzia.log", replace

/*---------------------------------------------------------------------------
 * Section 1: Data Loading and Exploration
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 1: DATA LOADING AND EXPLORATION"
di as txt _dup(70) "=" _n

// Load data
diddesign_data anzia2012, clear

// Basic data description
di as txt "Data Description:"
describe

// Summary statistics for key variables
di as txt _n "Summary Statistics for Key Variables:"
summarize lnavgsalary_cpi oncycle teachers_avg_yrs_exper

// Check panel structure
di as txt _n "Panel Structure:"
xtset district year
xtdescribe

// Treatment timing
di as txt _n "Treatment Distribution by Year:"
tabulate year oncycle

/*---------------------------------------------------------------------------
 * Section 2: Parallel Trends Assessment
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 2: PARALLEL TRENDS ASSESSMENT"
di as txt _dup(70) "=" _n

set seed 1234

// Check parallel trends with multiple lag periods
diddesign_check lnavgsalary_cpi, ///
    treatment(oncycle) ///
    id(district) time(year) ///
    lag(1 2 3) ///
    nboot(200)

// Display placebo test results
di as txt _n "Placebo Test Results:"
matrix list e(placebo), format(%9.4f)

// Visualize diagnostics
diddesign_plot, type(both) saving("anzia_diagnostics.gph", replace)

/*---------------------------------------------------------------------------
 * Section 3: Basic Double DID Estimation
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 3: BASIC DOUBLE DID ESTIMATION"
di as txt _dup(70) "=" _n

set seed 1234

diddesign lnavgsalary_cpi, ///
    treatment(oncycle) ///
    id(district) time(year) ///
    nboot(200)

// e(estimates) columns: [lead, estimate, std_error, ci_lo, ci_hi, weight]
// Row order: Double-DID (1), DID (2), sDID (3)
di as txt _n "Estimation Results:"
di as txt _dup(50) "-"
di as txt "  Double DID estimate: " %9.4f e(estimates)[1,2]
di as txt "  DID estimate:        " %9.4f e(estimates)[2,2]
di as txt "  sDID estimate:       " %9.4f e(estimates)[3,2]
di as txt _dup(50) "-"
di as txt "  Weight on DID:       " %9.4f e(weights)[1,1]
di as txt "  Weight on sDID:      " %9.4f e(weights)[1,2]
di as txt _dup(50) "-"

ereturn list


/*---------------------------------------------------------------------------
 * Section 4: Estimation with Covariates
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 4: ESTIMATION WITH COVARIATES"
di as txt _dup(70) "=" _n

set seed 1234

// Covariates: teacher experience, demographics
diddesign lnavgsalary_cpi ///
    teachers_avg_yrs_exper ami_pc asian_pc black_pc hisp_pc, ///
    treatment(oncycle) ///
    id(district) time(year) ///
    nboot(200)

di as txt _n "Results with Covariates:"
di as txt "  Double DID estimate: " %9.4f e(estimates)[1,2]
di as txt "  DID estimate:        " %9.4f e(estimates)[2,2]
di as txt "  sDID estimate:       " %9.4f e(estimates)[3,2]

/*---------------------------------------------------------------------------
 * Section 5: Dynamic Treatment Effects
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 5: DYNAMIC TREATMENT EFFECTS"
di as txt _dup(70) "=" _n

set seed 1234

// Estimate effects at multiple lead periods (0, 1, 2)
diddesign lnavgsalary_cpi, ///
    treatment(oncycle) ///
    id(district) time(year) ///
    lead(0 1 2) ///
    nboot(200)

di as txt _n "Dynamic Treatment Effects:"
matrix list e(estimates), format(%9.4f)

// Plot treatment effect estimates
diddesign_plot, type(estimates) saving("anzia_effects.gph", replace)

/*---------------------------------------------------------------------------
 * Section 6: Visualization
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 6: VISUALIZATION"
di as txt _dup(70) "=" _n

// Run diagnostic check to generate plot data
set seed 1234
diddesign_check lnavgsalary_cpi, ///
    treatment(oncycle) ///
    id(district) time(year) ///
    lag(1 2 3) ///
    nboot(200)

// Trends plot
diddesign_plot, type(trends) saving("anzia_trends.gph", replace)

// Placebo plot with equivalence CI
diddesign_plot, type(placebo) saving("anzia_placebo.gph", replace)

/*---------------------------------------------------------------------------
 * Section 7: Results Summary
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 7: RESULTS SUMMARY"
di as txt _dup(70) "=" _n

di as txt "Anzia (2012) - Effect of on-cycle elections on teacher salaries"
di as txt ""
di as txt "Key findings:"
di as txt "1. Placebo tests support the parallel trends assumption"
di as txt "2. Moving to on-cycle elections is associated with higher salaries"
di as txt "3. Results are robust across DID, sDID, and Double-DID estimators"
di as txt "4. Covariates do not substantially change the estimates"

/*---------------------------------------------------------------------------
 * End of Example
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "ANZIA EXAMPLE COMPLETED SUCCESSFULLY"
di as txt _dup(70) "=" _n

log close
