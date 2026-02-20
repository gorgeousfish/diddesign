/*===========================================================================
 * example_paglayan.do - Staggered Adoption Design Example
 *
 * Demonstrates the use of diddesign for staggered adoption (SA) designs,
 * where treatment timing varies across units. This example analyzes the
 * effect of collective bargaining requirements on education expenditure.
 *
 * Data Description:
 *   Panel data of US states from 1959-2000. States adopted collective
 *   bargaining requirements at different times (staggered adoption).
 *
 * Variables:
 *   - stateid: Numeric state identifier
 *   - state: State name (string)
 *   - year: Year (1959-2000)
 *   - treatment: Collective bargaining required (0/1)
 *   - pupil_expenditure: Per-pupil expenditure
 *   - teacher_salary: Average teacher salary
 *   - log_salary: Log of teacher salary (pre-computed)
 *
 * Key Notes:
 *   - States adopt treatment at different times (staggered)
 *   - Some states never adopt treatment (control group)
 *   - Use design(sa) for staggered adoption
 *   - SA design requires panel data (not RCS)
 *
 * Usage:
 *   After installing diddesign, download example files:
 *     . net get diddesign
 *   Then run from the examples/ directory:
 *     . cd examples
 *     . do example_paglayan.do
 *
 * Reference:
 *   Paglayan (2019). "Public-Sector Unions and the Size of Government."
 *   American Journal of Political Science.
 *===========================================================================*/

version 16
clear all
set more off
capture log close

log using "example_paglayan.log", replace

/*---------------------------------------------------------------------------
 * Section 1: Data Loading and Preparation
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 1: DATA LOADING AND PREPARATION"
di as txt _dup(70) "=" _n

// Load data using diddesign_data command
diddesign_data paglayan2019, clear

// Note: The R replication code excludes Wisconsin ("WI") and DC ("DC"):
//   paglayan2019 <- paglayan2019 %>% filter(!(state %in% c("WI", "DC")))
// The Stata dataset already excludes these states during data preparation.

// Basic data description
di as txt "Data Description:"
describe

// Summary statistics
di as txt _n "Summary Statistics:"
summarize treatment pupil_expenditure teacher_salary year

// Create log-transformed outcome variable
di as txt _n "Creating Log-Transformed Outcome Variables..."
gen log_expenditure = log(pupil_expenditure + 1)
label variable log_expenditure "Log per-pupil expenditure"

// Recompute log_salary to match R replication code: log(teacher_salary + 1)
// The dataset contains a pre-computed log_salary that is numerically equivalent,
// but we recompute for exact consistency with the R code.
drop log_salary
gen log_salary = log(teacher_salary + 1)
label variable log_salary "Log teacher salary"

// Note: stateid (numeric state identifier) already exists in dataset

// Summary of outcome variables
di as txt _n "Outcome Variables:"
summarize log_expenditure log_salary

/*---------------------------------------------------------------------------
 * Section 2: Treatment Pattern Exploration
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 2: TREATMENT PATTERN EXPLORATION"
di as txt _dup(70) "=" _n

// Treatment adoption summary
di as txt "Treatment Adoption by Year:"
tabulate year treatment

// First treatment year by state
preserve
collapse (min) first_treat_year = year if treatment == 1, by(state)
di as txt _n "First Treatment Year by State:"
list, sep(0)
restore

// Count of treated states over time
preserve
collapse (sum) n_treated = treatment, by(year)
di as txt _n "Number of Treated States by Year:"
list, sep(0)
restore

/*---------------------------------------------------------------------------
 * Section 3: Parallel Trends Assessment for SA Design
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 3: SA PARALLEL TRENDS ASSESSMENT"
di as txt _dup(70) "=" _n

// Set random seed for reproducibility
set seed 1234

// Check parallel trends with multiple lag periods
// Note: SA design aggregates across different treatment timing groups
di as txt "Running SA parallel trends diagnostic..."
di as txt "This may take several minutes due to bootstrap..." _n

// Note: treatment goes in treatment() option, not in the varlist
diddesign_check log_expenditure, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) ///
    lag(1 2 3 4 5) ///
    thres(1) ///
    nboot(200)

// Display results
di as txt _n "SA Placebo Test Results:"
matrix list e(placebo), format(%9.4f)

// Interpretation for SA design
di as txt _n "Interpretation for Staggered Adoption:"
di as txt "- Results are time-weighted averages across treatment timing groups"
di as txt "- Each lag represents pre-treatment period relative to adoption"
di as txt "- thres(1) includes groups with at least 1 treated unit"

/*---------------------------------------------------------------------------
 * Section 4: SA Double DID Estimation
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 4: SA DOUBLE DID ESTIMATION"
di as txt _dup(70) "=" _n

set seed 1234

// Main SA estimation
di as txt "Running SA Double DID estimation..."
diddesign log_expenditure, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) ///
    thres(1) ///
    nboot(200)

// Display detailed results
// e(estimates) columns: [lead, estimate, std_error, ci_lo, ci_hi, weight]
// Row order: SA-Double-DID (row 1), SA-DID (row 2), SA-sDID (row 3)
di as txt _n "SA Estimation Results:"
di as txt _dup(50) "-"
di as txt "  SA-Double DID estimate: " %9.4f e(estimates)[1,2]
di as txt "  SA-DID estimate:        " %9.4f e(estimates)[2,2]
di as txt "  SA-sDID estimate:       " %9.4f e(estimates)[3,2]
di as txt _dup(50) "-"
di as txt "  Weight on SA-DID:       " %9.4f e(weights)[1,1]
di as txt "  Weight on SA-sDID:      " %9.4f e(weights)[1,2]
di as txt _dup(50) "-"

ereturn list

/*---------------------------------------------------------------------------
 * Section 5: Dynamic Effects with Multiple Lead Periods
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 5: SA DYNAMIC TREATMENT EFFECTS"
di as txt _dup(70) "=" _n

set seed 1234

// Estimate effects at multiple lead periods (0 through 9)
di as txt "Estimating dynamic treatment effects (lead 0-9)..."
diddesign log_expenditure, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) ///
    lead(0 1 2 3 4 5 6 7 8 9) ///
    thres(1) ///
    nboot(200)

di as txt _n "Dynamic Treatment Effects (SA Design):"
matrix list e(estimates), format(%9.4f)

// Interpretation
di as txt _n "Dynamic Effects Interpretation:"
di as txt "- lead(0): Effect at time of treatment adoption"
di as txt "- lead(k): Effect k periods after adoption"
di as txt "- Pattern shows how treatment effect evolves over time"

/*---------------------------------------------------------------------------
 * Section 6: Sensitivity to Threshold Parameter
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 6: SENSITIVITY TO THRESHOLD PARAMETER"
di as txt _dup(70) "=" _n

di as txt "The thres() parameter controls minimum treated units per timing group"
di as txt ""

// thres = 1 (include all groups)
set seed 1234
di as txt "thres(1) - Include all timing groups:"
diddesign log_expenditure, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) thres(1) nboot(100)
local est_thres1 = e(estimates)[1,2]
di as txt "  Double DID: " %9.4f `est_thres1'

// thres = 2 (require at least 2 treated units)
set seed 1234
di as txt _n "thres(2) - Require >= 2 treated units per group:"
diddesign log_expenditure, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) thres(2) nboot(100)
local est_thres2 = e(estimates)[1,2]
di as txt "  Double DID: " %9.4f `est_thres2'

// thres = 3
set seed 1234
di as txt _n "thres(3) - Require >= 3 treated units per group:"
diddesign log_expenditure, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) thres(3) nboot(100)
local est_thres3 = e(estimates)[1,2]
di as txt "  Double DID: " %9.4f `est_thres3'

di as txt _n "Sensitivity Summary:"
di as txt "  thres(1): " %9.4f `est_thres1'
di as txt "  thres(2): " %9.4f `est_thres2'
di as txt "  thres(3): " %9.4f `est_thres3'

/*---------------------------------------------------------------------------
 * Section 7: Visualization for SA Design
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 7: SA VISUALIZATION"
di as txt _dup(70) "=" _n

// First run diagnostic check for plot data
set seed 1234
diddesign_check log_expenditure, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) lag(1 2 3 4 5) thres(1) nboot(200)

// Treatment pattern plot (unique to SA design)
di as txt "Generating treatment pattern plot..."
diddesign_plot, type(pattern) saving("paglayan_pattern.gph", replace)

// Placebo plot
di as txt "Generating SA placebo plot..."
diddesign_plot, type(placebo) saving("paglayan_placebo.gph", replace)

// Combined plot
di as txt "Generating combined plot..."
diddesign_plot, type(both) saving("paglayan_combined.gph", replace)

/*---------------------------------------------------------------------------
 * Section 8: Alternative Outcome - Teacher Salary
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 8: ALTERNATIVE OUTCOME - TEACHER SALARY"
di as txt _dup(70) "=" _n

set seed 1234

di as txt "Estimating effect on log teacher salary..."
diddesign log_salary, ///
    treatment(treatment) ///
    id(stateid) time(year) ///
    design(sa) ///
    thres(1) ///
    nboot(200)

di as txt _n "Effect on Log Teacher Salary:"
di as txt "  SA-Double DID: " %9.4f e(estimates)[1,2]
di as txt "  SA-DID:        " %9.4f e(estimates)[2,2]
di as txt "  SA-sDID:       " %9.4f e(estimates)[3,2]

/*---------------------------------------------------------------------------
 * Section 9: Results Interpretation for SA Design
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SECTION 9: SA DESIGN INTERPRETATION"
di as txt _dup(70) "=" _n

di as txt "STAGGERED ADOPTION INTERPRETATION GUIDE:"
di as txt ""
di as txt "1. What SA-DID Measures:"
di as txt "   - Time-weighted average treatment effect"
di as txt "   - Weights based on number of treated units at each timing"
di as txt "   - Aggregates across different treatment timing groups"
di as txt ""
di as txt "2. Control Group Definition:"
di as txt "   - For each treated unit, control group consists of:"
di as txt "     a) Never-treated units"
di as txt "     b) Not-yet-treated units (treated later)"
di as txt ""
di as txt "3. Time Weights:"
di as txt "   - Weights proportional to number of treated units"
di as txt "   - Larger groups contribute more to the estimate"
di as txt ""
di as txt "4. Threshold Parameter (thres):"
di as txt "   - Sets minimum treated units for inclusion"
di as txt "   - Higher thres = more precision, less coverage"
di as txt "   - Lower thres = more coverage, less precision"
di as txt ""
di as txt "5. Dynamic Effects:"
di as txt "   - lead(k) shows effect k periods after adoption"
di as txt "   - Can reveal delayed or building treatment effects"
di as txt ""
di as txt "6. Key Assumptions:"
di as txt "   - Parallel trends (tested via placebo)"
di as txt "   - No anticipation effects"
di as txt "   - SUTVA (no spillovers)"

/*---------------------------------------------------------------------------
 * End of Example
 *---------------------------------------------------------------------------*/

di as txt _n _dup(70) "="
di as txt "SA EXAMPLE COMPLETED SUCCESSFULLY"
di as txt _dup(70) "=" _n

di as txt "Key takeaways:"
di as txt "1. SA design handles staggered treatment adoption"
di as txt "2. Use design(sa) with id() and time() for SA analysis"
di as txt "3. thres() controls precision-coverage tradeoff"
di as txt "4. Dynamic effects via lead() trace post-adoption trajectories"

log close
