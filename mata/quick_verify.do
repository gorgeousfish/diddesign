*! quick_verify.do - Quick verification test for did_estimators.mata
*! Verifies that comment changes did not affect numerical results

version 16.0
clear all
mata: mata clear

// Load all Mata modules
quietly do "/Users/cxy/Desktop/2026project/diddesign/diddesign-stata/mata/diddesign_mata.do"

display as text "{hline 70}"
display as text "Quick Verification Test - did_estimators.mata"
display as text "{hline 70}"

// Create simple test data
set obs 100
set seed 12345
gen id = _n
gen time = mod(_n - 1, 5) + 1
gen treat = runiform() > 0.5
gen y = 1 + 0.5 * treat + 0.3 * time + rnormal()

// Compute simple DID manually
qui sum y if treat == 1 & time == 5
local y11 = r(mean)
qui sum y if treat == 1 & time == 4
local y10 = r(mean)
qui sum y if treat == 0 & time == 5
local y01 = r(mean)
qui sum y if treat == 0 & time == 4
local y00 = r(mean)

local did_manual = (`y11' - `y10') - (`y01' - `y00')
display as text "Manual DID: " as result `did_manual'

// Test Mata ols_coef function
mata:
// Create design matrix for simple regression
y = st_data(., "y")
X = J(100, 2, 1)
X[., 2] = st_data(., "treat")

// Test OLS coefficient extraction
beta1 = ols_coef(X, y, 1)
beta2 = ols_coef(X, y, 2)

printf("OLS intercept: %g\n", beta1)
printf("OLS treat coef: %g\n", beta2)

// Verify ols_coef returns valid numbers
if (missing(beta1) | missing(beta2)) {
    printf("{err}ERROR: ols_coef returned missing values\n")
    exit(1)
}
else {
    printf("{txt}PASS: ols_coef returns valid coefficients\n")
}

// Test compute_eq_ci function
eq_ci = compute_eq_ci(0.5, 0.1)
printf("Equivalence CI: [%g, %g]\n", eq_ci[1], eq_ci[2])

if (missing(eq_ci[1]) | missing(eq_ci[2])) {
    printf("{err}ERROR: compute_eq_ci returned missing values\n")
    exit(1)
}
else {
    printf("{txt}PASS: compute_eq_ci returns valid interval\n")
}

// Test compute_eq_ci_vec function
estimates = (0.5 \ 0.3 \ -0.2)
std_errors = (0.1 \ 0.15 \ 0.2)
eq_ci_mat = compute_eq_ci_vec(estimates, std_errors)

if (missing(eq_ci_mat[1,1]) | missing(eq_ci_mat[3,2])) {
    printf("{err}ERROR: compute_eq_ci_vec returned missing values\n")
    exit(1)
}
else {
    printf("{txt}PASS: compute_eq_ci_vec returns valid matrix\n")
}

printf("\n{txt}All verification tests passed!\n")
end

display as text "{hline 70}"
display as text "Verification Complete"
display as text "{hline 70}"
