*! did_estimators.mata - Core estimation functions for difference-in-differences
*!
*! Implements point estimation for the standard DID estimator (under the parallel
*! trends assumption) and the sequential DID estimator (under the weaker parallel
*! trends-in-trends assumption). Also provides data structures and equivalence
*! confidence interval computation for assessing pre-treatment trends.

version 16.0

mata:
mata set matastrict on

// ----------------------------------------------------------------------------
// DATA STRUCTURES
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * struct did_data - Data container for DID estimation
 *
 * Stores outcome, treatment, and derived variables for difference-in-
 * differences estimation. Supports panel and repeated cross-sectional designs.
 *---------------------------------------------------------------------------*/
struct did_data {
    // Original variables
    real colvector outcome       // Y_it: outcome
    real colvector treatment     // D_it: treatment indicator (0/1)
    real colvector id_unit       // Unit identifier i
    real colvector id_time       // Time period t (normalized)
    real matrix    covariates    // X_it: covariates (optional)
    real colvector cluster_var   // Cluster identifier (optional)
    
    // Derived variables
    real colvector Gi            // G_i: treatment group (1) vs control (0)
    real colvector It            // I_t: post-treatment indicator (0/1)
    real colvector id_time_std   // Standardized time (0 = treatment period)
    real colvector outcome_delta // ΔY_it: first-differenced outcome
    
    // Metadata
    real scalar    N             // Number of observations
    real scalar    n_units       // Number of units
    real scalar    n_periods     // Number of time periods
    real scalar    treat_year    // Treatment period (standardized)
    real scalar    is_panel      // 1 = panel, 0 = repeated cross-section
}

/*---------------------------------------------------------------------------
 * struct did_option - Estimation options
 *
 * User-specified options for bootstrap, variance estimation, and display.
 *---------------------------------------------------------------------------*/
struct did_option {
    real scalar    n_boot        // Bootstrap replications
    real scalar    parallel      // 1 = parallel computing enabled
    real scalar    se_boot       // 1 = bootstrap CI, 0 = analytical
    string scalar  id_cluster    // Cluster variable name
    real rowvector lead          // Post-treatment period indices
    real scalar    thres         // Staggered adoption threshold
    real rowvector lag           // Pre-treatment period indices
    real scalar    stdz          // 1 = standardize estimates
    real scalar    level         // Confidence level (percent)
    real scalar    seed          // Random seed (. if unset)
    string scalar  var_cluster_pre  // Internal: original cluster variable
    real scalar    quiet         // 1 = suppress progress display
}

// ----------------------------------------------------------------------------
// INITIALIZATION FUNCTIONS
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * init_did_option() - Initialize did_option with default values
 *
 * Returns:
 *   struct did_option: initialized with default values
 *---------------------------------------------------------------------------*/
struct did_option scalar init_did_option()
{
    struct did_option scalar opt
    
    opt.n_boot         = 30
    opt.parallel       = 1
    opt.se_boot        = 0
    opt.id_cluster     = ""
    opt.lead           = 0
    opt.thres          = 2
    opt.lag            = 1
    opt.stdz           = 1
    opt.level          = 95
    opt.seed           = .
    opt.var_cluster_pre = ""
    opt.quiet          = 0
    
    return(opt)
}

// ----------------------------------------------------------------------------
// OLS HELPER FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * ols_coef() - Extract OLS coefficient via normal equations
 *
 * Computes β = (X'X)⁻¹X'y and returns the coefficient at specified index.
 *
 * Arguments:
 *   X        : real matrix (n × p), design matrix
 *   y        : real colvector (n × 1), outcome
 *   coef_idx : real scalar, coefficient index (1-based)
 *
 * Returns:
 *   real scalar: coefficient value, or missing if singular/invalid
 *---------------------------------------------------------------------------*/
real scalar ols_coef(real matrix X, real colvector y, real scalar coef_idx)
{
    real matrix XtX, XtX_inv
    real colvector Xty, beta
    real scalar n, p, n_valid
    real colvector valid_idx
    
    n = rows(X)
    p = cols(X)
    
    // Check minimum observations
    if (n < p) {
        return(.)
    }
    
    // Listwise deletion for missing values
    valid_idx = selectindex(rowmissing(X) :== 0 :& y :< .)
    n_valid = length(valid_idx)
    
    // Need at least p observations for OLS
    if (n_valid < p) {
        return(.)
    }
    
    // Use cross() for efficient computation on valid observations only
    XtX = cross(X[valid_idx, .], X[valid_idx, .])
    Xty = cross(X[valid_idx, .], y[valid_idx])
    
    // Rank check for numerical stability
    if (rank(XtX) < p) {
        return(.)
    }
    
    // Symmetric matrix inversion
    XtX_inv = invsym(XtX)
    
    // Singularity check
    if (missing(XtX_inv[1,1])) {
        return(.)
    }
    if (any(diagonal(XtX_inv) :== 0)) {
        return(.)
    }
    
    // Compute coefficients
    beta = XtX_inv * Xty
    
    // Validate coefficient index
    if (coef_idx < 1 | coef_idx > length(beta)) {
        return(.)
    }
    
    return(beta[coef_idx])
}


// ----------------------------------------------------------------------------
// MAIN DID ESTIMATION FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * did_fit() - Compute DID and sequential DID point estimates
 *
 * The standard DID estimator τ̂_DID is consistent under the parallel trends
 * assumption. The sequential DID estimator τ̂_sDID requires only the weaker
 * parallel trends-in-trends assumption.
 *
 * Arguments:
 *   Y        : real colvector, outcome Y_it
 *   Y_delta  : real colvector, first-differenced outcome ΔY_it
 *   Gi       : real colvector, group indicator G_i (1 = treated)
 *   It       : real colvector, post-treatment indicator I_t
 *   X        : real matrix, covariates (or empty)
 *   time_std : real colvector, standardized time index
 *   lead     : real scalar, post-treatment period (0 = treatment year)
 *
 * Returns:
 *   real rowvector (1 × 2): (τ̂_DID, τ̂_sDID)
 *
 * Regression models:
 *   DID:  Y_it = α + β₁G_i + β₂I_t + τG_i×I_t + X'γ + ε
 *   sDID: ΔY_it = α + β₁G_i + β₂I_t + τG_i×I_t + X'γ + ε
 *---------------------------------------------------------------------------*/
real rowvector did_fit(real colvector Y, real colvector Y_delta,
                       real colvector Gi, real colvector It,
                       real matrix X, real colvector time_std,
                       real scalar lead)
{
    real colvector idx, idx_did, idx_sdid
    real colvector Y_sub, Yd_sub, Gi_sub, It_sub
    real colvector Gi_did, It_did, Gi_sdid, It_sdid
    real colvector valid_did, valid_sdid, cov_valid, gi_it_valid
    real matrix X_sub, X_did, X_sdid, design_did, design_sdid
    real rowvector result
    real scalar tau_did, tau_sdid, n_sub, k_cov
    
    // -------------------------------------------------------------------------
    // Input validation
    // -------------------------------------------------------------------------
    if (lead < 0) {
        errprintf("{err}Error: lead must be >= 0 (got %g)\n", lead)
        return((., .))
    }
    
    // Subset data to relevant time periods: t ∈ {-1, lead}
    idx = selectindex((time_std :== -1) :| (time_std :== lead))
    
    // Return missing if no observations in specified periods
    if (length(idx) == 0) {
        return((., .))
    }
    
    // Extract observations for the selected time periods
    Y_sub = Y[idx]
    Yd_sub = Y_delta[idx]
    Gi_sub = Gi[idx]
    It_sub = It[idx]
    n_sub = length(idx)
    
    // Handle covariates
    k_cov = cols(X)
    if (k_cov > 0) {
        X_sub = X[idx, .]
    }
    
    // Listwise deletion for missing values
    valid_did = (Y_sub :< .)
    valid_sdid = (Yd_sub :< .)
    
    // Exclude observations with missing covariates
    if (k_cov > 0) {
        cov_valid = (rowmissing(X_sub) :== 0)
        valid_did = valid_did :& cov_valid
        valid_sdid = valid_sdid :& cov_valid
    }
    
    // Exclude observations with missing group or time indicators
    gi_it_valid = (Gi_sub :< .) :& (It_sub :< .)
    valid_did = valid_did :& gi_it_valid
    valid_sdid = valid_sdid :& gi_it_valid
    
    // Obtain indices of valid observations
    idx_did = selectindex(valid_did)
    idx_sdid = selectindex(valid_sdid)
    
    // Standard DID estimation (requires at least 4 observations)
    if (length(idx_did) < 4) {
        tau_did = .
    }
    else {
        Gi_did = Gi_sub[idx_did]
        It_did = It_sub[idx_did]
        
        // Check for sufficient variation
        if (min(Gi_did) == max(Gi_did)) {
            tau_did = .
        }
        else if (min(It_did) == max(It_did)) {
            tau_did = .
        }
        else {
            // Construct design matrix: [1, G_i, I_t, G_i×I_t, X]
            design_did = J(length(idx_did), 1, 1), Gi_did, It_did, Gi_did :* It_did
            
            if (k_cov > 0) {
                X_did = X_sub[idx_did, .]
                design_did = design_did, X_did
            }
            
            // Extract coefficient on interaction term (position 4)
            tau_did = ols_coef(design_did, Y_sub[idx_did], 4)
        }
    }
    
    // Sequential DID estimation
    if (length(idx_sdid) < 4) {
        tau_sdid = .
    }
    else {
        Gi_sdid = Gi_sub[idx_sdid]
        It_sdid = It_sub[idx_sdid]
        
        // Check for sufficient variation
        if (min(Gi_sdid) == max(Gi_sdid)) {
            tau_sdid = .
        }
        else if (min(It_sdid) == max(It_sdid)) {
            tau_sdid = .
        }
        else {
            // Construct design matrix: [1, G_i, I_t, G_i×I_t, X]
            design_sdid = J(length(idx_sdid), 1, 1), Gi_sdid, It_sdid, Gi_sdid :* It_sdid
            
            if (k_cov > 0) {
                X_sdid = X_sub[idx_sdid, .]
                design_sdid = design_sdid, X_sdid
            }
            
            // Extract coefficient on interaction term (position 4)
            tau_sdid = ols_coef(design_sdid, Yd_sub[idx_sdid], 4)
        }
    }
    
    // Return point estimates
    result = (tau_did, tau_sdid)
    return(result)
}


/*---------------------------------------------------------------------------
 * did_fit_struct() - Wrapper for did_fit() using did_data structure
 *
 * Arguments:
 *   data : struct did_data, prepared data
 *   lead : real scalar, post-treatment period (default: 0)
 *
 * Returns:
 *   real rowvector (1 × 2): (τ̂_DID, τ̂_sDID)
 *---------------------------------------------------------------------------*/
real rowvector did_fit_struct(struct did_data scalar data, | real scalar lead)
{
    // Default lead = 0
    if (args() < 2) lead = 0
    
    // Extract data from structure and call did_fit()
    return(did_fit(
        data.outcome,
        data.outcome_delta,
        data.Gi,
        data.It,
        data.covariates,
        data.id_time_std,
        lead
    ))
}

// ----------------------------------------------------------------------------
// OPTION POPULATION FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _diddesign_populate_option() - Populate global did_option structure
 *
 * Internal function called from _diddesign_parse.ado to transfer parsed
 * command-line options to the global did_opt structure.
 *
 * Returns:
 *   0 on success
 *---------------------------------------------------------------------------*/
real scalar _diddesign_populate_option(
    real scalar n_boot,
    real scalar parallel,
    real scalar se_boot,
    string scalar id_cluster,
    real rowvector lead,
    real scalar thres,
    real rowvector lag,
    real scalar stdz,
    real scalar level,
    real scalar seed
)
{
    external struct did_option scalar did_opt
    
    did_opt = init_did_option()
    did_opt.n_boot         = n_boot
    did_opt.parallel       = parallel
    did_opt.se_boot        = se_boot
    did_opt.id_cluster     = id_cluster
    did_opt.lead           = lead
    did_opt.thres          = thres
    did_opt.lag            = lag
    did_opt.stdz           = stdz
    did_opt.level          = level
    did_opt.seed           = seed
    did_opt.var_cluster_pre = ""
    
    return(0)
}

// ----------------------------------------------------------------------------
// DATA POPULATION FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _diddesign_populate_data() - Populate global did_data structure
 *
 * Internal function called from _diddesign_prep.ado to transfer prepared
 * data from Stata variables to the global did_dat structure.
 *
 * Returns:
 *   0 on success
 *---------------------------------------------------------------------------*/
real scalar _diddesign_populate_data(
    string scalar outcome_var,
    string scalar treatment_var,
    string scalar id_var,
    string scalar id_time_var,
    string scalar covar_vars,
    string scalar cluster_var,
    string scalar Gi_var,
    string scalar It_var,
    string scalar id_time_std_var,
    string scalar outcome_delta_var,
    real scalar N,
    real scalar n_units,
    real scalar n_periods,
    real scalar treat_year,
    real scalar is_panel,
    string scalar touse_var
)
{
    external struct did_data scalar did_dat
    string rowvector covar_list
    real scalar k
    
    did_dat = did_data()
    
    // Original variables
    did_dat.outcome = st_data(., outcome_var, touse_var)
    did_dat.treatment = st_data(., treatment_var, touse_var)
    
    if (id_var != "") {
        did_dat.id_unit = st_data(., id_var, touse_var)
    }
    else {
        did_dat.id_unit = J(N, 1, .)
    }
    
    did_dat.id_time = st_data(., id_time_var, touse_var)
    
    // Covariates
    if (covar_vars != "") {
        covar_list = tokens(covar_vars)
        did_dat.covariates = st_data(., covar_list, touse_var)
    }
    else {
        did_dat.covariates = J(0, 0, .)
    }
    
    // Cluster variable
    if (cluster_var != "") {
        did_dat.cluster_var = st_data(., cluster_var, touse_var)
    }
    else {
        did_dat.cluster_var = J(0, 1, .)
    }
    
    // Derived variables
    did_dat.Gi = st_data(., Gi_var, touse_var)
    did_dat.It = st_data(., It_var, touse_var)
    did_dat.id_time_std = st_data(., id_time_std_var, touse_var)
    did_dat.outcome_delta = st_data(., outcome_delta_var, touse_var)
    
    // Metadata
    did_dat.N = N
    did_dat.n_units = n_units
    did_dat.n_periods = n_periods
    did_dat.treat_year = treat_year
    did_dat.is_panel = is_panel
    
    return(0)
}

// ----------------------------------------------------------------------------
// EQUIVALENCE CONFIDENCE INTERVAL FUNCTIONS
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * compute_eq_ci() - Compute 95% equivalence confidence interval
 *
 * Computes the equivalence CI for assessing parallel trends using the
 * Two One-Sided Tests (TOST) methodology.
 *
 * Arguments:
 *   estimate  : real scalar, point estimate (e.g., placebo DID)
 *   std_error : real scalar, standard error
 *
 * Returns:
 *   real rowvector (1 × 2): symmetric equivalence interval (-ν, ν)
 *
 * Method:
 *   1. Construct 90% CI: estimate ± z_{0.95} × SE
 *   2. Compute ν = max(|CI90_UB|, |CI90_LB|)
 *   3. Return (-ν, ν)
 *---------------------------------------------------------------------------*/
real rowvector compute_eq_ci(real scalar estimate, real scalar std_error)
{
    real scalar z_95, CI90_UB, CI90_LB, CI90_UB_ab, CI90_LB_ab, nu
    
    // Input validation: return missing for invalid inputs
    if (missing(estimate) | missing(std_error) | std_error <= 0) {
        return((., .))
    }
    
    // z critical value for 90% CI (corresponds to 95% equivalence test)
    // invnormal(0.95) ≈ 1.64485362695147
    z_95 = invnormal(0.95)
    
    // Compute 90% CI bounds
    CI90_UB = estimate + z_95 * std_error
    CI90_LB = estimate - z_95 * std_error
    
    // Take absolute values
    CI90_UB_ab = abs(CI90_UB)
    CI90_LB_ab = abs(CI90_LB)
    
    // Symmetric equivalence bound: max of absolute values
    nu = max((CI90_UB_ab, CI90_LB_ab))
    
    // Return symmetric equivalence CI: (-nu, nu)
    return((-nu, nu))
}


/*---------------------------------------------------------------------------
 * compute_eq_ci_vec() - Vectorized equivalence CI computation
 *
 * Batch computation of equivalence confidence intervals.
 *
 * Arguments:
 *   estimates  : real colvector (n × 1), point estimates
 *   std_errors : real colvector (n × 1), standard errors
 *
 * Returns:
 *   real matrix (n × 2): each row is (-ν, ν) equivalence interval
 *---------------------------------------------------------------------------*/
real matrix compute_eq_ci_vec(real colvector estimates, real colvector std_errors)
{
    real scalar n, z_95, i
    real colvector CI90_UB, CI90_LB, CI90_UB_ab, CI90_LB_ab, nu
    real colvector invalid_mask
    real matrix result
    
    n = rows(estimates)
    
    // Dimension check
    if (rows(std_errors) != n) {
        _error("estimates and std_errors must have the same number of rows")
    }
    
    if (n == 0) {
        return(J(0, 2, .))
    }
    
    // 90% CI corresponds to 95% equivalence test
    z_95 = invnormal(0.95)
    
    // Vectorized computation
    CI90_UB = estimates :+ z_95 :* std_errors
    CI90_LB = estimates :- z_95 :* std_errors
    CI90_UB_ab = abs(CI90_UB)
    CI90_LB_ab = abs(CI90_LB)
    nu = rowmax((CI90_UB_ab, CI90_LB_ab))
    
    // Mark invalid entries
    invalid_mask = (estimates :>= .) :| (std_errors :>= .) :| (std_errors :<= 0)
    
    // Set nu to missing for invalid rows
    for (i = 1; i <= n; i++) {
        if (invalid_mask[i]) {
            nu[i] = .
        }
    }
    
    // Build result matrix: (-nu, nu)
    result = (-nu, nu)
    
    return(result)
}


// ----------------------------------------------------------------------------
// MODULE VERIFICATION FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _did_estimators_loaded() - Verify module is loaded
 *---------------------------------------------------------------------------*/
void _did_estimators_loaded()
{
    printf("{txt}did_estimators.mata loaded successfully\n")
}

end
