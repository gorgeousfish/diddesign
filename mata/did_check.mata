*! did_check.mata - Parallel trends diagnostic functions
*!
*! Implements placebo tests for assessing the parallel trends assumption in
*! difference-in-differences designs. Supports both standard DID and staggered
*! adoption designs. Provides cluster-bootstrap inference and equivalence
*! confidence intervals for pre-treatment trend evaluation.

version 16.0

mata:
mata set matastrict on

// ============================================================================
// Parallel Trends Diagnostic Module
// ============================================================================
//
// This module implements placebo tests for evaluating the parallel trends
// assumption in difference-in-differences designs. Under parallel trends,
// pre-treatment DID estimates should be approximately zero; significant
// deviations indicate potential violations.
//
// Placebo test procedure:
//   For each pre-treatment lag l:
//     1. Subset data to periods {-l, -l-1}
//     2. Define pseudo-treatment indicator It = 1{time >= -l}
//     3. Estimate DID coefficient (expected zero under parallel trends)
//     4. Standardize by control group baseline standard deviation
//
// Core functionality:
//   - Standard DID placebo tests (did_placebo)
//   - Staggered adoption placebo tests (did_sad_placebo)
//   - Cluster-bootstrap standard errors (did_placebo_boot_full)
//   - Equivalence confidence intervals for trend assessment
//
// ============================================================================

// ----------------------------------------------------------------------------
// DATA STRUCTURES
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * struct placebo_result
 *
 * placebo test point estimates.
 *---------------------------------------------------------------------------*/
struct placebo_result {
    real colvector est           // Raw estimates
    real colvector est_std       // Standardized estimates
    real colvector lags          // Feasible lag values
}

/*---------------------------------------------------------------------------
 * struct placebo_boot_result
 *
 * bootstrap results including standard errors.
 *---------------------------------------------------------------------------*/
struct placebo_boot_result {
    real colvector se            // Bootstrap standard errors (raw)
    real colvector se_std        // Bootstrap standard errors (standardized)
    real scalar n_valid          // Number of successful iterations
    real matrix boot_est         // Bootstrap estimates (raw)
    real matrix boot_est_std     // Bootstrap estimates (standardized)
}

/*---------------------------------------------------------------------------
 * struct check_result
 *
 * Complete diagnostic check results including standard errors and
 * equivalence confidence intervals for parallel trends assessment.
 *---------------------------------------------------------------------------*/
struct check_result {
    real colvector lag           // Lag values
    real colvector estimate      // Standardized estimates
    real colvector estimate_raw  // Raw estimates
    real colvector std_error     // Standard errors (standardized)
    real colvector std_error_raw // Standard errors (raw)
    real colvector eq_ci_low     // Equivalence CI lower bounds
    real colvector eq_ci_high    // Equivalence CI upper bounds
}

// ----------------------------------------------------------------------------
// CORE FUNCTIONS
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * did_placebo() - Compute Placebo DID Estimates
 *
 * Computes placebo DID estimates for pre-treatment periods to assess
 * the parallel trends assumption. Under parallel trends, these estimates
 * should be approximately zero.
 *
 * Arguments:
 *   Y        : real colvector - outcome variable
 *   Gi       : real colvector - treatment group indicator 
 *   time_std : real colvector - standardized time (0 = treatment)
 *   X        : real matrix - covariates 
 *   lags     : real rowvector - lag periods to test 
 *   stdz     : real scalar - standardization flag
 *
 * Returns:
 *   struct placebo_result with raw and standardized estimates
 *
 * Algorithm:
 *   1. Filter infeasible lags: keep only lags < max_lag
 *      where max_lag = abs(min(time_std))
 *   
 *   For each valid lag l:
 *     2. Filter data to periods {-l, -l-1}
 *     3. Create pseudo-treatment indicator It = 1{time_std >= -l}
 *     4. Compute raw DID estimate via OLS
 *     5. Standardize by control group baseline SD:
 *        outcome_std = (outcome - mean(control)) / sd(control)
 *     6. Compute standardized DID estimate
 *---------------------------------------------------------------------------*/
struct placebo_result scalar did_placebo(real colvector Y,
                                         real colvector Gi,
                                         real colvector time_std,
                                         real matrix X,
                                         real rowvector lags,
                                         real scalar stdz)
{
    struct placebo_result scalar result
    real rowvector lags_abs, valid_lags
    real colvector idx, Y_use, Gi_use, time_use, It
    real colvector Y_std, ct_idx
    real matrix X_use, design
    real scalar max_lag, n_lags, i, lag, ct_mean, ct_sd
    real scalar est_raw, est_standardized, n_use, k_cov
    
    // Lag feasibility filtering
    lags_abs = abs(lags)
    max_lag = abs(min(time_std))
    
    // Warn if lag=0 specified
    if (any(lags :== 0)) {
        printf("{txt}Warning: lag=0 tests treatment period vs pre-period, not a true placebo test\n")
    }
    
    // Keep only feasible lags
    valid_lags = select(lags_abs, lags_abs :< max_lag)
    n_lags = cols(valid_lags)
    
    // Handle empty case
    if (n_lags == 0) {
        result.est = J(0, 1, .)
        result.est_std = J(0, 1, .)
        result.lags = J(0, 1, .)
        return(result)
    }
    
    // Initialize result containers
    result.est = J(n_lags, 1, .)
    result.est_std = J(n_lags, 1, .)
    result.lags = valid_lags'
    
    k_cov = cols(X)
    
    // -------------------------------------------------------------------------
    // Main loop over lags
    // -------------------------------------------------------------------------
    for (i = 1; i <= n_lags; i++) {
        lag = valid_lags[i]
        
        // Subset to periods {-lag, -lag-1}
        idx = selectindex((time_std :== -lag) :| (time_std :== -lag - 1))
        
        if (length(idx) == 0) {
            continue
        }
        
        // Extract subset
        Y_use = Y[idx]
        Gi_use = Gi[idx]
        time_use = time_std[idx]
        n_use = length(idx)
        
        if (k_cov > 0) {
            X_use = X[idx, .]
        }
        else {
            X_use = J(n_use, 0, .)
        }
        
        // Define pseudo-treatment indicator
        It = (time_use :>= -lag)
        
        // Estimate raw DID
        est_raw = _ols_did_coef(Y_use, Gi_use, It, X_use)
        result.est[i] = est_raw
        
        // Standardization
        if (stdz) {
            ct_idx = selectindex((It :== 0) :& (Gi_use :== 0))
            
            if (length(ct_idx) == 0) {
                result.est_std[i] = .
                continue
            }
            
            // Exclude missing values from control baseline
            real colvector ct_vals, ct_vals_valid
            ct_vals = Y_use[ct_idx]
            ct_vals_valid = select(ct_vals, ct_vals :< .)
            
            // Require at least 2 observations for sample variance
            if (length(ct_vals_valid) < 2) {
                result.est_std[i] = .
                continue
            }
            
            // Compute control baseline statistics
            ct_mean = mean(ct_vals_valid)
            ct_sd = sqrt(variance(ct_vals_valid))
            
            if (ct_sd == 0 | missing(ct_sd)) {
                result.est_std[i] = .
                continue
            }
            
            // Standardize outcome
            Y_std = (Y_use :- ct_mean) :/ ct_sd
            
            // Estimate standardized DID
            est_standardized = _ols_did_coef(Y_std, Gi_use, It, X_use)
            result.est_std[i] = est_standardized
        }
    }
    
    return(result)
}

/*---------------------------------------------------------------------------
 * _ols_did_coef() - Extract DID Coefficient via OLS
 *
 * Constructs design matrix [1, Gi, It, Gi*It, X] and extracts the
 * interaction coefficient. Listwise deletion handles missing values.
 *
 * Arguments:
 *   Y  : real colvector - outcome variable
 *   Gi : real colvector - group indicator
 *   It : real colvector - time indicator
 *   X  : real matrix - covariates 
 *
 * Returns:
 *   real scalar: coefficient on Gi*It interaction term
 *---------------------------------------------------------------------------*/
real scalar _ols_did_coef(real colvector Y, real colvector Gi, 
                          real colvector It, real matrix X)
{
    real matrix design, XtX, XtX_inv
    real colvector Xty, beta, valid_idx
    real scalar n, n_valid, p, k_cov
    
    n = rows(Y)
    k_cov = cols(X)
    
    // Listwise deletion: exclude missing observations
    valid_idx = selectindex((Y :< .) :& (Gi :< .) :& (It :< .))
    
    if (k_cov > 0 && length(valid_idx) > 0) {
        valid_idx = select(valid_idx, rowmissing(X[valid_idx, .]) :== 0)
    }
    
    n_valid = length(valid_idx)
    
    if (n_valid == 0) {
        return(.)
    }
    
    // Construct design matrix: [1, Gi, It, Gi*It, X]
    design = J(n_valid, 1, 1), Gi[valid_idx], It[valid_idx], Gi[valid_idx] :* It[valid_idx]
    
    if (k_cov > 0) {
        design = design, X[valid_idx, .]
    }
    
    p = cols(design)
    
    if (n_valid < p) {
        return(.)
    }
    
    // OLS estimation: beta = (X'X)^(-1) X'y
    XtX = cross(design, design)
    Xty = cross(design, Y[valid_idx])
    
    if (rank(XtX) < p) {
        return(.)
    }
    
    XtX_inv = invsym(XtX)
    
    if (missing(XtX_inv[1,1])) {
        return(.)
    }
    if (any(diagonal(XtX_inv) :== 0)) {
        return(.)
    }
    
    beta = XtX_inv * Xty
    
    // Return Gi*It coefficient
    return(beta[4])
}

/*---------------------------------------------------------------------------
 * did_placebo_boot() - Single Bootstrap Iteration for Placebo Tests
 *
 * Performs one cluster bootstrap iteration for placebo test inference.
 *
 * Arguments:
 *   data        : struct did_data - data structure
 *   cluster_ids : real colvector - cluster identifiers
 *   cluster_var : real colvector - cluster membership
 *   lags        : real rowvector - lag parameters
 *   is_panel    : real scalar - data type indicator
 *
 * Returns:
 *   struct placebo_result with bootstrap estimates
 *
 * Algorithm:
 *   1. Sample clusters with replacement
 *   2. Construct bootstrap dataset preserving within-cluster structure
 *   3. Renumber unit IDs in bootstrap sample
 *   4. Compute placebo estimates on bootstrap sample
 *---------------------------------------------------------------------------*/
struct placebo_result scalar did_placebo_boot(struct did_data scalar data,
                                              real colvector cluster_ids,
                                              real colvector cluster_var,
                                              real rowvector lags,
                                              real scalar is_panel)
{
    struct placebo_result scalar result
    real colvector id_boot, idx, new_id_unit
    real colvector Y_boot, Gi_boot, time_std_boot
    real matrix X_boot
    real scalar n_clusters, i, j, k, n_obs, k_cov
    
    n_clusters = rows(cluster_ids)
    k_cov = cols(data.covariates)
    
    // Sample clusters with replacement
    id_boot = cluster_ids[safe_sample_idx(n_clusters, n_clusters)]
    
    // Count total observations
    n_obs = 0
    for (j = 1; j <= n_clusters; j++) {
        idx = selectindex(cluster_var :== id_boot[j])
        n_obs = n_obs + length(idx)
    }
    
    // Allocate bootstrap arrays
    Y_boot = J(n_obs, 1, .)
    Gi_boot = J(n_obs, 1, .)
    time_std_boot = J(n_obs, 1, .)
    new_id_unit = J(n_obs, 1, .)
    if (k_cov > 0) {
        X_boot = J(n_obs, k_cov, .)
    }
    else {
        X_boot = J(0, 0, .)
    }
    
    // Second pass: fill bootstrap data
    k = 1
    for (j = 1; j <= n_clusters; j++) {
        idx = selectindex(cluster_var :== id_boot[j])
        
        if (length(idx) > 0) {
            // Copy data for this cluster
            Y_boot[|k \ k + length(idx) - 1|] = data.outcome[idx]
            Gi_boot[|k \ k + length(idx) - 1|] = data.Gi[idx]
            time_std_boot[|k \ k + length(idx) - 1|] = data.id_time_std[idx]
            
            // Renumber id_unit for bootstrap sample
            new_id_unit[|k \ k + length(idx) - 1|] = J(length(idx), 1, j)
            
            // Copy covariates if present
            if (k_cov > 0) {
                X_boot[|k, 1 \ k + length(idx) - 1, k_cov|] = data.covariates[idx, .]
            }
            
            k = k + length(idx)
        }
    }
    
    // Compute placebo estimates on bootstrap sample
    result = did_placebo(Y_boot, Gi_boot, time_std_boot, X_boot, lags, 1)
    
    return(result)
}

/*---------------------------------------------------------------------------
 * did_placebo_boot_full() - Complete Bootstrap SE Computation
 *
 * Performs n_boot bootstrap iterations and computes standard errors.
 * Failed iterations are tracked and excluded from variance computation.
 *
 * Arguments:
 *   data     : struct did_data - data structure
 *   lags     : real rowvector - lag parameters
 *   n_boot   : real scalar - number of bootstrap iterations
 *   is_panel : real scalar - data type indicator
 *   cluster  : string scalar - cluster variable name (optional)
 *
 * Returns:
 *   struct placebo_boot_result with standard errors and estimates
 *---------------------------------------------------------------------------*/
struct placebo_boot_result scalar did_placebo_boot_full(
    struct did_data scalar data,
    real rowvector lags,
    real scalar n_boot,
    real scalar is_panel,
    string scalar cluster)
{
    struct placebo_boot_result scalar result
    struct placebo_result scalar boot_est
    real matrix boot_est_mat, boot_est_std_mat
    real colvector valid_idx, cluster_ids, cluster_var
    real colvector valid_rows, col_data
    real scalar b, n_valid, n_lags, i, lag_idx
    
    // Handle cluster variable
    if (cluster == "" & is_panel) {
        cluster_var = data.id_unit
    }
    else if (cluster == "") {
        cluster_var = (1::rows(data.outcome))
    }
    else {
        cluster_var = data.cluster_var
    }
    
    cluster_ids = uniqrows(cluster_var)
    
    // Pre-allocate result matrices
    n_lags = cols(lags)
    boot_est_mat = J(n_boot, n_lags, .)
    boot_est_std_mat = J(n_boot, n_lags, .)
    valid_idx = J(n_boot, 1, 0)
    
    // Bootstrap loop
    for (b = 1; b <= n_boot; b++) {
        boot_est = did_placebo_boot(data, cluster_ids, cluster_var, lags, is_panel)
        
        if (rows(boot_est.lags) > 0) {
            for (i = 1; i <= rows(boot_est.lags); i++) {
                lag_idx = _find_lag_position(lags, boot_est.lags[i])
                if (lag_idx > 0) {
                    boot_est_mat[b, lag_idx] = boot_est.est[i]
                    boot_est_std_mat[b, lag_idx] = boot_est.est_std[i]
                }
            }
            valid_idx[b] = 1
        }
    }
    
    // Remove invalid iterations
    n_valid = sum(valid_idx)
    
    if (n_valid > 0) {
        valid_rows = selectindex(valid_idx)
        boot_est_mat = boot_est_mat[valid_rows, .]
        boot_est_std_mat = boot_est_std_mat[valid_rows, .]
    }
    else {
        // All iterations failed
        boot_est_mat = J(0, n_lags, .)
        boot_est_std_mat = J(0, n_lags, .)
    }
    
    // Compute standard errors (using n-1 denominator)
    result.se = J(n_lags, 1, .)
    result.se_std = J(n_lags, 1, .)
    
    if (n_valid > 1) {
        for (i = 1; i <= n_lags; i++) {
            col_data = boot_est_mat[., i]
            col_data = select(col_data, col_data :< .)  // Remove missing values
            if (rows(col_data) > 1) {
                result.se[i] = sqrt(variance(col_data))
            }
            
            col_data = boot_est_std_mat[., i]
            col_data = select(col_data, col_data :< .)  // Remove missing values
            if (rows(col_data) > 1) {
                result.se_std[i] = sqrt(variance(col_data))
            }
        }
    }
    
    result.n_valid = n_valid
    result.boot_est = boot_est_mat
    result.boot_est_std = boot_est_std_mat
    
    return(result)
}

/*---------------------------------------------------------------------------
 * _find_lag_position() - Find position of lag value in lags vector
 *
 * Maps bootstrap results to lag positions.
 *
 * Arguments:
 *   lags    : real rowvector, original lag values
 *   lag_val : real scalar, lag value to find
 *
 * Returns:
 *   real scalar: position (1-indexed), or 0 if not found
 *---------------------------------------------------------------------------*/
real scalar _find_lag_position(real rowvector lags, real scalar lag_val)
{
    real scalar i, n, tol
    
    // Use tolerance comparison for floating-point robustness
    tol = 1e-10
    
    n = cols(lags)
    for (i = 1; i <= n; i++) {
        if (abs(abs(lags[i]) - abs(lag_val)) < tol) {
            return(i)
        }
    }
    return(0)
}

/*---------------------------------------------------------------------------
 * did_sad_placebo() - Staggered Adoption Design Placebo Tests
 *
 * Computes time-weighted placebo estimates aggregated across treatment
 * cohorts for staggered adoption designs. Infeasible periods are excluded
 * and weights are renormalized to sum to unity.
 *
 * Arguments:
 *   data   : struct did_data - panel data structure
 *   option : struct did_option - estimation options
 *
 * Returns:
 *   struct sa_placebo_result containing:
 *     - estimates[,1]: standardized placebo estimates
 *     - estimates[,2]: raw placebo estimates
 *     - Gmat: treatment timing matrix
 *
 * Algorithm:
 *   1. Construct Gmat (treatment timing indicator matrix)
 *   2. Identify valid treatment periods via threshold criterion
 *   3. For each valid period t:
 *      a. Subset to units treated at t and their controls
 *      b. Compute period-specific placebo estimates
 *   4. Aggregate using time weights with renormalization
 *---------------------------------------------------------------------------*/
struct sa_placebo_result scalar did_sad_placebo(struct did_data scalar data,
                                                 struct did_option scalar option)
{
    struct sa_placebo_result scalar result
    struct placebo_result scalar placebo_tmp
    real matrix Gmat, est_did, est_did_std
    real colvector id_time_use, time_weight
    pointer vector id_subj_use
    real scalar n_periods, n_lags, i, j, lag_idx
    real colvector tmp, tmp_std, w_use, valid_idx
    real colvector Y_use, Gi_use, time_std_use
    real matrix X_use
    real colvector idx, idx_subj
    real scalar t, n_use
    
    // Initialize result
    n_lags = cols(option.lag)
    result.estimates = J(n_lags, 2, .)
    result.valid_lags = option.lag
    
    // Create Gmat (group indicator matrix)
    Gmat = create_gmat(data.id_unit, data.id_time, data.treatment)
    result.Gmat = Gmat
    
    if (rows(Gmat) == 0 || cols(Gmat) == 0) {
        return(result)
    }
    
    // Get valid periods
    id_time_use = get_periods(Gmat, option.thres)
    
    if (rows(id_time_use) == 0) {
        return(result)
    }
    
    // Get valid subjects for each period
    id_subj_use = get_subjects(Gmat, id_time_use)
    
    // Get time weights
    time_weight = get_time_weight(Gmat, id_time_use)
    
    n_periods = rows(id_time_use)
    
    // Initialize period-specific estimate matrices
    est_did = J(n_periods, n_lags, .)
    est_did_std = J(n_periods, n_lags, .)
    
    // For each valid period, compute placebo estimates
    for (i = 1; i <= n_periods; i++) {
        t = id_time_use[i]
        idx_subj = *id_subj_use[i]
        
        // Subset data: units in id_subj_use[i], times <= t
        idx = _sa_placebo_subset_idx(data, idx_subj, t)
        
        if (length(idx) == 0) {
            continue
        }
        
        Y_use = data.outcome[idx]
        n_use = rows(Y_use)
        
        // Compute Gi and id_time_std for subset
        _sa_placebo_compute_Gi_time_std(data, idx, idx_subj, t, Gmat, &Gi_use, &time_std_use)
        
        if (cols(data.covariates) > 0) {
            X_use = data.covariates[idx, .]
        }
        else {
            X_use = J(n_use, 0, .)
        }
        
        // Run placebo regression
        placebo_tmp = did_placebo(Y_use, Gi_use, time_std_use, X_use, option.lag, option.stdz)
        
        // Store results (handle infeasible lags via lag name matching)
        for (j = 1; j <= rows(placebo_tmp.lags); j++) {
            lag_idx = _find_lag_position(option.lag, placebo_tmp.lags[j])
            if (lag_idx > 0) {
                est_did[i, lag_idx] = placebo_tmp.est[j]
                est_did_std[i, lag_idx] = placebo_tmp.est_std[j]
            }
        }
    }
    
    // Compute time-weighted average with missing value handling
    // Standardized and raw estimates use separate valid indices
    // (missing patterns may differ when standardization fails but raw estimation succeeds)
    real colvector valid_idx_std, valid_idx_raw, w_use_std, w_use_raw
    real scalar w_sum_std, w_sum_raw
    
    for (j = 1; j <= n_lags; j++) {
        tmp = est_did[., j]
        tmp_std = est_did_std[., j]
        
        // Standardized estimates
        valid_idx_std = selectindex(tmp_std :< .)
        if (rows(valid_idx_std) > 0) {
            w_use_std = time_weight[valid_idx_std]
            w_sum_std = sum(w_use_std)
            if (w_sum_std > 0) {
                w_use_std = w_use_std / w_sum_std
                result.estimates[j, 1] = sum(tmp_std[valid_idx_std] :* w_use_std)
            }
        }
        
        // Raw estimates
        valid_idx_raw = selectindex(tmp :< .)
        if (rows(valid_idx_raw) > 0) {
            w_use_raw = time_weight[valid_idx_raw]
            w_sum_raw = sum(w_use_raw)
            if (w_sum_raw > 0) {
                w_use_raw = w_use_raw / w_sum_raw
                result.estimates[j, 2] = sum(tmp[valid_idx_raw] :* w_use_raw)
            }
        }
    }
    
    return(result)
}

/*---------------------------------------------------------------------------
 * _sa_placebo_subset_idx() - Subset Indices for Staggered Adoption
 *
 * Returns observation indices satisfying:
 *   - Unit is in the valid subject set for this period
 *   - Time is at or before the current period t
 *
 * Arguments:
 *   data     : struct did_data - full panel data
 *   idx_subj : real colvector - valid unit indices (rows in Gmat)
 *   t        : real scalar - current period (column in Gmat)
 *
 * Returns:
 *   real colvector: observation row indices
 *---------------------------------------------------------------------------*/
real colvector _sa_placebo_subset_idx(struct did_data scalar data,
                                       real colvector idx_subj,
                                       real scalar t)
{
    real colvector units, valid_units, idx
    real scalar N, i
    real colvector is_valid_unit, is_valid_time
    transmorphic scalar valid_set
    
    N = rows(data.outcome)
    
    // Get unique unit IDs
    units = uniqrows(data.id_unit)
    
    // Validate idx_subj bounds before array access
    if (rows(idx_subj) > 0) {
        if (max(idx_subj) > rows(units) || min(idx_subj) < 1) {
            errprintf("Error: _sa_placebo_subset_idx(): idx_subj contains out-of-bounds indices\n")
            errprintf("       idx_subj range: [%g, %g], units count: %g\n", 
                      min(idx_subj), max(idx_subj), rows(units))
            return(J(0, 1, .))
        }
    }
    
    // Get valid unit IDs
    valid_units = units[idx_subj]
    
    // Build valid_units set for O(1) lookup
    valid_set = asarray_create("real", 1)
    for (i = 1; i <= rows(valid_units); i++) {
        asarray(valid_set, valid_units[i], 1)
    }
    
    // Create indicator for valid units in O(N)
    is_valid_unit = J(N, 1, 0)
    for (i = 1; i <= N; i++) {
        if (asarray_contains(valid_set, data.id_unit[i])) {
            is_valid_unit[i] = 1
        }
    }
    
    // Create indicator for valid times (time <= t)
    is_valid_time = (data.id_time :<= t)
    
    // Return indices where both conditions are met
    idx = selectindex(is_valid_unit :& is_valid_time)
    
    return(idx)
}

/*---------------------------------------------------------------------------
 * _sa_placebo_compute_Gi_time_std() - Compute Group and Time Indicators
 *
 * For the data subset, computes treatment group indicator and
 * standardized time relative to treatment period.
 *
 * Arguments:
 *   data         : struct did_data - full panel data
 *   idx          : real colvector - observation row indices
 *   idx_subj     : real colvector - valid unit indices (rows in Gmat)
 *   t            : real scalar - current treatment period
 *   Gmat         : real matrix - treatment timing indicator matrix
 *   Gi           : pointer(real colvector) - output group indicator
 *   id_time_std  : pointer(real colvector) - output standardized time
 *
 * Output:
 *   Gi = 1 if unit is newly treated at t, 0 if control
 *   id_time_std = id_time - t (time relative to treatment)
 *---------------------------------------------------------------------------*/
void _sa_placebo_compute_Gi_time_std(struct did_data scalar data,
                                      real colvector idx,
                                      real colvector idx_subj,
                                      real scalar t,
                                      real matrix Gmat,
                                      pointer(real colvector) scalar Gi,
                                      pointer(real colvector) scalar id_time_std)
{
    real scalar n_obs, i, u, unit_idx
    real colvector units, valid_units
    transmorphic scalar unit_idx_map
    
    n_obs = rows(idx)
    *Gi = J(n_obs, 1, .)
    *id_time_std = J(n_obs, 1, .)
    
    // Get unique units and valid units
    units = uniqrows(data.id_unit)
    valid_units = units[idx_subj]
    
    // Build unit index map for O(1) lookup
    unit_idx_map = asarray_create("real", 1)
    for (i = 1; i <= rows(units); i++) {
        asarray(unit_idx_map, units[i], i)
    }
    
    for (i = 1; i <= n_obs; i++) {
        u = data.id_unit[idx[i]]
        
        // Find unit index in Gmat using asarray (O(1) lookup)
        if (asarray_contains(unit_idx_map, u)) {
            unit_idx = asarray(unit_idx_map, u)
            
            if (unit_idx > 0 && unit_idx <= rows(Gmat)) {
                // Gi = 1 if Gmat[unit, t] == 1 (newly treated at t)
                // Gi = 0 if Gmat[unit, t] == 0 (control)
                (*Gi)[i] = (Gmat[unit_idx, t] == 1) ? 1 : 0
                
                // id_time_std = id_time - t
                (*id_time_std)[i] = data.id_time[idx[i]] - t
            }
        }
    }
}

/*---------------------------------------------------------------------------
 * did_sad_placebo_boot() - Bootstrap SE for Staggered Adoption Placebo
 *
 * Computes cluster-bootstrap standard errors for staggered adoption
 * placebo tests. Failed iterations are excluded from variance computation.
 *
 * Arguments:
 *   data   : struct did_data - panel data structure
 *   option : struct did_option - estimation options including n_boot
 *
 * Returns:
 *   struct sa_placebo_boot_result containing:
 *     - se_std, se_orig: bootstrap standard errors
 *     - boot_est_std, boot_est_orig: bootstrap estimate matrices
 *     - n_valid: count of successful iterations
 *
 * Algorithm:
 *   1. For each bootstrap iteration:
 *      a. Sample units with replacement
 *      b. Compute staggered adoption placebo estimates
 *      c. Validate and store results
 *   2. Compute SE as sample standard deviation of valid estimates
 *---------------------------------------------------------------------------*/
struct sa_placebo_boot_result scalar did_sad_placebo_boot(
    struct did_data scalar data,
    struct did_option scalar option)
{
    struct sa_placebo_boot_result scalar result
    struct sa_placebo_result scalar boot_result
    struct did_data scalar boot_data
    real matrix boot_est_std, boot_est_orig
    real scalar n_boot, n_lags, b, j
    real colvector col_data
    
    n_boot = option.n_boot
    n_lags = cols(option.lag)
    
    // Initialize result
    result.n_boot = n_boot
    result.se_std = J(n_lags, 1, .)
    result.se_orig = J(n_lags, 1, .)
    
    // Pre-allocate bootstrap estimate matrices
    boot_est_std = J(n_boot, n_lags, .)
    boot_est_orig = J(n_boot, n_lags, .)
    
    // Bootstrap loop with validation
    real scalar valid_count, progress_freq
    valid_count = 0
    
    // Progress display frequency
    progress_freq = max((1, floor(n_boot / 10)))
    
    for (b = 1; b <= n_boot; b++) {
        // Progress display (controlled by quiet option)
        if (option.quiet == 0 && mod(b, progress_freq) == 0) {
            printf("{txt}Bootstrap: %g/%g (%g%%)\n", b, n_boot, round(100*b/n_boot))
            displayflush()
        }
        
        // Sample panel data with replacement (unit-level)
        boot_data = sample_panel(data)
        
        // Compute staggered adoption placebo estimates on bootstrap sample
        boot_result = did_sad_placebo(boot_data, option)
        
        // Validate result: check ALL lags have valid estimates
        real scalar all_lags_valid
        all_lags_valid = 1
        if (rows(boot_result.estimates) == n_lags && cols(boot_result.estimates) >= 2) {
            for (j = 1; j <= n_lags; j++) {
                if (missing(boot_result.estimates[j, 1]) || 
                    missing(boot_result.estimates[j, 2])) {
                    all_lags_valid = 0
                    break
                }
            }
        }
        else {
            all_lags_valid = 0
        }
        
        if (all_lags_valid) {
            boot_est_std[b, .] = boot_result.estimates[., 1]'  // Standardized
            boot_est_orig[b, .] = boot_result.estimates[., 2]' // Original
            valid_count++
        }
    }
    
    // Final progress display
    if (option.quiet == 0) {
        printf("{txt}Bootstrap: %g/%g (100%%)\n", n_boot, n_boot)
        displayflush()
    }
    
    // Store valid count and warn if some iterations failed
    result.n_valid = valid_count
    
    if (valid_count < n_boot & option.quiet == 0) {
        printf("{txt}Warning: %g of %g staggered adoption placebo bootstrap iterations failed\n",
               n_boot - valid_count, n_boot)
    }
    
    // Remove invalid rows from bootstrap matrices
    if (valid_count > 0 && valid_count < n_boot) {
        real colvector valid_idx, valid_rows
        real scalar row_valid
        valid_idx = J(n_boot, 1, 0)
        for (b = 1; b <= n_boot; b++) {
            row_valid = 1
            for (j = 1; j <= n_lags; j++) {
                if (boot_est_std[b, j] >= . || boot_est_orig[b, j] >= .) {
                    row_valid = 0
                    break
                }
            }
            valid_idx[b] = row_valid
        }
        valid_rows = selectindex(valid_idx)
        if (rows(valid_rows) > 0) {
            boot_est_std = boot_est_std[valid_rows, .]
            boot_est_orig = boot_est_orig[valid_rows, .]
        }
        else {
            boot_est_std = J(0, n_lags, .)
            boot_est_orig = J(0, n_lags, .)
        }
    }
    else if (valid_count == 0) {
        // All iterations failed
        boot_est_std = J(0, n_lags, .)
        boot_est_orig = J(0, n_lags, .)
    }
    
    // Store bootstrap estimates
    result.boot_est_std = boot_est_std
    result.boot_est_orig = boot_est_orig
    
    // Compute standard errors (using n-1 denominator)
    for (j = 1; j <= n_lags; j++) {
        col_data = boot_est_std[., j]
        col_data = select(col_data, col_data :< .)  // Remove missing values
        if (rows(col_data) > 1) {
            result.se_std[j] = sqrt(variance(col_data))
        }
        
        col_data = boot_est_orig[., j]
        col_data = select(col_data, col_data :< .)  // Remove missing values
        if (rows(col_data) > 1) {
            result.se_orig[j] = sqrt(variance(col_data))
        }
    }
    
    return(result)
}

// ----------------------------------------------------------------------------
// MODULE VERIFICATION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _did_check_loaded() - Module Load Verification
 *---------------------------------------------------------------------------*/
void _did_check_loaded()
{
    printf("{txt}did_check.mata loaded successfully\n")
}

/*---------------------------------------------------------------------------
 * _test_did_placebo() - Test Wrapper for did_placebo()
 *
 * Wrapper function for interactive testing that returns a matrix
 * instead of a struct for easier inspection.
 *
 * Arguments:
 *   Y        : real colvector - outcome variable
 *   Gi       : real colvector - group indicator
 *   time_std : real colvector - standardized time index
 *   lags     : real rowvector - lag periods to test
 *
 * Returns:
 *   real matrix (n_lags x 3): columns are [lag, est, est_std]
 *---------------------------------------------------------------------------*/
real matrix _test_did_placebo(real colvector Y, real colvector Gi,
                              real colvector time_std, real rowvector lags)
{
    struct placebo_result scalar res
    real matrix output
    real scalar i, n
    
    res = did_placebo(Y, Gi, time_std, J(rows(Y), 0, .), lags, 1)
    
    n = rows(res.lags)
    if (n == 0) {
        return(J(0, 3, .))
    }
    
    output = J(n, 3, .)
    for (i = 1; i <= n; i++) {
        output[i, 1] = res.lags[i]
        output[i, 2] = res.est[i]
        output[i, 3] = res.est_std[i]
    }
    
    return(output)
}

/*---------------------------------------------------------------------------
 * _test_did_placebo_boot() - Test Wrapper for Bootstrap Functions
 *
 * Wrapper function for interactive testing of bootstrap SE computation.
 * Returns a matrix for easier inspection.
 *
 * Arguments:
 *   Y        : real colvector - outcome variable
 *   Gi       : real colvector - group indicator
 *   time_std : real colvector - standardized time index
 *   id_unit  : real colvector - unit identifier
 *   lags     : real rowvector - lag periods to test
 *   n_boot   : real scalar - number of bootstrap iterations
 *
 * Returns:
 *   real matrix (n_lags x 5): columns are [lag, est, est_std, se, se_std]
 *---------------------------------------------------------------------------*/
real matrix _test_did_placebo_boot(real colvector Y, real colvector Gi,
                                   real colvector time_std, real colvector id_unit,
                                   real rowvector lags, real scalar n_boot)
{
    struct did_data scalar data
    struct placebo_result scalar point_res
    struct placebo_boot_result scalar boot_res
    real matrix output
    real scalar i, n
    
    // Populate data structure
    data.outcome = Y
    data.Gi = Gi
    data.id_time_std = time_std
    data.id_unit = id_unit
    data.covariates = J(0, 0, .)
    data.cluster_var = J(0, 1, .)
    data.is_panel = 1
    
    // Get point estimates
    point_res = did_placebo(Y, Gi, time_std, J(rows(Y), 0, .), lags, 1)
    
    // Get bootstrap SE
    boot_res = did_placebo_boot_full(data, lags, n_boot, 1, "")
    
    n = rows(point_res.lags)
    if (n == 0) {
        return(J(0, 5, .))
    }
    
    output = J(n, 5, .)
    for (i = 1; i <= n; i++) {
        output[i, 1] = point_res.lags[i]
        output[i, 2] = point_res.est[i]
        output[i, 3] = point_res.est_std[i]
        if (i <= rows(boot_res.se)) {
            output[i, 4] = boot_res.se[i]
            output[i, 5] = boot_res.se_std[i]
        }
    }
    
    return(output)
}

/*---------------------------------------------------------------------------
 * _diddesign_check_main() - Main Entry Point for diddesign_check
 *
 * Called from diddesign_check.ado to perform parallel trends diagnostics.
 * Reads data from Stata, computes placebo estimates and bootstrap SE,
 * and stores results in external global variables.
 *
 * Arguments:
 *   depvar      : string scalar - outcome variable name
 *   treatment   : string scalar - treatment variable name
 *   id_var      : string scalar - unit identifier variable name
 *   time_var    : string scalar - time variable name
 *   post_var    : string scalar - post-treatment indicator (RCS only)
 *   covars      : string scalar - covariate names (optional)
 *   cluster_var : string scalar - cluster variable name
 *   touse       : string scalar - sample marker variable name
 *   design      : string scalar - design type ("did" or "sa")
 *   lags        : real rowvector - lag values for placebo tests
 *   n_boot      : real scalar - number of bootstrap iterations
 *   thres       : real scalar - staggered adoption threshold
 *   is_panel    : real scalar - data type indicator
 *   quiet       : real scalar - suppress progress (1=yes, 0=no)
 *
 * Side Effects:
 *   Populates external globals: _check_placebo, _check_trends, _check_Gmat,
 *   _check_n_lags, _check_n_boot_valid, _check_filtered_lags
 *---------------------------------------------------------------------------*/
void _diddesign_check_main(
    string scalar depvar,
    string scalar treatment,
    string scalar id_var,
    string scalar time_var,
    string scalar post_var,
    string scalar covars,
    string scalar cluster_var,
    string scalar touse,
    string scalar design,
    real rowvector lags,
    real scalar n_boot,
    real scalar thres,
    real scalar is_panel,
    real scalar quiet
)
{
    // Declare external global result variables
    external real matrix _check_placebo
    external real matrix _check_trends
    external real matrix _check_Gmat
    external real scalar _check_n_lags
    external real scalar _check_n_boot_valid
    external string scalar _check_filtered_lags
    
    struct did_data scalar data
    struct placebo_result scalar point_res
    struct placebo_boot_result scalar boot_res
    real matrix eq_ci, trends_data
    real colvector Y, D, id_unit, id_time, Gi, id_time_std
    real colvector cluster_col, It_post
    real matrix X
    string rowvector covar_list
    real scalar N, n_lags, i, max_lag
    real rowvector valid_lags, filtered_lags_vec
    string scalar filtered_str
    
    // Read data from Stata
    Y = st_data(., depvar, touse)
    D = st_data(., treatment, touse)
    id_time = st_data(., time_var, touse)
    N = rows(Y)
    
    // Handle panel vs RCS data
    if (is_panel) {
        id_unit = st_data(., id_var, touse)
    }
    else {
        // RCS data: use observation row numbers as pseudo-id
        id_unit = (1::N)
    }
    
    // Read post-treatment indicator for RCS data
    if (!is_panel && post_var != "") {
        It_post = st_data(., post_var, touse)
    }
    else {
        It_post = J(0, 1, .)
    }
    
    // Read covariates if specified
    if (covars != "") {
        covar_list = tokens(covars)
        X = st_data(., covar_list, touse)
    }
    else {
        X = J(N, 0, .)
    }
    
    // Read cluster variable if specified
    if (cluster_var != "") {
        cluster_col = st_data(., cluster_var, touse)
    }
    else {
        cluster_col = J(0, 1, .)
    }
    
    // Step 1.5: Normalize id_time to consecutive integers (1, 2, 3, ...)
    // This ensures lag filtering works correctly for any time scale
    {
        real colvector unique_times, id_time_norm
        real scalar n_times_uniq
        transmorphic scalar time_map
        
        unique_times = uniqrows(id_time)
        n_times_uniq = rows(unique_times)
        
        // Build hash map for O(1) lookup
        time_map = asarray_create("real", 1)
        for (i = 1; i <= n_times_uniq; i++) {
            asarray(time_map, unique_times[i], i)
        }
        
        // Map to normalized integers, preserving missing values
        id_time_norm = J(N, 1, .)
        for (i = 1; i <= N; i++) {
            if (id_time[i] >= .) {
                id_time_norm[i] = .
            }
            else {
                id_time_norm[i] = asarray(time_map, id_time[i])
            }
        }
        
        id_time = id_time_norm
    }
    
    // Step 1.6: Normalize id_unit for staggered adoption design (required for Gmat)
    if (design == "sa") {
        real colvector unique_units, id_unit_norm
        real scalar n_units_uniq
        transmorphic scalar unit_map
        
        unique_units = uniqrows(id_unit)
        n_units_uniq = rows(unique_units)
        
        unit_map = asarray_create("real", 1)
        for (i = 1; i <= n_units_uniq; i++) {
            asarray(unit_map, unique_units[i], i)
        }
        
        // Map to normalized integers, preserving missing values
        id_unit_norm = J(N, 1, .)
        for (i = 1; i <= N; i++) {
            if (id_unit[i] >= .) {
                id_unit_norm[i] = .
            }
            else {
                id_unit_norm[i] = asarray(unit_map, id_unit[i])
            }
        }
        
        id_unit = id_unit_norm
    }
    
    // Compute Gi (group indicator) and id_time_std (standardized time)
    if (is_panel) {
        // Panel: Gi = 1 if unit ever treated
        // id_time_std = time relative to treatment
        _compute_Gi_and_time_std(Y, D, id_unit, id_time, &Gi, &id_time_std)
    }
    else {
        // RCS: Gi = treatment indicator, id_time_std = time relative to treatment year
        _compute_Gi_and_time_std_rcs(D, id_time, It_post, &Gi, &id_time_std)
    }
    
    // -------------------------------------------------------------------------
    // Populate did_data structure
    // -------------------------------------------------------------------------
    data.outcome = Y
    data.treatment = D
    data.id_unit = id_unit
    data.id_time = id_time
    data.covariates = X
    data.Gi = Gi
    data.id_time_std = id_time_std
    data.N = N
    data.is_panel = is_panel
    
    if (rows(cluster_col) > 0) {
        data.cluster_var = cluster_col
    }
    else {
        data.cluster_var = J(0, 1, .)
    }
    
    // -------------------------------------------------------------------------
    // Filter lags and track filtered ones
    // -------------------------------------------------------------------------
    max_lag = abs(min(id_time_std))
    valid_lags = select(lags, lags :< max_lag)
    filtered_lags_vec = select(lags, lags :>= max_lag)
    
    // Build filtered lags string for warning
    filtered_str = ""
    if (cols(filtered_lags_vec) > 0) {
        for (i = 1; i <= cols(filtered_lags_vec); i++) {
            if (i > 1) filtered_str = filtered_str + " "
            filtered_str = filtered_str + strofreal(filtered_lags_vec[i])
        }
    }
    _check_filtered_lags = filtered_str
    
    n_lags = cols(valid_lags)
    _check_n_lags = n_lags
    
    // Handle case with no valid lags
    if (n_lags == 0) {
        _check_placebo = J(0, 7, .)
        _check_trends = J(0, 5, .)
        _check_Gmat = J(0, 0, .)
        _check_n_boot_valid = 0
        return
    }
    
    // Branch by design type
    if (design == "did") {
        _check_std_did(data, valid_lags, n_boot, cluster_var, quiet)
    }
    else if (design == "sa") {
        _check_sa_did(data, valid_lags, n_boot, thres, cluster_var, quiet)
    }
    else {
        errprintf("Error: Invalid design type '%s'. Expected 'did' or 'sa'.\n", design)
        _check_placebo = J(0, 7, .)
        _check_trends = J(0, 5, .)
        _check_Gmat = J(0, 0, .)
        _check_n_boot_valid = 0
        return
    }
    
    // -------------------------------------------------------------------------
    // Compute trends data
    // -------------------------------------------------------------------------
    _check_trends = _compute_trends(Y, Gi, id_time_std)
}

/*---------------------------------------------------------------------------
 * _compute_Gi_and_time_std() - Compute Gi and Standardized Time (Panel)
 *
 * For panel data, computes group indicator and time relative to treatment.
 *
 * Arguments:
 *   Y           : real colvector - outcome variable
 *   D           : real colvector - treatment indicator
 *   id_unit     : real colvector - unit identifier
 *   id_time     : real colvector - time identifier
 *   Gi          : pointer(real colvector) - output group indicator
 *   id_time_std : pointer(real colvector) - output standardized time
 *
 * Output:
 *   Gi = 1 if unit ever treated, 0 otherwise
 *   id_time_std: time relative to treatment (0 = treatment)
 *     - Treated units: relative to own treatment time
 *     - Control units: relative to max treatment time across treated
 *---------------------------------------------------------------------------*/
void _compute_Gi_and_time_std(
    real colvector Y,
    real colvector D,
    real colvector id_unit,
    real colvector id_time,
    pointer(real colvector) scalar Gi,
    pointer(real colvector) scalar id_time_std
)
{
    // Optimized algorithm using asarray for O(N + n_units) complexity
    real scalar N, n_units, i, u, u_idx, treat_time, max_treat_time
    real colvector units, unit_treat_time, unit_Gi
    real colvector idx
    transmorphic scalar unit_idx_map, unit_to_pos
    
    N = rows(Y)
    units = uniqrows(id_unit)
    n_units = rows(units)
    
    // Build hash map: unit -> position in 'units' array
    unit_to_pos = asarray_create("real")
    for (i = 1; i <= n_units; i++) {
        asarray(unit_to_pos, units[i], i)
    }
    
    // Build observation index lists for each unit in O(N) time
    unit_idx_map = asarray_create("real")
    for (i = 1; i <= N; i++) {
        u = id_unit[i]
        u_idx = asarray(unit_to_pos, u)
        if (asarray_contains(unit_idx_map, u_idx)) {
            asarray(unit_idx_map, u_idx, asarray(unit_idx_map, u_idx) \ i)
        }
        else {
            asarray(unit_idx_map, u_idx, i)
        }
    }
    
    // Initialize unit-level arrays
    unit_treat_time = J(n_units, 1, .)  // Treatment time for each unit (. if never treated)
    unit_Gi = J(n_units, 1, 0)          // Group indicator for each unit
    
    // Find treatment time for each unit using pre-built index lists
    for (i = 1; i <= n_units; i++) {
        idx = asarray(unit_idx_map, i)
        
        // Check if unit ever treated
        if (any(D[idx] :== 1)) {
            // Find first treatment time
            treat_time = min(select(id_time[idx], D[idx] :== 1))
            unit_treat_time[i] = treat_time
            unit_Gi[i] = 1
        }
    }
    
    // Find max treatment time among treated units
    max_treat_time = max(select(unit_treat_time, unit_Gi :== 1))
    
    // Handle case where no units are treated
    if (missing(max_treat_time)) {
        max_treat_time = max(id_time)
    }
    
    // Compute Gi and id_time_std for each observation
    // Using pre-built index lists - no additional selectindex() calls
    *Gi = J(N, 1, .)
    *id_time_std = J(N, 1, .)
    
    for (i = 1; i <= n_units; i++) {
        idx = asarray(unit_idx_map, i)
        
        // Set Gi
        (*Gi)[idx] = J(length(idx), 1, unit_Gi[i])
        
        // Set id_time_std
        if (unit_Gi[i] == 1) {
            // Treated unit: relative to own treatment time
            (*id_time_std)[idx] = id_time[idx] :- unit_treat_time[i]
        }
        else {
            // Control unit: relative to max treatment time
            (*id_time_std)[idx] = id_time[idx] :- max_treat_time
        }
    }
}

/*---------------------------------------------------------------------------
 * _compute_Gi_and_time_std_rcs() - Compute Gi and Standardized Time (RCS)
 *
 * For repeated cross-section data, computes group indicator and
 * time relative to treatment period.
 *
 * Arguments:
 *   D           : real colvector - treatment/group indicator
 *   id_time     : real colvector - time identifier
 *   It_post     : real colvector - post-treatment indicator 
 *   Gi          : pointer(real colvector) - output group indicator
 *   id_time_std : pointer(real colvector) - output standardized time
 *
 * Output:
 *   Gi = D (treatment variable is the group indicator for RCS)
 *   id_time_std = normalized_time - treat_year
 *   where treat_year = min(time where It_post == 1)
 *---------------------------------------------------------------------------*/
void _compute_Gi_and_time_std_rcs(
    real colvector D,
    real colvector id_time,
    real colvector It_post,
    pointer(real colvector) scalar Gi,
    pointer(real colvector) scalar id_time_std
)
{
    real scalar N, i, treat_year
    real colvector unique_times, id_time_norm
    real colvector post_times
    transmorphic scalar time_map
    
    N = rows(D)
    
    // Gi = D (for RCS, treatment variable IS the group indicator)
    *Gi = D
    
    // Normalize id_time to sequential integers (1, 2, 3, ...)
    unique_times = uniqrows(id_time)
    time_map = asarray_create("real", 1)
    for (i = 1; i <= rows(unique_times); i++) {
        asarray(time_map, unique_times[i], i)
    }
    
    // Map to normalized integers, preserving missing values
    id_time_norm = J(N, 1, .)
    for (i = 1; i <= N; i++) {
        if (id_time[i] >= .) {
            id_time_norm[i] = .
        }
        else {
            id_time_norm[i] = asarray(time_map, id_time[i])
        }
    }
    
    // Find treat_year = min(id_time where It_post == 1)
    post_times = select(id_time_norm, It_post :== 1)
    if (rows(post_times) > 0) {
        treat_year = min(post_times)
    }
    else {
        treat_year = max(id_time_norm)
    }
    
    // id_time_std = id_time - treat_year
    *id_time_std = id_time_norm :- treat_year
}

/*---------------------------------------------------------------------------
 * _check_std_did() - Standard DID Placebo Tests
 *
 * Computes placebo estimates and bootstrap SE for standard DID design.
 *
 * Arguments:
 *   data        : struct did_data - data structure
 *   lags        : real rowvector - lag values to test
 *   n_boot      : real scalar - number of bootstrap iterations
 *   cluster_var : string scalar - cluster variable name
 *   quiet       : real scalar - suppress progress (1=yes)
 *
 * Side Effects:
 *   Populates external globals _check_placebo, _check_n_boot_valid
 *---------------------------------------------------------------------------*/
void _check_std_did(
    struct did_data scalar data,
    real rowvector lags,
    real scalar n_boot,
    string scalar cluster_var,
    real scalar quiet
)
{
    external real matrix _check_placebo
    external real scalar _check_n_boot_valid
    
    struct placebo_result scalar point_res
    struct placebo_boot_result scalar boot_res
    real matrix eq_ci
    real scalar n_lags, i
    
    // -------------------------------------------------------------------------
    // Compute point estimates
    // -------------------------------------------------------------------------
    point_res = did_placebo(data.outcome, data.Gi, data.id_time_std, 
                            data.covariates, lags, 1)
    
    n_lags = rows(point_res.lags)
    
    // Handle empty result
    if (n_lags == 0) {
        _check_placebo = J(0, 7, .)
        _check_n_boot_valid = 0
        return
    }
    
    // -------------------------------------------------------------------------
    // Compute bootstrap standard errors
    // -------------------------------------------------------------------------
    boot_res = did_placebo_boot_full(data, lags, n_boot, data.is_panel, cluster_var)
    _check_n_boot_valid = boot_res.n_valid
    
    // Compute equivalence CIs (with dimension check)
    if (rows(point_res.est_std) != rows(boot_res.se_std)) {
        printf("{err}Warning: Dimension mismatch between point estimates (%g) and bootstrap SE (%g)\n",
               rows(point_res.est_std), rows(boot_res.se_std))
        printf("{err}Using minimum dimension for equivalence CI computation\n")
        real scalar min_dim
        min_dim = min((rows(point_res.est_std), rows(boot_res.se_std)))
        eq_ci = J(n_lags, 2, .)
        if (min_dim > 0) {
            eq_ci[1::min_dim, .] = compute_eq_ci_vec(point_res.est_std[1::min_dim], boot_res.se_std[1::min_dim])
        }
    }
    else {
        eq_ci = compute_eq_ci_vec(point_res.est_std, boot_res.se_std)
    }
    
    // -------------------------------------------------------------------------
    // Build result matrix
    // Columns: lag, estimate, std_error, estimate_orig, std_error_orig, EqCI95_LB, EqCI95_UB
    // -------------------------------------------------------------------------
    _check_placebo = J(n_lags, 7, .)
    
    for (i = 1; i <= n_lags; i++) {
        _check_placebo[i, 1] = point_res.lags[i]           // lag
        _check_placebo[i, 2] = point_res.est_std[i]        // estimate (standardized)
        _check_placebo[i, 3] = boot_res.se_std[i]          // std_error (standardized)
        _check_placebo[i, 4] = point_res.est[i]            // estimate_orig (raw)
        _check_placebo[i, 5] = boot_res.se[i]              // std_error_orig (raw)
        _check_placebo[i, 6] = eq_ci[i, 1]                 // EqCI95_LB
        _check_placebo[i, 7] = eq_ci[i, 2]                 // EqCI95_UB
    }
}

/*---------------------------------------------------------------------------
 * _check_sa_did() - Staggered Adoption Placebo Tests
 *
 * Computes placebo estimates and bootstrap SE for staggered adoption design.
 *
 * Arguments:
 *   data        : struct did_data - data structure
 *   lags        : real rowvector - lag values to test
 *   n_boot      : real scalar - number of bootstrap iterations
 *   thres       : real scalar - minimum treated fraction threshold
 *   cluster_var : string scalar - cluster variable name
 *   quiet       : real scalar - suppress progress (1=yes)
 *
 * Side Effects:
 *   Populates external globals _check_placebo, _check_Gmat, _check_n_boot_valid
 *---------------------------------------------------------------------------*/
void _check_sa_did(
    struct did_data scalar data,
    real rowvector lags,
    real scalar n_boot,
    real scalar thres,
    string scalar cluster_var,
    real scalar quiet
)
{
    external real matrix _check_placebo
    external real matrix _check_Gmat
    external real scalar _check_n_boot_valid
    
    struct did_option scalar option
    struct sa_placebo_result scalar point_res
    struct sa_placebo_boot_result scalar boot_res
    real matrix eq_ci
    real scalar n_lags, i
    
    // -------------------------------------------------------------------------
    // Setup option structure
    // -------------------------------------------------------------------------
    option = init_did_option()
    option.lag = lags
    option.n_boot = n_boot
    option.thres = thres
    option.stdz = 1
    option.quiet = quiet
    
    // Handle cluster variable (default to id_unit for panel data)
    if (cluster_var == "" & data.is_panel) {
        option.id_cluster = "id_unit"
    }
    else {
        option.id_cluster = cluster_var
    }
    
    // Compute point estimates
    point_res = did_sad_placebo(data, option)
    
    n_lags = rows(point_res.estimates)
    
    // Store Gmat only if staggered adoption succeeded (valid Gmat has >1 row and >1 column)
    if (n_lags > 0 && rows(point_res.Gmat) > 1 && cols(point_res.Gmat) > 1) {
        _check_Gmat = point_res.Gmat
    }
    else {
        // Set to empty matrix (0x0) instead of placeholder
        // This will cause ado level to skip e(Gmat) storage
        _check_Gmat = J(0, 0, .)
    }
    
    // Handle empty result
    if (n_lags == 0) {
        _check_placebo = J(0, 7, .)
        _check_n_boot_valid = 0
        return
    }
    
    // Compute bootstrap standard errors
    boot_res = did_sad_placebo_boot(data, option)
    _check_n_boot_valid = boot_res.n_valid
    
    // Compute equivalence CIs (with dimension check)
    if (rows(point_res.estimates) != rows(boot_res.se_std)) {
        printf("{err}Warning: Dimension mismatch in staggered adoption placebo results\n")
        printf("{err}Point estimates: %g rows, Bootstrap SE: %g rows\n",
               rows(point_res.estimates), rows(boot_res.se_std))
        eq_ci = J(n_lags, 2, .)
    }
    else {
        eq_ci = compute_eq_ci_vec(point_res.estimates[., 1], boot_res.se_std)
    }
    
    // Build result matrix
    // Columns: lag, estimate, std_error, estimate_orig, std_error_orig, EqCI95_LB, EqCI95_UB
    _check_placebo = J(n_lags, 7, .)
    
    for (i = 1; i <= n_lags; i++) {
        _check_placebo[i, 1] = option.lag[i]                    // lag
        _check_placebo[i, 2] = point_res.estimates[i, 1]        // estimate (standardized)
        _check_placebo[i, 3] = boot_res.se_std[i]               // std_error (standardized)
        _check_placebo[i, 4] = point_res.estimates[i, 2]        // estimate_orig (raw)
        _check_placebo[i, 5] = boot_res.se_orig[i]              // std_error_orig (raw)
        _check_placebo[i, 6] = eq_ci[i, 1]                      // EqCI95_LB
        _check_placebo[i, 7] = eq_ci[i, 2]                      // EqCI95_UB
    }
}

/*---------------------------------------------------------------------------
 * _compute_trends() - Compute Trends Data for Visualization
 *
 * Computes group-period summary statistics for parallel trends visualization.
 *
 * Arguments:
 *   Y           : real colvector - outcome variable
 *   Gi          : real colvector - group indicator
 *   id_time_std : real colvector - standardized time
 *
 * Returns:
 *   real matrix (n_rows x 5) with columns:
 *     [id_time_std, Gi, outcome_mean, outcome_sd, n_obs]
 *   Rows are ordered by (time, group) and exclude empty cells.
 *   outcome_sd is sample standard deviation (n-1 denominator)
 *---------------------------------------------------------------------------*/
real matrix _compute_trends(
    real colvector Y,
    real colvector Gi,
    real colvector id_time_std
)
{
    real matrix result
    real colvector times, groups, idx
    real scalar n_times, n_groups, t, g, row, n_obs
    real scalar y_mean, y_sd
    
    // Get unique times and groups
    times = uniqrows(id_time_std)
    groups = uniqrows(Gi)
    n_times = rows(times)
    n_groups = rows(groups)
    
    // Allocate result matrix
    result = J(n_times * n_groups, 5, .)
    
    row = 1
    for (t = 1; t <= n_times; t++) {
        for (g = 1; g <= n_groups; g++) {
            // Find observations for this time-group combination
            idx = selectindex((id_time_std :== times[t]) :& (Gi :== groups[g]))
            
            // Filter out missing outcome values before counting
            if (length(idx) > 0) {
                idx = select(idx, Y[idx] :< .)
            }
            
            n_obs = length(idx)
            
            if (n_obs > 0) {
                y_mean = mean(Y[idx])
                // SD uses n-1 denominator (sample standard deviation)
                y_sd = sqrt(variance(Y[idx]))
            }
            else {
                y_mean = .
                y_sd = .
            }
            
            result[row, 1] = times[t]      // id_time_std
            result[row, 2] = groups[g]     // Gi
            result[row, 3] = y_mean        // outcome_mean
            result[row, 4] = y_sd          // outcome_sd (SD, not SE)
            result[row, 5] = n_obs         // n_obs
            
            row++
        }
    }
    
    // Filter out rows with n_obs == 0
    real colvector valid_rows
    valid_rows = selectindex(result[., 5] :> 0)
    if (rows(valid_rows) > 0) {
        result = result[valid_rows, .]
    }
    else {
        // If all rows are empty, return empty matrix with correct dimensions
        result = J(0, 5, .)
    }
    
    return(result)
}

// ============================================================================
// GLOBAL VARIABLE INITIALIZATION
// ============================================================================
// External variables for communication between Mata functions and ado file.
// Initialized at module load time.
// ============================================================================
void _diddesign_check_init_globals()
{
    external real matrix _check_placebo
    external real matrix _check_trends
    external real matrix _check_Gmat
    external real scalar _check_n_lags
    external real scalar _check_n_boot_valid
    external string scalar _check_filtered_lags
    
    // Initialize with empty/default values
    _check_placebo = J(0, 7, .)
    _check_trends = J(0, 5, .)
    _check_Gmat = J(0, 0, .)
    _check_n_lags = 0
    _check_n_boot_valid = 0
    _check_filtered_lags = ""
}

// Call initialization function immediately when module is loaded
_diddesign_check_init_globals()

end
