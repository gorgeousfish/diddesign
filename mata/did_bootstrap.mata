*! did_bootstrap.mata - Bootstrap variance estimation for DID estimators
*!
*! Implements block bootstrap for variance-covariance estimation of the
*! standard DID and sequential DID estimators. The VCOV matrix serves as
*! the optimal weight matrix in the GMM framework for the double DID estimator.
*!
*! Main functions:
*!   did_boot_std()    - Block bootstrap for standard DID design
*!   sample_panel()    - Unit-level bootstrap for staggered adoption design
*!   bootstrap_se()    - Extract standard errors from bootstrap VCOV
*!   bootstrap_ci()    - Compute percentile confidence intervals

version 16.0

mata:
mata set matastrict on

// ============================================================================
// Block Bootstrap for Variance Estimation
// ============================================================================
// The variance-covariance matrix of (tau_DID, tau_sDID) is computed via
// cluster bootstrap:
//
//   Var(tau_DID, tau_sDID) = (1/B) * sum_b (tau^{(b)} - tau_bar)(tau^{(b)} - tau_bar)'
//
// This VCOV matrix is used as the optimal weight matrix W* in GMM estimation
// of the double DID estimator:
//
//   tau_ddid = argmin (tau - tau_DID, tau - tau_sDID)' W* (tau - tau_DID, tau - tau_sDID)
//
// Sampling is performed at the cluster level to preserve within-cluster
// correlation structure required for valid inference.
// ============================================================================

// ----------------------------------------------------------------------------
// Bootstrap Configuration Constants
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _BOOT_FAIL_WARN_PCT() - Warning threshold for bootstrap failure rate
 *---------------------------------------------------------------------------*/
real scalar _BOOT_FAIL_WARN_PCT()
{
    return(0.10)
}

/*---------------------------------------------------------------------------
 * _BOOT_MIN_SUCCESS() - Minimum successful iterations for reliable variance
 *---------------------------------------------------------------------------*/
real scalar _BOOT_MIN_SUCCESS()
{
    return(10)
}

// ----------------------------------------------------------------------------
// Data Structures
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * struct boot_result - Bootstrap result container
 *
 * Stores bootstrap estimates and metadata for variance computation.
 *---------------------------------------------------------------------------*/
struct boot_result {
    real matrix    estimates     // Bootstrap estimates (n_successful x 2*n_lead)
    pointer vector vcov          // VCOV matrix for each lead (pointer to 2x2 matrix)
    real scalar    n_successful  // Number of successful iterations
    real scalar    n_failed      // Number of failed iterations
}

// ----------------------------------------------------------------------------
// Cluster Sampling
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _boot_sample_clusters() - Sample cluster IDs with replacement
 *
 * Bootstrap sampling is performed at the cluster level to preserve
 * within-cluster correlation structure required for valid inference.
 *
 * Arguments:
 *   id_cluster_vec : real colvector - unique cluster identifiers
 *   n_clusters     : real scalar - number of clusters
 *
 * Returns:
 *   real colvector - sampled cluster IDs (same length as input)
 *---------------------------------------------------------------------------*/
real colvector _boot_sample_clusters(real colvector id_cluster_vec,
                                     real scalar n_clusters)
{
    real colvector random_idx, id_boot
    
    random_idx = safe_sample_idx(n_clusters, n_clusters)
    id_boot = id_cluster_vec[random_idx]
    
    return(id_boot)
}

// ----------------------------------------------------------------------------
// Bootstrap Dataset Construction
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _boot_create_dataset() - Construct bootstrap dataset from sampled clusters
 *
 * A new dataset is created by selecting observations belonging to sampled
 * clusters and reassigning sequential unit IDs. For panel data without an
 * explicit cluster variable, clustering is performed by id_unit.
 *
 * Arguments:
 *   data           : struct did_data - original data
 *   id_boot        : real colvector - sampled cluster IDs
 *   id_cluster_vec : real colvector - unique cluster identifiers
 *
 * Returns:
 *   struct did_data - bootstrap sample with reassigned unit IDs
 *---------------------------------------------------------------------------*/
struct did_data scalar _boot_create_dataset(struct did_data scalar data,
                                            real colvector id_boot,
                                            real colvector id_cluster_vec)
{
    struct did_data scalar dat_boot
    real scalar n_clusters, j, n_obs, start_idx, has_cluster, n_j, use_id_unit
    real colvector idx, cluster_col
    
    n_clusters = rows(id_boot)
    has_cluster = (rows(data.cluster_var) > 0 && cols(data.cluster_var) > 0)
    
    // Determine clustering column: explicit cluster, id_unit, or row-level
    use_id_unit = (!has_cluster && data.is_panel && rows(data.id_unit) > 0)
    
    if (has_cluster) {
        cluster_col = data.cluster_var
    }
    else if (use_id_unit) {
        cluster_col = data.id_unit
    }
    else {
        cluster_col = J(0, 1, .)
    }
    
    // First pass: count total observations
    n_obs = 0
    for (j = 1; j <= n_clusters; j++) {
        if (rows(cluster_col) == 0) {
            n_obs = n_obs + 1
        }
        else {
            n_obs = n_obs + sum(cluster_col :== id_boot[j])
        }
    }
    
    if (n_obs == 0) {
        dat_boot.N = 0
        return(dat_boot)
    }
    
    // Initialize bootstrap data structure
    dat_boot.outcome = J(n_obs, 1, .)
    dat_boot.outcome_delta = J(n_obs, 1, .)
    dat_boot.treatment = J(n_obs, 1, .)
    dat_boot.id_unit = J(n_obs, 1, .)
    dat_boot.id_time = J(n_obs, 1, .)
    dat_boot.id_time_std = J(n_obs, 1, .)
    dat_boot.Gi = J(n_obs, 1, .)
    dat_boot.It = J(n_obs, 1, .)
    
    if (cols(data.covariates) > 0) {
        dat_boot.covariates = J(n_obs, cols(data.covariates), .)
    }
    else {
        dat_boot.covariates = J(0, 0, .)
    }
    
    if (has_cluster) {
        dat_boot.cluster_var = J(n_obs, 1, .)
    }
    else {
        dat_boot.cluster_var = J(0, 1, .)
    }
    
    // Second pass: fill data
    start_idx = 1
    for (j = 1; j <= n_clusters; j++) {
        if (rows(cluster_col) == 0) {
            dat_boot.outcome[start_idx] = data.outcome[id_boot[j]]
            dat_boot.treatment[start_idx] = data.treatment[id_boot[j]]
            dat_boot.id_time[start_idx] = data.id_time[id_boot[j]]
            dat_boot.id_unit[start_idx] = j
            
            if (cols(data.covariates) > 0) {
                dat_boot.covariates[start_idx, .] = data.covariates[id_boot[j], .]
            }
            
            start_idx = start_idx + 1
        }
        else {
            idx = selectindex(cluster_col :== id_boot[j])
            n_j = rows(idx)
            if (n_j == 0) continue
            
            if (n_j == 1) {
                dat_boot.outcome[start_idx] = data.outcome[idx[1]]
                dat_boot.treatment[start_idx] = data.treatment[idx[1]]
                dat_boot.id_time[start_idx] = data.id_time[idx[1]]
                dat_boot.id_unit[start_idx] = j
                
                if (has_cluster) {
                    dat_boot.cluster_var[start_idx] = data.cluster_var[idx[1]]
                }
                
                if (cols(data.covariates) > 0) {
                    dat_boot.covariates[start_idx, .] = data.covariates[idx[1], .]
                }
            }
            else {
                dat_boot.outcome[start_idx::(start_idx+n_j-1)] = data.outcome[idx]
                dat_boot.treatment[start_idx::(start_idx+n_j-1)] = data.treatment[idx]
                dat_boot.id_time[start_idx::(start_idx+n_j-1)] = data.id_time[idx]
                dat_boot.id_unit[start_idx::(start_idx+n_j-1)] = J(n_j, 1, j)
                
                if (has_cluster) {
                    dat_boot.cluster_var[start_idx::(start_idx+n_j-1)] = data.cluster_var[idx]
                }
                
                if (cols(data.covariates) > 0) {
                    dat_boot.covariates[start_idx::(start_idx+n_j-1), .] = data.covariates[idx, .]
                }
            }
            
            start_idx = start_idx + n_j
        }
    }
    
    // Set metadata
    dat_boot.N = n_obs
    dat_boot.n_units = n_clusters
    dat_boot.n_periods = data.n_periods
    dat_boot.is_panel = data.is_panel
    dat_boot.treat_year = data.treat_year
    
    return(dat_boot)
}

// ----------------------------------------------------------------------------
// Outcome Delta Computation
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _compute_outcome_delta() - Outcome transformation for sequential DID
 *
 * The outcome transformation required for the sequential DID estimator
 * is computed. This estimator is consistent under the parallel trends-in-
 * trends assumption:
 *
 *   DeltaY_{it} = Y_{it} - E[Y_{i,t-1} | G_i]
 *
 * where E[Y_{i,t-1} | G_i] is the group mean outcome in the previous period.
 *
 * The sequential DID estimator subtracts the pre-treatment DID from the
 * standard DID to remove bias when trends are not parallel but the change
 * in trends is the same across groups (parallel trends-in-trends).
 *
 * Arguments:
 *   data : struct did_data - data with outcome, Gi, id_time_std
 *
 * Returns:
 *   real colvector - transformed outcome (Y - lagged group mean)
 *---------------------------------------------------------------------------*/
real colvector _compute_outcome_delta(struct did_data scalar data)
{
    real scalar n, i, g, t
    real colvector outcome_delta, Ymean
    string scalar key, lag_key
    transmorphic scalar group_sum_map, group_count_map, group_has_na_map
    real scalar sum_y, count_y
    
    n = data.N
    outcome_delta = J(n, 1, .)
    Ymean = J(n, 1, .)
    
    // Hash tables for O(N) complexity
    group_sum_map = asarray_create("string", 1)
    group_count_map = asarray_create("string", 1)
    group_has_na_map = asarray_create("string", 1)
    
    // First pass: accumulate sums and counts for each (Gi, id_time_std) group
    for (i = 1; i <= n; i++) {
        g = data.Gi[i]
        t = data.id_time_std[i]
        
        if (missing(g) || missing(t)) continue
        
        key = strofreal(g) + "_" + strofreal(t)
        
        if (missing(data.outcome[i])) {
            asarray(group_has_na_map, key, 1)
            if (!asarray_contains(group_sum_map, key)) {
                asarray(group_sum_map, key, 0)
                asarray(group_count_map, key, 0)
            }
            continue
        }
        
        if (asarray_contains(group_sum_map, key)) {
            asarray(group_sum_map, key, asarray(group_sum_map, key) + data.outcome[i])
            asarray(group_count_map, key, asarray(group_count_map, key) + 1)
        }
        else {
            asarray(group_sum_map, key, data.outcome[i])
            asarray(group_count_map, key, 1)
        }
    }
    
    // Second pass: assign lag group means using O(1) lookups
    for (i = 1; i <= n; i++) {
        g = data.Gi[i]
        t = data.id_time_std[i]
        
        if (missing(g) || missing(t)) continue
        
        lag_key = strofreal(g) + "_" + strofreal(t - 1)
        
        if (asarray_contains(group_sum_map, lag_key)) {
            if (asarray_contains(group_has_na_map, lag_key)) {
                Ymean[i] = .
            }
            else {
                sum_y = asarray(group_sum_map, lag_key)
                count_y = asarray(group_count_map, lag_key)
                if (count_y > 0) {
                    Ymean[i] = sum_y / count_y
                }
            }
        }
    }
    
    outcome_delta = data.outcome - Ymean
    
    return(outcome_delta)
}

// ----------------------------------------------------------------------------
// Bootstrap Data Preparation
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _normalize_time() - Normalize time index to consecutive integers
 *
 * Time values are compressed to consecutive integers 1, 2, 3, ..., preserving
 * the temporal ordering. This is required for bootstrap samples that may have
 * non-consecutive time periods due to resampling.
 *
 * Arguments:
 *   time_values : real colvector - original time values
 *
 * Returns:
 *   real colvector - normalized time values (1, 2, 3, ...)
 *---------------------------------------------------------------------------*/
real colvector _normalize_time(real colvector time_values)
{
    real colvector unique_times, result
    real scalar n, n_unique, i
    transmorphic scalar time_map
    
    n = rows(time_values)
    if (n == 0) return(J(0, 1, .))
    
    unique_times = uniqrows(time_values)
    n_unique = rows(unique_times)
    
    result = J(n, 1, .)
    
    // Hash map for O(1) lookup
    time_map = asarray_create("real", 1)
    for (i = 1; i <= n_unique; i++) {
        asarray(time_map, unique_times[i], i)
    }
    
    for (i = 1; i <= n; i++) {
        if (!missing(time_values[i])) {
            result[i] = asarray(time_map, time_values[i])
        }
    }
    
    return(result)
}

/*---------------------------------------------------------------------------
 * _boot_panel_prep() - Re-run panel data preparation for bootstrap sample
 *
 * Derived variables (Gi, It, id_time_std, outcome_delta) are recomputed
 * for a bootstrap sample. Time indices are re-normalized to consecutive
 * integers to handle gaps from bootstrap sampling.
 *
 * DID identification requires both treated and control units, and
 * observations in both pre- and post-treatment periods. An empty Gi
 * vector signals identification failure.
 *
 * Arguments:
 *   dat_boot : struct did_data - bootstrap sample
 *
 * Returns:
 *   struct did_data - bootstrap sample with recomputed derived variables
 *---------------------------------------------------------------------------*/
struct did_data scalar _boot_panel_prep(struct did_data scalar dat_boot)
{
    real scalar n, i, treat_year, max_treat
    real colvector unit_ids, unit_idx, treated_times
    
    n = dat_boot.N
    if (n == 0) return(dat_boot)
    
    // Compute Gi = max(treatment) by id_unit
    unit_ids = uniqrows(dat_boot.id_unit)
    dat_boot.Gi = J(n, 1, .)
    
    for (i = 1; i <= rows(unit_ids); i++) {
        unit_idx = selectindex(dat_boot.id_unit :== unit_ids[i])
        max_treat = max(dat_boot.treatment[unit_idx])
        dat_boot.Gi[unit_idx] = J(rows(unit_idx), 1, max_treat)
    }
    
    // Re-normalize time index to consecutive integers
    dat_boot.id_time = _normalize_time(dat_boot.id_time)
    
    // Identify treatment time
    treated_times = select(dat_boot.id_time, dat_boot.treatment :== 1)
    
    if (rows(treated_times) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    
    treat_year = min(treated_times)
    dat_boot.treat_year = treat_year
    
    // Validate identification: require both treated and control units
    if (sum(dat_boot.Gi :== 0) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    if (sum(dat_boot.Gi :== 1) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    
    // Compute standardized time index
    dat_boot.id_time_std = dat_boot.id_time :- treat_year
    
    // Compute post-treatment indicator
    dat_boot.It = (dat_boot.id_time :>= treat_year)
    
    // Validate identification: require both pre and post periods
    if (sum(dat_boot.It :== 0) == 0 || sum(dat_boot.It :== 1) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    
    // Compute outcome_delta for sequential DID
    dat_boot.outcome_delta = _compute_outcome_delta(dat_boot)
    
    return(dat_boot)
}

/*---------------------------------------------------------------------------
 * _boot_rcs_prep() - Re-run RCS data preparation for bootstrap sample
 *
 * For repeated cross-sectional data, Gi and It are input variables (not
 * computed from treatment timing). This function copies Gi and It from
 * original data, re-normalizes time indices, and recomputes id_time_std
 * and outcome_delta for the sequential DID estimator.
 *
 * Arguments:
 *   dat_boot  : struct did_data - bootstrap sample
 *   data_orig : struct did_data - original data (source of Gi, It)
 *   id_boot   : real colvector - sampled indices
 *
 * Returns:
 *   struct did_data - bootstrap sample with recomputed time indices
 *---------------------------------------------------------------------------*/
struct did_data scalar _boot_rcs_prep(struct did_data scalar dat_boot,
                                      struct did_data scalar data_orig,
                                      real colvector id_boot)
{
    real scalar n, j, start_idx, n_j, has_cluster, treat_year
    real scalar n_orig, boot_idx
    real colvector idx, treated_times
    
    n = dat_boot.N
    if (n == 0) return(dat_boot)
    
    has_cluster = (rows(data_orig.cluster_var) > 0 && cols(data_orig.cluster_var) > 0)
    n_orig = rows(data_orig.Gi)
    
    // Copy Gi and It from original data
    dat_boot.Gi = J(n, 1, .)
    dat_boot.It = J(n, 1, .)
    
    start_idx = 1
    for (j = 1; j <= rows(id_boot); j++) {
        if (!has_cluster) {
            boot_idx = id_boot[j]
            if (boot_idx < 1 || boot_idx > n_orig) {
                continue
            }
            if (start_idx > n) {
                break
            }
            
            dat_boot.Gi[start_idx] = data_orig.Gi[boot_idx]
            dat_boot.It[start_idx] = data_orig.It[boot_idx]
            
            start_idx = start_idx + 1
        }
        else {
            idx = selectindex(data_orig.cluster_var :== id_boot[j])
            n_j = rows(idx)
            if (n_j == 0) continue
            
            if (start_idx + n_j - 1 > n) {
                n_j = n - start_idx + 1
                if (n_j <= 0) break
            }
            
            if (n_j == 1) {
                dat_boot.Gi[start_idx] = data_orig.Gi[idx[1]]
                dat_boot.It[start_idx] = data_orig.It[idx[1]]
            }
            else {
                dat_boot.Gi[start_idx::(start_idx+n_j-1)] = data_orig.Gi[idx]
                dat_boot.It[start_idx::(start_idx+n_j-1)] = data_orig.It[idx]
            }
            
            start_idx = start_idx + n_j
        }
    }
    
    // Verify index consistency
    if (start_idx - 1 != n) {
        errprintf("{txt}Warning: _boot_rcs_prep index mismatch: expected %g rows, filled %g\n", 
                  n, start_idx - 1)
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    
    // Re-normalize time index to consecutive integers
    dat_boot.id_time = _normalize_time(dat_boot.id_time)
    
    // Recompute treat_year from normalized time
    treated_times = select(dat_boot.id_time, dat_boot.It :== 1)
    
    if (rows(treated_times) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    treat_year = min(treated_times)
    dat_boot.treat_year = treat_year
    
    // Validate identification
    if (sum(dat_boot.Gi :== 0) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    if (sum(dat_boot.Gi :== 1) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    if (sum(dat_boot.It :== 0) == 0) {
        dat_boot.Gi = J(0, 1, .)
        return(dat_boot)
    }
    
    // Recompute id_time_std and outcome_delta
    dat_boot.id_time_std = dat_boot.id_time :- treat_year
    dat_boot.outcome_delta = _compute_outcome_delta(dat_boot)
    
    return(dat_boot)
}

// ----------------------------------------------------------------------------
// Main Bootstrap Function
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * did_boot_std() - Block bootstrap for standard DID design
 *
 * The variance-covariance matrix Sigma of (tau_DID, tau_sDID) is estimated
 * via cluster bootstrap. This VCOV matrix is used to construct the optimal
 * weight matrix in the GMM framework for the double DID estimator:
 *
 *   tau_dDID = argmin (tau - tau_DID, tau - tau_sDID)' W (tau - tau_DID, tau - tau_sDID)
 *
 * where W = Sigma^{-1} is the precision matrix (optimal GMM weight).
 * The double DID estimator combines the standard DID and sequential DID
 * to achieve efficiency under parallel trends and robustness under the
 * weaker parallel trends-in-trends assumption.
 *
 * Sampling is performed at the cluster level with replacement to preserve
 * within-cluster correlation structure for valid inference.
 *
 * Arguments:
 *   data   : struct did_data - prepared data
 *   lead   : real rowvector - post-treatment period indices
 *   n_boot : real scalar - number of bootstrap iterations
 *   seed   : real scalar - random seed (optional)
 *
 * Returns:
 *   struct boot_result - bootstrap estimates and VCOV matrices
 *---------------------------------------------------------------------------*/
struct boot_result scalar did_boot_std(struct did_data scalar data,
                                       real rowvector lead,
                                       real scalar n_boot,
                                       | real scalar seed)
{
    struct boot_result scalar result
    struct did_data scalar dat_boot
    real colvector id_cluster_vec, id_boot, valid_idx
    real scalar n_clusters, n_lead, b, l, n_successful, n_failed
    real matrix boot_est
    real rowvector est_b
    
    if (args() >= 4 && !missing(seed)) {
        rseed(seed)
    }
    
    n_lead = cols(lead)
    
    // Handle n_boot=0 case
    if (n_boot == 0) {
        printf("{txt}Note: n_boot=0, no bootstrap inference available\n")
        result.estimates = J(0, 2 * n_lead, .)
        result.vcov = J(n_lead, 1, NULL)
        result.n_successful = 0
        result.n_failed = 0
        return(result)
    }
    
    // Determine cluster structure for bootstrap sampling
    if (rows(data.cluster_var) > 0 && cols(data.cluster_var) > 0) {
        real colvector valid_cluster_mask, valid_clusters
        real scalar n_missing_clusters
        
        valid_cluster_mask = (data.cluster_var :< .)
        valid_clusters = select(data.cluster_var, valid_cluster_mask)
        n_missing_clusters = rows(data.cluster_var) - rows(valid_clusters)
        
        if (n_missing_clusters > 0) {
            printf("{txt}Warning: cluster variable contains %g missing values (excluded from bootstrap)\n",
                   n_missing_clusters)
        }
        
        if (rows(valid_clusters) > 0) {
            id_cluster_vec = uniqrows(valid_clusters)
        }
        else {
            errprintf("Error: cluster variable contains only missing values\n")
            result.estimates = J(0, 2 * n_lead, .)
            result.vcov = J(n_lead, 1, NULL)
            result.n_successful = 0
            result.n_failed = n_boot
            return(result)
        }
    }
    else if (data.is_panel && rows(data.id_unit) > 0) {
        id_cluster_vec = uniqrows(data.id_unit)
    }
    else {
        id_cluster_vec = (1::data.N)
    }
    n_clusters = rows(id_cluster_vec)
    
    // Initialize bootstrap estimates matrix
    boot_est = J(n_boot, 2 * n_lead, .)
    valid_idx = J(n_boot, 1, 0)
    
    // Bootstrap loop
    for (b = 1; b <= n_boot; b++) {
        
        // Sample clusters with replacement
        id_boot = _boot_sample_clusters(id_cluster_vec, n_clusters)
        
        // Construct bootstrap dataset
        dat_boot = _boot_create_dataset(data, id_boot, id_cluster_vec)
        
        if (dat_boot.N == 0) {
            continue
        }
        
        // Re-run data preparation
        if (data.is_panel) {
            dat_boot = _boot_panel_prep(dat_boot)
        }
        else {
            dat_boot = _boot_rcs_prep(dat_boot, data, id_boot)
        }
        
        if (rows(dat_boot.Gi) == 0 || missing(dat_boot.Gi[1])) {
            continue
        }
        
        // Compute estimates for each lead
        est_b = J(1, 2 * n_lead, .)
        for (l = 1; l <= n_lead; l++) {
            est_b[1, (2*l-1)..(2*l)] = did_fit(
                dat_boot.outcome,
                dat_boot.outcome_delta,
                dat_boot.Gi,
                dat_boot.It,
                dat_boot.covariates,
                dat_boot.id_time_std,
                lead[l]
            )
        }
        
        // Check if at least one lead has valid estimates
        real scalar has_valid_est, ll_check
        has_valid_est = 0
        for (ll_check = 1; ll_check <= n_lead; ll_check++) {
            if (!missing(est_b[1, 2*ll_check-1]) || !missing(est_b[1, 2*ll_check])) {
                has_valid_est = 1
                break
            }
        }
        if (has_valid_est) {
            boot_est[b, .] = est_b
            valid_idx[b] = 1
        }
    }
    
    // Summarize results
    n_successful = sum(valid_idx)
    n_failed = n_boot - n_successful
    
    if (n_successful > 0) {
        boot_est = select(boot_est, valid_idx)
    }
    else {
        boot_est = J(0, 2 * n_lead, .)
    }
    
    // Warn if many iterations failed
    if (n_failed > _BOOT_FAIL_WARN_PCT() * n_boot) {
        real scalar pct_failed
        pct_failed = 100 * n_failed / n_boot
        errprintf("Warning: %g / %g bootstrap iterations failed (%g percent)\n",
                  n_failed, n_boot, pct_failed)
    }
    
    if (n_successful < _BOOT_MIN_SUCCESS() && n_successful > 0) {
        errprintf("Warning: Only %g bootstrap iterations succeeded, results may be unreliable\n",
                  n_successful)
    }
    
    // Compute VCOV for each lead
    result.estimates = boot_est
    result.vcov = J(n_lead, 1, NULL)
    result.n_successful = n_successful
    result.n_failed = n_failed
    
    if (n_successful >= 2) {
        for (l = 1; l <= n_lead; l++) {
            real colvector valid_lead_idx
            real matrix boot_est_lead
            real scalar n_valid_lead
            
            valid_lead_idx = selectindex(rowmissing(boot_est[., (2*l-1)..(2*l)]) :== 0)
            n_valid_lead = rows(valid_lead_idx)
            
            if (n_valid_lead >= 2) {
                boot_est_lead = boot_est[valid_lead_idx, (2*l-1)..(2*l)]
                result.vcov[l] = &(compute_vcov(boot_est_lead))
            }
            else {
                result.vcov[l] = &(J(2, 2, .))
            }
        }
    }
    
    return(result)
}

// ----------------------------------------------------------------------------
// Standard Error and Confidence Interval Functions
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * bootstrap_se() - Extract standard errors from bootstrap VCOV
 *
 * Standard errors are computed as square roots of diagonal elements:
 *   SE(tau_DID)  = sqrt(Var[1,1])
 *   SE(tau_sDID) = sqrt(Var[2,2])
 *
 * Arguments:
 *   vcov : real matrix (2 x 2) - variance-covariance matrix
 *
 * Returns:
 *   real rowvector (1 x 2) - (SE_DID, SE_sDID)
 *---------------------------------------------------------------------------*/
real rowvector bootstrap_se(real matrix vcov)
{
    real rowvector se
    
    if (rows(vcov) != 2 || cols(vcov) != 2) {
        return((., .))
    }
    
    if (missing(vcov[1,1]) || missing(vcov[2,2])) {
        return((., .))
    }
    
    se = (sqrt(vcov[1,1]), sqrt(vcov[2,2]))
    
    return(se)
}

/*---------------------------------------------------------------------------
 * bootstrap_ci() - Compute percentile bootstrap confidence intervals
 *
 * Confidence intervals are computed using the percentile method:
 *   ci_low  = quantile(boot_est, alpha/2)
 *   ci_high = quantile(boot_est, 1 - alpha/2)
 * where alpha = 1 - level/100.
 *
 * Arguments:
 *   boot_est : real colvector - bootstrap estimates for one parameter
 *   level    : real scalar - confidence level (e.g., 95)
 *
 * Returns:
 *   real rowvector (1 x 2) - (ci_low, ci_high)
 *---------------------------------------------------------------------------*/
real rowvector bootstrap_ci(real colvector boot_est, real scalar level)
{
    real scalar alpha, p_low, p_high
    real rowvector ci
    
    if (level <= 0 || level >= 100) {
        return((., .))
    }
    
    alpha = 1 - level / 100
    p_low = alpha / 2
    p_high = 1 - alpha / 2
    
    ci = (quantile_sorted(boot_est, p_low), quantile_sorted(boot_est, p_high))
    
    return(ci)
}

// ============================================================================
// Panel Bootstrap for Staggered Adoption Design
// ============================================================================
// Unit-level bootstrap is implemented for the staggered adoption (SA) design
// where treatment timing varies across units. Unlike the standard DID bootstrap
// which samples clusters, the SA design requires sampling entire units with
// replacement to preserve the time series structure needed for:
//
//   1. Reconstructing the treatment timing matrix (Gmat) encoding adoption times
//   2. Computing period-specific SA-ATT estimates: tau_DID(t), tau_sDID(t)
//   3. Aggregating via time weights: tau_bar = sum_t pi_t * tau(t)
//
// The SA double DID extends the basic double DID framework by applying the
// GMM combination of DID and sequential DID estimators at each treatment time,
// then aggregating across time periods using appropriate weights.
// ============================================================================

/*---------------------------------------------------------------------------
 * sample_panel() - Unit-level bootstrap for staggered adoption design
 *
 * Entire units are sampled with replacement, preserving the time series
 * structure required for period-specific SA-ATT estimation. Each sampled
 * unit is assigned a new sequential ID to handle duplicate units.
 *
 * Derived fields (Gi, It, id_time_std, outcome_delta) are NOT populated here.
 * These must be recomputed downstream after the treatment timing matrix (Gmat)
 * is rebuilt from the bootstrap sample.
 *
 * Arguments:
 *   data : struct did_data - panel data structure
 *
 * Returns:
 *   struct did_data - bootstrap sample with reassigned unit IDs
 *---------------------------------------------------------------------------*/
struct did_data scalar sample_panel(struct did_data scalar data)
{
    struct did_data scalar boot_data
    real colvector id_vec, id_boot, idx, new_id_unit
    real matrix boot_outcome, boot_treatment, boot_time, boot_covariates
    real colvector boot_cluster
    real scalar n_units, i, j, n_obs_total, n_obs_i
    real scalar has_covariates, has_cluster
    
    // Validate panel data structure
    if (!data.is_panel) {
        errprintf("sample_panel(): Requires panel data (is_panel=1)\n")
        errprintf("               For RCS data, use cluster bootstrap instead\n")
        boot_data.outcome = J(0, 1, .)
        boot_data.treatment = J(0, 1, .)
        boot_data.id_unit = J(0, 1, .)
        boot_data.id_time = J(0, 1, .)
        boot_data.covariates = J(0, 0, .)
        boot_data.cluster_var = J(0, 1, .)
        boot_data.Gi = J(0, 1, .)
        boot_data.It = J(0, 1, .)
        boot_data.id_time_std = J(0, 1, .)
        boot_data.outcome_delta = J(0, 1, .)
        boot_data.N = 0
        boot_data.n_units = 0
        boot_data.n_periods = 0
        boot_data.treat_year = .
        boot_data.is_panel = 0
        return(boot_data)
    }
    
    // Handle empty dataset
    if (data.N == 0 || data.n_units == 0) {
        boot_data.outcome = J(0, 1, .)
        boot_data.treatment = J(0, 1, .)
        boot_data.id_unit = J(0, 1, .)
        boot_data.id_time = J(0, 1, .)
        boot_data.covariates = J(0, 0, .)
        boot_data.cluster_var = J(0, 1, .)
        boot_data.Gi = J(0, 1, .)
        boot_data.It = J(0, 1, .)
        boot_data.id_time_std = J(0, 1, .)
        boot_data.outcome_delta = J(0, 1, .)
        boot_data.N = 0
        boot_data.n_units = 0
        boot_data.n_periods = data.n_periods
        boot_data.treat_year = .
        boot_data.is_panel = data.is_panel
        return(boot_data)
    }
    
    // Sample unit IDs with replacement
    id_vec = uniqrows(data.id_unit)
    n_units = rows(id_vec)
    id_boot = id_vec[safe_sample_idx(n_units, n_units)]
    
    has_covariates = (cols(data.covariates) > 0 && rows(data.covariates) == data.N)
    has_cluster = (rows(data.cluster_var) == data.N && cols(data.cluster_var) > 0)
    
    // Cache selectindex() results for efficiency
    transmorphic idx_cache
    idx_cache = asarray_create("real", 1)
    asarray_notfound(idx_cache, J(0, 1, .))
    
    // First pass: count total observations
    n_obs_total = 0
    for (i = 1; i <= n_units; i++) {
        idx = selectindex(data.id_unit :== id_boot[i])
        asarray(idx_cache, i, idx)
        n_obs_total = n_obs_total + rows(idx)
    }
    
    // Allocate result matrices
    boot_outcome = J(n_obs_total, 1, .)
    boot_treatment = J(n_obs_total, 1, .)
    boot_time = J(n_obs_total, 1, .)
    new_id_unit = J(n_obs_total, 1, .)
    
    if (has_covariates) {
        boot_covariates = J(n_obs_total, cols(data.covariates), .)
    }
    else {
        boot_covariates = J(0, 0, .)
    }
    
    if (has_cluster) {
        boot_cluster = J(n_obs_total, 1, .)
    }
    else {
        boot_cluster = J(0, 1, .)
    }
    
    // Second pass: fill data using cached indices
    j = 1
    for (i = 1; i <= n_units; i++) {
        idx = asarray(idx_cache, i)
        n_obs_i = rows(idx)
        
        if (n_obs_i == 0) continue
        
        if (n_obs_i == 1) {
            boot_outcome[j] = data.outcome[idx[1]]
            boot_treatment[j] = data.treatment[idx[1]]
            boot_time[j] = data.id_time[idx[1]]
            new_id_unit[j] = i
            
            if (has_covariates) {
                boot_covariates[j, .] = data.covariates[idx[1], .]
            }
            
            if (has_cluster) {
                boot_cluster[j] = data.cluster_var[idx[1]]
            }
        }
        else {
            boot_outcome[j::(j+n_obs_i-1)] = data.outcome[idx]
            boot_treatment[j::(j+n_obs_i-1)] = data.treatment[idx]
            boot_time[j::(j+n_obs_i-1)] = data.id_time[idx]
            new_id_unit[j::(j+n_obs_i-1)] = J(n_obs_i, 1, i)
            
            if (has_covariates) {
                boot_covariates[j::(j+n_obs_i-1), .] = data.covariates[idx, .]
            }
            
            if (has_cluster) {
                boot_cluster[j::(j+n_obs_i-1)] = data.cluster_var[idx]
            }
        }
        
        j = j + n_obs_i
    }
    
    // Populate result structure
    boot_data.outcome = boot_outcome
    boot_data.treatment = boot_treatment
    boot_data.id_time = boot_time
    boot_data.id_unit = new_id_unit
    boot_data.covariates = boot_covariates
    boot_data.cluster_var = boot_cluster
    
    // Update metadata
    boot_data.N = n_obs_total
    boot_data.n_units = n_units
    boot_data.n_periods = data.n_periods
    boot_data.is_panel = data.is_panel
    
    // Derived fields are NOT populated - must be computed downstream
    boot_data.Gi = J(0, 1, .)
    boot_data.It = J(0, 1, .)
    boot_data.id_time_std = J(0, 1, .)
    boot_data.outcome_delta = J(0, 1, .)
    boot_data.treat_year = .
    
    return(boot_data)
}

// ----------------------------------------------------------------------------
// Module Verification Function
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _did_bootstrap_loaded() - Verify module is loaded
 *---------------------------------------------------------------------------*/
void _did_bootstrap_loaded()
{
    printf("{txt}did_bootstrap.mata loaded successfully\n")
}

end
