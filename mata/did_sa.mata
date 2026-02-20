*! did_sa.mata - Staggered adoption design estimation
*!
*! Extends the double DID framework to settings where treatment timing varies
*! across units. Period-specific estimates are aggregated via time-weighted
*! averaging, with variance computed through panel bootstrap.

version 16.0

mata:
mata set matastrict on

// ============================================================================
// STAGGERED ADOPTION DESIGN ESTIMATION
// ============================================================================
//
// The staggered adoption (SA) design allows different units to receive treatment
// at different time periods. The SA-ATT at time t is defined as:
//   tau^SA(t) = E[Y_it(1) - Y_it(0) | G_it = 1]
//
// The time-average SA-ATT aggregates period-specific effects:
//   tau_bar^SA = sum_t pi_t * tau^SA(t)
// where pi_t = n_{1t} / sum_t' n_{1t'} (proportion treated at time t)
//
// Algorithm:
//   1. Construct treatment timing matrix G_it in {-1, 0, 1}
//   2. Identify valid periods with n_treated >= threshold
//   3. Compute time weights pi_t proportional to treated units
//   4. Estimate period-specific tau_DID(t) and tau_sDID(t) using {t-2, t-1, t}
//   5. Aggregate via time-weighted average
//   6. Compute variance via panel bootstrap
//
// ============================================================================


// ----------------------------------------------------------------------------
// DATA STRUCTURES
// ----------------------------------------------------------------------------

// Note: struct sa_point is defined in did_gmm.mata to resolve dependency order.

/*---------------------------------------------------------------------------
 * struct sa_data - Staggered Adoption Context
 *
 * Contains treatment timing information computed from panel data. This
 * structure enables code reuse between estimation and placebo tests.
 *---------------------------------------------------------------------------*/
struct sa_data {
    real matrix    Gmat          // Treatment timing matrix (N_units x T)
                                 // G_it: -1 = previously treated, 0 = control, 1 = newly treated
    real colvector id_time_use   // Valid period indices where n_treated >= threshold
    pointer vector id_subj_use   // Valid unit indices per period
    real colvector time_weight   // Time weights pi_t, normalized to sum to 1
}

/*---------------------------------------------------------------------------
 * struct sa_placebo_result - Placebo Test Results for SA Design
 *
 * Contains time-weighted placebo estimates for assessing the parallel
 * trends assumption in staggered adoption settings.
 *---------------------------------------------------------------------------*/
struct sa_placebo_result {
    real matrix    estimates     // n_lags x 2: (standardized, original scale)
    real matrix    Gmat          // Treatment pattern matrix for visualization
    real rowvector valid_lags    // Lag values included in estimation
}

/*---------------------------------------------------------------------------
 * struct sa_placebo_boot_result - Bootstrap Results for SA Placebo Tests
 *
 * Contains bootstrap inference results including standard errors and
 * bootstrap distributions for SA placebo tests.
 *---------------------------------------------------------------------------*/
struct sa_placebo_boot_result {
    real colvector se_std        // Standard errors for standardized estimates
    real colvector se_orig       // Standard errors for original-scale estimates
    real scalar    n_boot        // Number of bootstrap iterations requested
    real scalar    n_valid       // Number of successful bootstrap iterations
    real matrix    boot_est_std  // Bootstrap estimates (standardized): n_valid x n_lags
    real matrix    boot_est_orig // Bootstrap estimates (original): n_valid x n_lags
}


// ----------------------------------------------------------------------------
// SA POINT ESTIMATION FUNCTIONS
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * sa_double_did() - SA Point Estimation Coordinator
 *
 * Coordinates the staggered adoption estimation process. Treatment timing
 * structures are constructed and period-specific estimation is delegated
 * to sa_compute_did().
 *
 * Arguments:
 *   data   : struct did_data - panel data structure
 *   option : struct did_option - estimation options (thres, lead)
 *
 * Returns:
 *   struct sa_point containing time-weighted (tau_DID, tau_sDID) for each lead
 *
 * Algorithm:
 *   1. Gmat is constructed from treatment timing
 *   2. Valid periods are identified where n_treated >= threshold
 *   3. Valid units per period are selected (not previously treated)
 *   4. Time weights are computed proportional to newly treated units
 *   5. Period-specific estimation is delegated to sa_compute_did()
 *---------------------------------------------------------------------------*/
struct sa_point scalar sa_double_did(struct did_data scalar data,
                                      struct did_option scalar option)
{
    struct sa_point scalar result
    real matrix Gmat
    real colvector id_time_use, time_weight
    pointer vector id_subj_use
    real scalar n_lead
    
    // Initialize result with missing values
    n_lead = length(option.lead)
    result.DID = J(1, n_lead, .)
    result.sDID = J(1, n_lead, .)
    
    // Step 1: Create group matrix
    Gmat = create_gmat(data.id_unit, data.id_time, data.treatment)
    
    // Handle empty Gmat
    if (rows(Gmat) == 0 || cols(Gmat) == 0) {
        errprintf("sa_double_did(): Failed to create Gmat\n")
        return(result)
    }
    
    // Step 2: Get valid periods
    id_time_use = get_periods(Gmat, option.thres)
    
    // Handle no valid periods
    if (rows(id_time_use) == 0) {
        errprintf("sa_double_did(): No valid periods found\n")
        return(result)
    }
    
    // Step 3: Get valid subjects for each period
    id_subj_use = get_subjects(Gmat, id_time_use)
    
    // Step 4: Get time weights
    time_weight = get_time_weight(Gmat, id_time_use)
    
    // Verify weights sum to 1.0
    if (abs(sum(time_weight) - 1.0) > 1e-10) {
        errprintf("sa_double_did(): Time weights do not sum to 1.0\n")
        return(result)
    }
    
    // Step 5: Compute period-specific estimates and time-weighted average
    result = sa_compute_did(data, id_time_use, id_subj_use, time_weight, option.lead)
    
    return(result)
}


/*---------------------------------------------------------------------------
 * sa_compute_did() - Period-Specific DID/sDID Estimation
 *
 * Period-specific estimates are computed and aggregated via time-weighted
 * average. This function implements the core SA estimation algorithm.
 *
 * Arguments:
 *   data        : struct did_data - panel data
 *   id_time_use : real colvector - valid period indices
 *   id_subj_use : pointer vector - valid subject indices for each period
 *   time_weight : real colvector - time weights pi_t
 *   lead        : real rowvector - lead parameters for dynamic effects
 *   min_time    : real scalar - minimum time requirement (default=3)
 *
 * Returns:
 *   struct sa_point containing time-weighted (tau_DID, tau_sDID) for each lead
 *
 * Algorithm:
 *   For each valid period t where (t >= min_time) and (t + max_lead <= T):
 *     1. Data is subset to valid units and times {t-2, t-1, t}
 *     2. Future outcome observations are added if lead > 0
 *     3. Derived variables are prepared for DID estimation
 *     4. tau_DID(t) and tau_sDID(t) are computed for each lead
 *     5. Corresponding time weight is stored
 *   
 *   For each lead, time-weighted average is computed:
 *     - Missing estimates are excluded (separately for DID and sDID)
 *     - Weights are renormalized: w_norm = w / sum(w)
 *     - tau_DID = sum(w_norm * tau_DID(t))
 *     - tau_sDID = sum(w_norm * tau_sDID(t))
 *---------------------------------------------------------------------------*/
struct sa_point scalar sa_compute_did(struct did_data scalar data,
                                       real colvector id_time_use,
                                       pointer vector id_subj_use,
                                       real colvector time_weight,
                                       real rowvector lead,
                                       | real scalar min_time)
{
    struct sa_point scalar result
    struct did_data scalar dat_use, dat_did
    real matrix est_did, est_sdid
    real colvector time_weight_new, idx_subj
    real colvector valid_idx_did, valid_idx_sdid, weight_norm_did, weight_norm_sdid
    real rowvector est_t
    real scalar n_periods, max_time, i, t, iter, ll, n_lead
    real scalar n_valid, max_lead
    real scalar sum_w_did, sum_w_sdid  // For weight sum validation
    
    // Default min_time = 3 (requires t-2, t-1, t periods for sDID calculation)
    if (args() < 6) min_time = 3
    
    n_periods = rows(id_time_use)
    n_lead = length(lead)
    max_time = max(data.id_time)
    max_lead = max(lead)
    
    // Initialize result with missing values
    result.DID = J(1, n_lead, .)
    result.sDID = J(1, n_lead, .)
    
    // Handle empty input
    if (n_periods == 0) {
        errprintf("sa_compute_did(): No valid periods for SA estimation\n")
        return(result)
    }
    
    // Pre-allocate for maximum possible valid periods
    est_did = J(n_periods, n_lead, .)
    est_sdid = J(n_periods, n_lead, .)
    time_weight_new = J(n_periods, 1, .)
    
    // Initialize iter counter before the loop
    iter = 1
    
    // Loop over all periods
    for (i = 1; i <= n_periods; i++) {
        t = id_time_use[i]
        
        // Key condition checks: need min_time periods before and max_lead periods after
        if ((t >= min_time) && (t + max_lead <= max_time)) {
            
            // 1. Subset data: subjects in id_subj_use[i], times {t-2, t-1, t}
            idx_subj = *id_subj_use[i]
            
            // Skip if no valid subjects for this period
            if (rows(idx_subj) == 0) {
                continue
            }
            
            dat_use = subset_data_sa(data, idx_subj, t)
            
            // 2. Handle lead > 0 (add future outcome observations)
            if (max_lead > 0) {
                dat_use = add_lead_outcomes(dat_use, data, idx_subj, t, lead)
            }
            
            // 3. Prepare data for DID estimation
            dat_did = sa_prepare_did_data(dat_use, t)
            
            // 4. Compute estimates for each lead
            for (ll = 1; ll <= n_lead; ll++) {
                est_t = did_fit(dat_did.outcome, dat_did.outcome_delta,
                               dat_did.Gi, dat_did.It,
                               dat_did.covariates, dat_did.id_time_std,
                               lead[ll])
                est_did[iter, ll] = est_t[1]   // tau_DID(t)
                est_sdid[iter, ll] = est_t[2]  // tau_sDID(t)
            }
            
            // 5. Store time weight for this period
            time_weight_new[iter] = time_weight[i]
            
            // 6. Increment iteration counter
            iter++
        }
    }
    
    // Trim to actual number of valid periods
    n_valid = iter - 1
    
    // Handle no valid periods
    if (n_valid == 0) {
        errprintf("sa_compute_did(): No valid periods for SA estimation\n")
        return(result)
    }
    
    est_did = est_did[1::n_valid, .]
    est_sdid = est_sdid[1::n_valid, .]
    time_weight_new = time_weight_new[1::n_valid]
    
    // Compute time-weighted averages (handling missing values separately for DID and sDID)
    for (ll = 1; ll <= n_lead; ll++) {
        // Handle DID missing values separately
        // Use :< . to correctly handle extended missing values (.a-.z)
        valid_idx_did = selectindex(est_did[., ll] :< .)
        
        if (length(valid_idx_did) > 0) {
            // Renormalize weights for DID
            weight_norm_did = time_weight_new[valid_idx_did]
            // Check for zero sum before normalization
            sum_w_did = sum(weight_norm_did)
            if (sum_w_did > 0 & !missing(sum_w_did)) {
                weight_norm_did = weight_norm_did / sum_w_did
                // Time-weighted average for DID
                result.DID[ll] = sum(est_did[valid_idx_did, ll] :* weight_norm_did)
            }
        }
        
        // Handle sDID missing values separately (may differ from DID)
        // Use :< . to correctly handle extended missing values (.a-.z)
        valid_idx_sdid = selectindex(est_sdid[., ll] :< .)
        
        if (length(valid_idx_sdid) > 0) {
            // Renormalize weights for sDID
            weight_norm_sdid = time_weight_new[valid_idx_sdid]
            // Check for zero sum before normalization
            sum_w_sdid = sum(weight_norm_sdid)
            if (sum_w_sdid > 0 & !missing(sum_w_sdid)) {
                weight_norm_sdid = weight_norm_sdid / sum_w_sdid
                // Time-weighted average for sDID
                result.sDID[ll] = sum(est_sdid[valid_idx_sdid, ll] :* weight_norm_sdid)
            }
        }
    }
    
    return(result)
}



/*---------------------------------------------------------------------------
 * subset_data_sa() - Subset Data for Period-Specific Estimation
 *
 * Panel data subset is extracted for specified units and times {t-2, t-1, t}.
 * The subset is prepared for period-specific DID estimation.
 *
 * Arguments:
 *   data     : struct did_data - full panel data
 *   idx_subj : real colvector - valid unit indices (Gmat row indices)
 *   t        : real scalar - treatment period
 *
 * Returns:
 *   struct did_data containing subset with valid units and times
 *
 * Algorithm:
 *   1. Gmat row indices are mapped to unit IDs
 *   2. Observations are filtered for valid units and times {t-2, t-1, t}
 *   3. All relevant fields are copied to result structure
 *---------------------------------------------------------------------------*/
struct did_data scalar subset_data_sa(struct did_data scalar data,
                                       real colvector idx_subj,
                                       real scalar t)
{
    struct did_data scalar result
    real colvector units, valid_units, idx, time_filter, unit_filter, combined_filter
    real scalar i, N_orig, N_sub
    
    // Initialize result structure
    result = did_data()
    
    // Handle empty input
    if (rows(idx_subj) == 0) {
        result.N = 0
        return(result)
    }
    
    N_orig = rows(data.id_unit)
    
    // Step 1: Get unique unit IDs from Gmat row indices
    // idx_subj contains row indices into Gmat (1-based)
    // Map Gmat row indices to actual unit IDs
    units = uniqrows(data.id_unit)
    
    // Validate idx_subj bounds before array access
    if (rows(idx_subj) > 0) {
        if (max(idx_subj) > rows(units) || min(idx_subj) < 1) {
            errprintf("Error: subset_data_sa(): idx_subj contains out-of-bounds indices\n")
            errprintf("       idx_subj range: [%g, %g], units count: %g\n", 
                      min(idx_subj), max(idx_subj), rows(units))
            result.N = 0
            return(result)
        }
    }
    
    valid_units = units[idx_subj]
    
    // Step 2: Create filter for valid units
    transmorphic scalar valid_set
    valid_set = asarray_create("real", 1)
    for (i = 1; i <= rows(valid_units); i++) {
        asarray(valid_set, valid_units[i], 1)
    }
    
    unit_filter = J(N_orig, 1, 0)
    for (i = 1; i <= N_orig; i++) {
        if (asarray_contains(valid_set, data.id_unit[i])) {
            unit_filter[i] = 1
        }
    }
    
    // Step 3: Create filter for valid times {t-2, t-1, t}
    time_filter = (data.id_time :== t) :| (data.id_time :== (t-1)) :| (data.id_time :== (t-2))
    
    // Step 4: Combine filters
    combined_filter = unit_filter :& time_filter
    idx = selectindex(combined_filter)
    
    // Handle empty result
    N_sub = length(idx)
    if (N_sub == 0) {
        result.N = 0
        return(result)
    }
    
    // Step 5: Subset all fields
    result.outcome = data.outcome[idx]
    result.treatment = data.treatment[idx]
    result.id_unit = data.id_unit[idx]
    result.id_time = data.id_time[idx]
    
    // Covariates (if any)
    if (cols(data.covariates) > 0 && rows(data.covariates) > 0) {
        result.covariates = data.covariates[idx, .]
    }
    else {
        result.covariates = J(0, 0, .)
    }
    
    // Cluster variable (if any)
    if (rows(data.cluster_var) > 0) {
        result.cluster_var = data.cluster_var[idx]
    }
    else {
        result.cluster_var = J(0, 1, .)
    }
    
    // Derived variables (if available)
    if (rows(data.Gi) > 0) {
        result.Gi = data.Gi[idx]
    }
    if (rows(data.It) > 0) {
        result.It = data.It[idx]
    }
    if (rows(data.id_time_std) > 0) {
        result.id_time_std = data.id_time_std[idx]
    }
    if (rows(data.outcome_delta) > 0) {
        result.outcome_delta = data.outcome_delta[idx]
    }
    
    // Metadata
    result.N = N_sub
    result.n_units = rows(uniqrows(result.id_unit))
    result.n_periods = rows(uniqrows(result.id_time))
    result.treat_year = t
    result.is_panel = data.is_panel
    
    return(result)
}


/*---------------------------------------------------------------------------
 * add_lead_outcomes() - Append Future Outcome Observations
 *
 * Outcome observations from future periods (t+1 to t+max_lead) are appended
 * when lead > 0. Treatment status from time t is preserved.
 *
 * Arguments:
 *   dat_use  : struct did_data - current subset data (times {t-2, t-1, t})
 *   data     : struct did_data - full panel data
 *   idx_subj : real colvector - valid unit indices
 *   t        : real scalar - treatment period
 *   lead     : real rowvector - lead parameters for dynamic effects
 *
 * Returns:
 *   struct did_data with lead outcome observations appended
 *
 * Algorithm:
 *   1. Treatment status at time t is extracted for valid units
 *   2. Observations are filtered for lead time periods
 *   3. Treatment info is joined by unit ID
 *   4. Lead observations are appended to existing data structure
 *---------------------------------------------------------------------------*/
struct did_data scalar add_lead_outcomes(struct did_data scalar dat_use,
                                          struct did_data scalar data,
                                          real colvector idx_subj,
                                          real scalar t,
                                          real rowvector lead)
{
    struct did_data scalar result
    real colvector units, valid_units, idx_treat, idx_lead
    real colvector unit_filter, time_filter, combined_filter
    real colvector treat_info_unit, treat_info_val
    real colvector lead_outcome, lead_treatment, lead_id_unit, lead_id_time
    real matrix lead_covariates
    real colvector lead_cluster
    real scalar N_orig, N_use, N_lead, N_total
    real scalar min_lead_time, max_lead_time, i, j
    real scalar max_lead, min_lead
    
    // Get lead bounds
    max_lead = max(lead)
    min_lead = min(lead)
    
    // If max_lead <= 0, no lead data needed
    if (max_lead <= 0) {
        return(dat_use)
    }
    
    N_orig = rows(data.id_unit)
    N_use = dat_use.N
    
    // Step 1: Get unique unit IDs from Gmat row indices
    units = uniqrows(data.id_unit)
    valid_units = units[idx_subj]
    
    // Step 2: Get treatment info at time t for valid subjects
    {
        transmorphic scalar valid_set, treat_map
        string scalar key
        
        // Build valid_units set for efficient lookup
        valid_set = asarray_create("real", 1)
        for (i = 1; i <= rows(valid_units); i++) {
            asarray(valid_set, valid_units[i], 1)
        }
        
        // Filter observations to valid units
        unit_filter = J(N_orig, 1, 0)
        for (i = 1; i <= N_orig; i++) {
            if (asarray_contains(valid_set, data.id_unit[i])) {
                unit_filter[i] = 1
            }
        }
        
        time_filter = (data.id_time :== t)
        combined_filter = unit_filter :& time_filter
        idx_treat = selectindex(combined_filter)
        
        // Build treatment info map for efficient lookup
        treat_map = asarray_create("real", 1)
        for (i = 1; i <= rows(idx_treat); i++) {
            asarray(treat_map, data.id_unit[idx_treat[i]], data.treatment[idx_treat[i]])
        }
        
        // Step 3: Get outcome observations for lead times
        min_lead_time = t + (min_lead > 1 ? min_lead : 1)
        max_lead_time = t + max_lead
        
        // Explicitly exclude missing values in time filter
        // Required because Mata treats missing >= x as true
        time_filter = (data.id_time :< .) :& (data.id_time :>= min_lead_time) :& (data.id_time :<= max_lead_time)
        combined_filter = unit_filter :& time_filter
        idx_lead = selectindex(combined_filter)
        
        // Handle no lead observations
        N_lead = length(idx_lead)
        if (N_lead == 0) {
            return(dat_use)
        }
        
        // Step 4: Extract lead data and join treatment info
        lead_outcome = data.outcome[idx_lead]
        lead_id_unit = data.id_unit[idx_lead]
        lead_id_time = data.id_time[idx_lead]
        
        // Join treatment info by unit ID
        lead_treatment = J(N_lead, 1, .)
        for (i = 1; i <= N_lead; i++) {
            if (asarray_contains(treat_map, lead_id_unit[i])) {
                lead_treatment[i] = asarray(treat_map, lead_id_unit[i])
            }
        }
    }
    
    // Handle covariates
    if (cols(data.covariates) > 0 && rows(data.covariates) > 0) {
        lead_covariates = data.covariates[idx_lead, .]
    }
    else {
        lead_covariates = J(0, 0, .)
    }
    
    // Handle cluster variable
    if (rows(data.cluster_var) > 0) {
        lead_cluster = data.cluster_var[idx_lead]
    }
    else {
        lead_cluster = J(0, 1, .)
    }
    
    // Step 5: Append lead data to dat_use
    result = did_data()
    N_total = N_use + N_lead
    
    result.outcome = dat_use.outcome \ lead_outcome
    result.treatment = dat_use.treatment \ lead_treatment
    result.id_unit = dat_use.id_unit \ lead_id_unit
    result.id_time = dat_use.id_time \ lead_id_time
    
    // Covariates
    if (cols(dat_use.covariates) > 0 && cols(lead_covariates) > 0) {
        result.covariates = dat_use.covariates \ lead_covariates
    }
    else if (cols(dat_use.covariates) > 0) {
        result.covariates = dat_use.covariates \ J(N_lead, cols(dat_use.covariates), .)
    }
    else if (cols(lead_covariates) > 0) {
        result.covariates = J(N_use, cols(lead_covariates), .) \ lead_covariates
    }
    else {
        result.covariates = J(0, 0, .)
    }
    
    // Cluster variable
    if (rows(dat_use.cluster_var) > 0 && rows(lead_cluster) > 0) {
        result.cluster_var = dat_use.cluster_var \ lead_cluster
    }
    else if (rows(dat_use.cluster_var) > 0) {
        result.cluster_var = dat_use.cluster_var \ J(N_lead, 1, .)
    }
    else if (rows(lead_cluster) > 0) {
        result.cluster_var = J(N_use, 1, .) \ lead_cluster
    }
    else {
        result.cluster_var = J(0, 1, .)
    }
    
    // Derived variables (will be recomputed by sa_prepare_did_data)
    result.Gi = J(0, 1, .)
    result.It = J(0, 1, .)
    result.id_time_std = J(0, 1, .)
    result.outcome_delta = J(0, 1, .)
    
    // Metadata
    result.N = N_total
    result.n_units = rows(uniqrows(result.id_unit))
    result.n_periods = rows(uniqrows(result.id_time))
    result.treat_year = t
    result.is_panel = dat_use.is_panel
    
    return(result)
}


/*---------------------------------------------------------------------------
 * sa_prepare_did_data() - Compute Derived Variables for DID Estimation
 *
 * Derived variables required for DID estimation are computed:
 *   - G_i: group indicator (treatment status)
 *   - I_t: post-treatment indicator
 *   - id_time_std: standardized time relative to treatment
 *   - Delta_Y: outcome change from lagged group mean
 *
 * Arguments:
 *   data : struct did_data - subset data from previous processing
 *   t    : real scalar - treatment period (for time standardization)
 *
 * Returns:
 *   struct did_data with derived variables computed
 *
 * Algorithm:
 *   1. G_i = max(treatment) per unit is computed (0 = control, 1 = treated)
 *   2. I_t = 1{id_time >= t} is computed
 *   3. id_time_std = id_time - t is computed
 *   4. Delta_Y = Y - Y_bar(G_i, id_time_std - 1) is computed
 *      where Y_bar is the lagged group mean
 *---------------------------------------------------------------------------*/
struct did_data scalar sa_prepare_did_data(struct did_data scalar data,
                                            real scalar t)
{
    struct did_data scalar result
    real colvector Gi, It, id_time_std, outcome_delta
    real scalar N, j, g, ts
    real scalar mean_y
    transmorphic scalar unit_map, group_mean_map, group_has_na_map
    real colvector idx
    string scalar key, lag_key
    real scalar sum_y, n_obs
    
    // Copy input data
    result = data
    N = data.N
    
    // Handle empty data
    if (N == 0) {
        result.Gi = J(0, 1, .)
        result.It = J(0, 1, .)
        result.id_time_std = J(0, 1, .)
        result.outcome_delta = J(0, 1, .)
        return(result)
    }
    
    // Step 1: Compute Gi (group indicator)
    // Gi = 1 if unit is ever treated, 0 otherwise
    unit_map = asarray_create("real", 1)
    for (j = 1; j <= N; j++) {
        if (asarray_contains(unit_map, data.id_unit[j])) {
            if (data.treatment[j] > asarray(unit_map, data.id_unit[j])) {
                asarray(unit_map, data.id_unit[j], data.treatment[j])
            }
        }
        else {
            asarray(unit_map, data.id_unit[j], data.treatment[j])
        }
    }
    
    // Map Gi back to observations
    Gi = J(N, 1, .)
    for (j = 1; j <= N; j++) {
        Gi[j] = asarray(unit_map, data.id_unit[j])
    }
    
    // Step 2: Compute It (post-treatment indicator)
    // Explicitly exclude missing values (required for correct comparison)
    It = (data.id_time :< .) :& (data.id_time :>= t)
    
    // Step 3: Compute id_time_std (standardized time relative to treatment)
    id_time_std = data.id_time :- t
    
    // Step 4: Compute outcome_delta (Delta_Y = Y - lag_group_mean(Y))
    // Group mean propagates missing values (if any missing in group, mean is missing)
    
    // Build (Gi, id_time_std) -> (sum, count) map using hash table
    group_mean_map = asarray_create("string", 1)
    group_has_na_map = asarray_create("string", 1)  // Track missing value presence per group
    
    for (j = 1; j <= N; j++) {
        // Skip observations with missing key components
        if (!missing(Gi[j]) && !missing(id_time_std[j])) {
            key = strofreal(Gi[j]) + "_" + strofreal(id_time_std[j])
            
            // Track if group contains any missing value (propagate to group mean)
            if (missing(data.outcome[j])) {
                // Mark this group as having missing value
                asarray(group_has_na_map, key, 1)
                // Still need to initialize group if not exists
                if (!asarray_contains(group_mean_map, key)) {
                    asarray(group_mean_map, key, (0, 0))
                }
            }
            else {
                // Accumulate sum and count for non-missing values
                if (asarray_contains(group_mean_map, key)) {
                    idx = asarray(group_mean_map, key)
                    idx[1] = idx[1] + data.outcome[j]  // sum
                    idx[2] = idx[2] + 1                 // count
                    asarray(group_mean_map, key, idx)
                }
                else {
                    asarray(group_mean_map, key, (data.outcome[j], 1))
                }
            }
        }
    }
    
    // Compute outcome_delta: Delta_Y = Y - Y_bar(G_i, t-1)
    outcome_delta = J(N, 1, .)
    
    for (j = 1; j <= N; j++) {
        // Skip observations with missing Gi or id_time_std
        if (missing(Gi[j]) || missing(id_time_std[j])) {
            continue
        }
        
        g = Gi[j]
        ts = id_time_std[j] - 1  // Lag by 1 period
        
        // Look up mean for (g, ts)
        lag_key = strofreal(g) + "_" + strofreal(ts)
        if (asarray_contains(group_mean_map, lag_key)) {
            // If lag group has any missing value, mean is missing (outcome_delta remains missing)
            if (asarray_contains(group_has_na_map, lag_key)) {
                // Group contains missing value, outcome_delta remains missing
            }
            else {
                idx = asarray(group_mean_map, lag_key)
                if (idx[2] > 0) {
                    mean_y = idx[1] / idx[2]  // sum / count
                    
                    if (!missing(data.outcome[j])) {
                        outcome_delta[j] = data.outcome[j] - mean_y
                    }
                }
            }
        }
    }
    
    // Store derived variables
    result.Gi = Gi
    result.It = It
    result.id_time_std = id_time_std
    result.outcome_delta = outcome_delta
    result.treat_year = t
    
    return(result)
}


// ----------------------------------------------------------------------------
// SA MAIN ESTIMATION FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * sa_estimate() - SA Design Main Entry Point
 *
 * The complete SA estimation workflow is orchestrated: data validation,
 * point estimation, bootstrap variance estimation, and GMM-based double
 * DID aggregation.
 *
 * Arguments:
 *   data   : struct did_data - panel data (normalized unit/time indices)
 *   option : struct did_option - estimation options:
 *            - thres: minimum treated units per period (default=2)
 *            - lead: lead values for dynamic effects (default=0)
 *            - n_boot: bootstrap iterations (default=30)
 *            - level: confidence level (default=95)
 *            - quiet: suppress progress display (default=0)
 *
 * Returns:
 *   struct sa_ddid_result containing:
 *     - estimate: SA-Double-DID estimates (tau_dDID)
 *     - tau_did, tau_sdid: SA-DID and SA-sDID estimates
 *     - variance, std_error: bootstrap variance and standard errors
 *     - ci_low, ci_high: bootstrap percentile confidence intervals
 *     - w_did, w_sdid: GMM optimal weights
 *     - W_matrices: precision matrices
 *
 * Algorithm:
 *   1. Data normalization is validated (indices start from 1, consecutive)
 *   2. Gmat is constructed and valid periods are validated
 *   3. Point estimation is performed via sa_double_did()
 *   4. Panel bootstrap: resample units, recompute full estimation
 *   5. GMM aggregation via sa_to_ddid(point_est, boot_est, lead)
 *---------------------------------------------------------------------------*/
struct sa_ddid_result scalar sa_estimate(struct did_data scalar data,
                                          struct did_option scalar option)
{
    struct sa_ddid_result scalar result
    struct sa_point scalar point_est
    struct sa_point scalar boot_pt
    pointer(struct sa_point scalar) vector boot_est
    struct did_data scalar dat_boot
    
    real matrix Gmat
    real colvector id_time_use, time_weight
    pointer vector id_subj_use
    
    real scalar n_boot, n_lead, b, n_boot_success
    real scalar progress_freq
    real colvector unique_times, unique_units  // For consecutive validation
    real scalar has_valid_estimate, ll  // For multi-lead validation
    
    // Step 0: Initialize and validate
    n_boot = option.n_boot
    n_lead = cols(option.lead)
    
    // Validate input data is panel format
    if (data.is_panel != 1) {
        errprintf("sa_estimate(): SA design requires panel data (not RCS)\n")
        result = sa_ddid_result()  // Initialize for error return
        return(result)
    }
    
    // Validate data has required fields
    if (data.N == 0) {
        errprintf("sa_estimate(): No observations in data\n")
        result = sa_ddid_result()  // Initialize for error return
        return(result)
    }
    
    // Step 1: Data preparation verification
    // Verify id_time is normalized (should start from 1)
    if (min(data.id_time) != 1) {
        errprintf("sa_estimate(): id_time should be normalized to start from 1\n")
        errprintf("               Found min(id_time) = %g\n", min(data.id_time))
        result = sa_ddid_result()  // Initialize for error return
        return(result)
    }
    
    // Verify id_unit is normalized (should start from 1)
    if (min(data.id_unit) != 1) {
        errprintf("sa_estimate(): id_unit should be normalized to start from 1\n")
        errprintf("               Found min(id_unit) = %g\n", min(data.id_unit))
        result = sa_ddid_result()  // Initialize for error return
        return(result)
    }
    
    // Verify id_time is consecutive integers 1, 2, ..., T
    unique_times = uniqrows(data.id_time)
    if (rows(unique_times) != max(data.id_time)) {
        errprintf("sa_estimate(): id_time must be consecutive integers 1, 2, ..., T\n")
        errprintf("               Found %g unique values but max = %g\n", 
                  rows(unique_times), max(data.id_time))
        result = sa_ddid_result()
        return(result)
    }
    
    // Also verify id_unit is consecutive
    unique_units = uniqrows(data.id_unit)
    if (rows(unique_units) != max(data.id_unit)) {
        errprintf("sa_estimate(): id_unit must be consecutive integers 1, 2, ..., N\n")
        errprintf("               Found %g unique values but max = %g\n",
                  rows(unique_units), max(data.id_unit))
        result = sa_ddid_result()
        return(result)
    }
    
    // Step 2: Compute Gmat and validate
    Gmat = create_gmat(data.id_unit, data.id_time, data.treatment)
    
    if (rows(Gmat) == 0 || cols(Gmat) == 0) {
        errprintf("sa_estimate(): Failed to create Gmat\n")
        result = sa_ddid_result()  // Initialize for error return
        return(result)
    }
    
    id_time_use = get_periods(Gmat, option.thres)
    
    if (rows(id_time_use) == 0) {
        errprintf("sa_estimate(): No valid periods found with threshold = %g\n", option.thres)
        result = sa_ddid_result()
        return(result)
    }
    
    // Get time weights for validation
    time_weight = get_time_weight(Gmat, id_time_use)
    
    // Verify weights sum to 1.0
    if (abs(sum(time_weight) - 1.0) > 1e-10) {
        errprintf("sa_estimate(): Time weights do not sum to 1.0 (sum = %g)\n", sum(time_weight))
        result = sa_ddid_result()
        return(result)
    }
    
    // Step 3: Point estimation on original data
    point_est = sa_double_did(data, option)
    
    // Check if ANY lead has valid estimates (not just the first one)
    has_valid_estimate = 0
    for (ll = 1; ll <= n_lead; ll++) {
        if (!missing(point_est.DID[ll]) || !missing(point_est.sDID[ll])) {
            has_valid_estimate = 1
            break
        }
    }
    if (!has_valid_estimate) {
        errprintf("sa_estimate(): Point estimation failed for all leads\n")
        result = sa_ddid_result()
        return(result)
    }
    
    // Step 4: Bootstrap loop for variance estimation
    // Each iteration fully recomputes Gmat, periods, subjects, and weights
    
    // Handle n_boot = 0 case (skip bootstrap)
    if (n_boot == 0) {
        errprintf("sa_estimate(): n_boot = 0, bootstrap skipped (confidence intervals unavailable)\n")
        // Return point estimates only
        result.tau_did = point_est.DID
        result.tau_sdid = point_est.sDID
        return(result)
    }
    
    // Pre-allocate pointer vector for bootstrap results
    boot_est = J(n_boot, 1, NULL)
    n_boot_success = 0
    progress_freq = 10
    
    for (b = 1; b <= n_boot; b++) {
        // Progress display
        if (option.quiet == 0 && mod(b, progress_freq) == 0) {
            printf("Bootstrap: %g/%g (%g%%)\n", b, n_boot, round(100*b/n_boot))
            displayflush()
        }
        
        // Resample panel data (unit-level resampling)
        dat_boot = sample_panel(data)
        
        // Compute SA estimates on bootstrap sample
        // (sa_double_did() internally recomputes Gmat, periods, subjects, weights)
        boot_pt = sa_double_did(dat_boot, option)
        
        // Store result if valid (check ALL leads, not just first one)
        has_valid_estimate = 0
        for (ll = 1; ll <= n_lead; ll++) {
            if (!missing(boot_pt.DID[ll]) || !missing(boot_pt.sDID[ll])) {
                has_valid_estimate = 1
                break
            }
        }
        if (has_valid_estimate) {
            // Allocate new struct and store pointer
            boot_est[b] = &(sa_point())
            (*boot_est[b]).DID = boot_pt.DID
            (*boot_est[b]).sDID = boot_pt.sDID
            n_boot_success++
        }
    }
    
    // Final progress display
    if (option.quiet == 0) {
        printf("Bootstrap: %g/%g (100%%)\n", n_boot, n_boot)
        displayflush()
    }
    
    // Check for sufficient bootstrap samples
    if (n_boot_success < 2) {
        errprintf("sa_estimate(): Bootstrap failed - only %g of %g iterations succeeded\n",
                  n_boot_success, n_boot)
        // Return point estimates only (initialize result for this error path)
        result = sa_ddid_result()
        result.tau_did = point_est.DID
        result.tau_sdid = point_est.sDID
        return(result)
    }
    
    // Warn if some bootstrap iterations failed
    if (n_boot_success < n_boot) {
        printf("Warning: %g of %g bootstrap iterations failed\n", 
               n_boot - n_boot_success, n_boot)
    }
    
    // Step 5: Compute Double DID with GMM weights
    result = sa_to_ddid(point_est, boot_est, option.lead, option.level)
    
    return(result)
}


/*---------------------------------------------------------------------------
 * _did_sa_main() - SA Estimation Entry Point (Ado Interface)
 *
 * This function serves as the entry point for SA estimation called from
 * _diddesign_sa.ado. The full SA estimation pipeline is implemented and
 * global result matrices are populated for retrieval by the ado file.
 *
 * Note: For direct Mata usage, sa_estimate() is preferred as it returns
 * a structured result without side effects.
 *
 * Arguments:
 *   lead    : real rowvector - lead values for dynamic effects
 *   n_boot  : real scalar - number of bootstrap iterations
 *   thres   : real scalar - threshold for valid periods
 *   level   : real scalar - confidence level (default: 95)
 *   quiet   : real scalar - suppress progress display (default: 0)
 *
 * Side Effects:
 *   Global matrices are populated: _sa_b, _sa_V, _sa_estimates, _sa_weights, etc.
 *---------------------------------------------------------------------------*/
real scalar _did_sa_main(real rowvector lead, real scalar n_boot, 
                          real scalar thres, real scalar level,
                          | real scalar quiet)
{
    // Declare external global data structure
    external struct did_data scalar did_dat
    external struct did_option scalar did_opt
    
    // Declare external global result variables for SA
    external real rowvector _sa_b
    external real matrix _sa_V
    external real matrix _sa_estimates
    external real rowvector _sa_lead_values
    external real matrix _sa_weights
    external real matrix _sa_W
    external real matrix _sa_vcov_gmm
    external real scalar _sa_n_boot_success
    external real scalar _sa_n_periods_valid
    external real matrix _sa_time_weights
    
    struct sa_point scalar point_est
    struct sa_point scalar boot_pt
    pointer(struct sa_point scalar) vector boot_est
    struct sa_ddid_result scalar ddid_res
    struct did_data scalar dat_boot
    
    real matrix Gmat
    real colvector id_time_use, time_weight
    pointer vector id_subj_use
    
    real scalar n_lead, l, row, b
    real scalar tau_did, tau_sdid, tau_ddid
    real scalar var_did, var_sdid, var_ddid
    real scalar se_did, se_sdid, se_ddid
    real scalar ci_lo_did, ci_hi_did, ci_lo_sdid, ci_hi_sdid
    real scalar ci_lo_ddid, ci_hi_ddid
    real scalar w_did, w_sdid
    real scalar n_boot_success
    real scalar progress_freq
    real scalar has_valid_estimate, boot_has_valid, ll_check
    
    // Step 0: Set default parameters
    if (args() < 4) level = 95
    if (args() < 5) quiet = 0
    
    // Store options in did_opt for use by sa_double_did
    did_opt.thres = thres
    did_opt.lead = lead
    did_opt.n_boot = n_boot
    did_opt.level = level
    did_opt.quiet = quiet
    
    // Step 1: Initialize result storage
    n_lead = cols(lead)
    
    _sa_b = J(1, 3 * n_lead, .)
    _sa_V = J(3 * n_lead, 3 * n_lead, 0)
    _sa_estimates = J(3 * n_lead, 6, .)
    _sa_lead_values = lead
    _sa_weights = J(n_lead, 2, .)
    _sa_W = J(n_lead, 4, .)           // Flattened 2x2 matrices
    _sa_vcov_gmm = J(n_lead, 4, .)    // Flattened 2x2 matrices
    _sa_n_boot_success = n_boot
    _sa_n_periods_valid = 0
    _sa_time_weights = J(0, 1, .)
    
    // Step 2: Compute Gmat and valid periods
    Gmat = create_gmat(did_dat.id_unit, did_dat.id_time, did_dat.treatment)
    
    if (rows(Gmat) == 0 || cols(Gmat) == 0) {
        errprintf("_did_sa_main(): Failed to create Gmat\n")
        return(1)
    }
    
    id_time_use = get_periods(Gmat, thres)
    _sa_n_periods_valid = rows(id_time_use)
    
    if (_sa_n_periods_valid == 0) {
        errprintf("_did_sa_main(): No valid periods found with threshold = %g\n", thres)
        return(2)
    }
    
    // Get time weights for reporting
    time_weight = get_time_weight(Gmat, id_time_use)
    _sa_time_weights = time_weight
    
    // Step 3: Obtain point estimates
    point_est = sa_double_did(did_dat, did_opt)
    
    // Check if ANY lead has valid estimates (not just the first one)
    has_valid_estimate = 0
    for (ll_check = 1; ll_check <= n_lead; ll_check++) {
        if (!missing(point_est.DID[ll_check]) || !missing(point_est.sDID[ll_check])) {
            has_valid_estimate = 1
            break
        }
    }
    if (!has_valid_estimate) {
        errprintf("_did_sa_main(): Point estimation failed for all leads\n")
        return(3)
    }
    
    // Step 4: Bootstrap loop for variance estimation
    boot_est = J(n_boot, 1, NULL)
    n_boot_success = 0
    progress_freq = max((1, floor(n_boot / 10)))
    
    for (b = 1; b <= n_boot; b++) {
        // Progress display
        if (quiet == 0 && mod(b, progress_freq) == 0) {
            printf("{txt}Bootstrap: %g/%g (%g%%)\n", b, n_boot, round(100*b/n_boot))
            displayflush()
        }
        
        // Sample panel data (block bootstrap by unit)
        dat_boot = sample_panel(did_dat)
        
        // Run SA estimation on bootstrap sample
        boot_pt = sa_double_did(dat_boot, did_opt)
        
        // Check if ANY lead has valid estimates
        boot_has_valid = 0
        for (ll_check = 1; ll_check <= n_lead; ll_check++) {
            if (!missing(boot_pt.DID[ll_check]) || !missing(boot_pt.sDID[ll_check])) {
                boot_has_valid = 1
                break
            }
        }
        if (boot_has_valid) {
            // Allocate new struct and store pointer
            boot_est[b] = &(sa_point())
            (*boot_est[b]).DID = boot_pt.DID
            (*boot_est[b]).sDID = boot_pt.sDID
            n_boot_success++
        }
    }
    
    // Final progress display
    if (quiet == 0) {
        printf("{txt}Bootstrap: %g/%g (100%%)\n", n_boot, n_boot)
        displayflush()
    }
    
    _sa_n_boot_success = n_boot_success
    
    // Check for sufficient bootstrap samples
    if (n_boot_success < 2) {
        errprintf("_did_sa_main(): Bootstrap failed - only %g of %g iterations succeeded\n",
                  n_boot_success, n_boot)
        return(4)
    }
    
    // Step 5: Compute Double DID via GMM
    ddid_res = sa_to_ddid(point_est, boot_est, lead, level)
    
    // Step 6: Store results in global matrices
    row = 1
    for (l = 1; l <= n_lead; l++) {
        
        // Extract results for this lead
        tau_ddid = ddid_res.estimate[l]
        tau_did = ddid_res.tau_did[l]
        tau_sdid = ddid_res.tau_sdid[l]
        
        var_ddid = ddid_res.variance[l]
        var_did = ddid_res.var_did[l]
        var_sdid = ddid_res.var_sdid[l]
        
        se_ddid = ddid_res.std_error[l]
        se_did = sqrt(var_did)
        se_sdid = sqrt(var_sdid)
        
        w_did = ddid_res.w_did[l]
        w_sdid = ddid_res.w_sdid[l]
        
        ci_lo_ddid = ddid_res.ci_low[l]
        ci_hi_ddid = ddid_res.ci_high[l]
        
        // Store weights
        _sa_weights[l, .] = (w_did, w_sdid)
        
        // Store W matrix (flattened) if available
        if (ddid_res.W_matrices[l] != NULL) {
            _sa_W[l, .] = vec(*ddid_res.W_matrices[l])'
        }
        
        // Store VCOV matrix (flattened) if available
        if (ddid_res.VCOV_matrices[l] != NULL) {
            _sa_vcov_gmm[l, .] = vec(*ddid_res.VCOV_matrices[l])'
        }
        else if (!missing(var_did) && !missing(var_sdid)) {
            // Fallback: Store only diagonal elements if VCOV not available
            _sa_vcov_gmm[l, 1] = var_did
            _sa_vcov_gmm[l, 2] = .  // Cov(DID, sDID) - not available
            _sa_vcov_gmm[l, 3] = .  // Cov(DID, sDID) - not available
            _sa_vcov_gmm[l, 4] = var_sdid
        }
        
        // Compute bootstrap percentile CIs for DID and sDID
        // (sa_to_ddid() returns only dDID CI; DID/sDID CIs computed separately)
        ci_lo_did = .
        ci_hi_did = .
        ci_lo_sdid = .
        ci_hi_sdid = .
        
        _sa_compute_bootstrap_ci(boot_est, l, level, &ci_lo_did, &ci_hi_did, 
                                 &ci_lo_sdid, &ci_hi_sdid)
        
        // ---------------------------------------------------------------------
        // Store results in matrices
        // ---------------------------------------------------------------------
        
        // e(b): coefficient vector [SA_dDID, SA_DID, SA_sDID] for each lead
        _sa_b[1, 3*(l-1)+1] = tau_ddid
        _sa_b[1, 3*(l-1)+2] = tau_did
        _sa_b[1, 3*(l-1)+3] = tau_sdid
        
        // e(V): variance-covariance matrix (diagonal)
        _sa_V[3*(l-1)+1, 3*(l-1)+1] = var_ddid
        _sa_V[3*(l-1)+2, 3*(l-1)+2] = var_did
        _sa_V[3*(l-1)+3, 3*(l-1)+3] = var_sdid
        
        // e(estimates): full results table
        // Row order: SA_dDID, SA_DID, SA_sDID for each lead
        // Columns: lead, estimate, std.error, ci_lo, ci_hi, weight
        
        // SA-Double-DID row
        _sa_estimates[row, 1] = lead[l]
        _sa_estimates[row, 2] = tau_ddid
        _sa_estimates[row, 3] = se_ddid
        _sa_estimates[row, 4] = ci_lo_ddid
        _sa_estimates[row, 5] = ci_hi_ddid
        _sa_estimates[row, 6] = .  // No weight for dDID
        row = row + 1
        
        // SA-DID row
        _sa_estimates[row, 1] = lead[l]
        _sa_estimates[row, 2] = tau_did
        _sa_estimates[row, 3] = se_did
        _sa_estimates[row, 4] = ci_lo_did
        _sa_estimates[row, 5] = ci_hi_did
        _sa_estimates[row, 6] = w_did
        row = row + 1
        
        // SA-sDID row
        _sa_estimates[row, 1] = lead[l]
        _sa_estimates[row, 2] = tau_sdid
        _sa_estimates[row, 3] = se_sdid
        _sa_estimates[row, 4] = ci_lo_sdid
        _sa_estimates[row, 5] = ci_hi_sdid
        _sa_estimates[row, 6] = w_sdid
        row = row + 1
    }
    
    return(0)  // Success
}


/*---------------------------------------------------------------------------
 * _sa_compute_bootstrap_ci() - Compute Bootstrap Percentile CIs
 *
 * Bootstrap percentile confidence intervals are computed for SA-DID and
 * SA-sDID estimators at a specified lead index.
 *
 * Arguments:
 *   boot_est   : pointer vector - bootstrap sa_point structures
 *   lead_idx   : real scalar - lead index (1-based)
 *   level      : real scalar - confidence level (e.g., 95)
 *   ci_lo_did  : pointer(real scalar) - output: lower CI for DID
 *   ci_hi_did  : pointer(real scalar) - output: upper CI for DID
 *   ci_lo_sdid : pointer(real scalar) - output: lower CI for sDID
 *   ci_hi_sdid : pointer(real scalar) - output: upper CI for sDID
 *---------------------------------------------------------------------------*/
void _sa_compute_bootstrap_ci(pointer(struct sa_point scalar) vector boot_est,
                               real scalar lead_idx,
                               real scalar level,
                               pointer(real scalar) scalar ci_lo_did,
                               pointer(real scalar) scalar ci_hi_did,
                               pointer(real scalar) scalar ci_lo_sdid,
                               pointer(real scalar) scalar ci_hi_sdid)
{
    struct sa_point scalar pt
    real colvector boot_did, boot_sdid
    real scalar n_boot, b, alpha
    
    n_boot = length(boot_est)
    alpha = 1 - level / 100
    
    // Initialize with missing
    *ci_lo_did = .
    *ci_hi_did = .
    *ci_lo_sdid = .
    *ci_hi_sdid = .
    
    if (n_boot < 2) {
        return
    }
    
    // Extract bootstrap estimates
    boot_did = J(n_boot, 1, .)
    boot_sdid = J(n_boot, 1, .)
    
    for (b = 1; b <= n_boot; b++) {
        if (boot_est[b] == NULL) {
            continue
        }
        
        pt = *boot_est[b]
        
        // Check lead_idx bounds before accessing vectors
        if (lead_idx < 1 || lead_idx > cols(pt.DID) || lead_idx > cols(pt.sDID)) {
            continue
        }
        
        boot_did[b] = pt.DID[lead_idx]
        boot_sdid[b] = pt.sDID[lead_idx]
    }
    
    // Compute quantiles for DID
    *ci_lo_did = quantile_sorted(boot_did, alpha / 2)
    *ci_hi_did = quantile_sorted(boot_did, 1 - alpha / 2)
    
    // Compute quantiles for sDID
    *ci_lo_sdid = quantile_sorted(boot_sdid, alpha / 2)
    *ci_hi_sdid = quantile_sorted(boot_sdid, 1 - alpha / 2)
}


/*---------------------------------------------------------------------------
 * _did_sa_prepare_data() - Load and Normalize Data for SA Estimation
 *
 * Data is loaded from Stata into the global did_dat structure and unit/time
 * identifiers are normalized to consecutive integers (1, 2, 3, ...).
 * This function is called from _diddesign_sa.ado before _did_sa_main().
 *
 * Arguments:
 *   outcome_var   : string - outcome variable name
 *   treatment_var : string - treatment indicator variable name
 *   id_var        : string - unit identifier variable name
 *   time_var      : string - time identifier variable name
 *   cluster_var   : string - cluster variable name (or empty)
 *   covariates    : string - space-separated covariate names
 *   touse_var     : string - sample marker variable name
 *
 * Side Effects:
 *   Global did_dat structure is populated with normalized data
 *---------------------------------------------------------------------------*/
real scalar _did_sa_prepare_data(string scalar outcome_var,
                                  string scalar treatment_var,
                                  string scalar id_var,
                                  string scalar time_var,
                                  string scalar cluster_var,
                                  string scalar covariates,
                                  string scalar touse_var)
{
    // Declare external global data structure
    external struct did_data scalar did_dat
    external struct did_option scalar did_opt
    
    real colvector outcome, treatment, id_unit, id_time, cluster
    real matrix covars
    real colvector unique_units, unique_times
    real colvector id_unit_norm, id_time_norm
    real scalar N, i, n_units, n_times, j
    string rowvector covar_names
    real scalar n_covars
    
    // Step 0: Initialize global structures
    did_dat = did_data()
    did_opt = init_did_option()
    
    // Step 1: Load data from Stata
    
    // Load main variables (only for touse == 1)
    outcome = st_data(., outcome_var, touse_var)
    treatment = st_data(., treatment_var, touse_var)
    id_unit = st_data(., id_var, touse_var)
    id_time = st_data(., time_var, touse_var)
    
    N = rows(outcome)
    
    if (N == 0) {
        errprintf("_did_sa_prepare_data(): No observations selected\n")
        return(1)
    }
    
    // Load cluster variable
    if (cluster_var != "") {
        cluster = st_data(., cluster_var, touse_var)
    }
    else {
        cluster = id_unit  // Default to unit ID
    }
    
    // Load covariates if specified
    if (covariates != "") {
        covar_names = tokens(covariates)
        n_covars = cols(covar_names)
        covars = st_data(., covariates, touse_var)
    }
    else {
        n_covars = 0
        covars = J(0, 0, .)
    }
    
    // Step 2: Normalize id_unit and id_time to sequential integers (1, 2, ...)
    unique_units = uniqrows(id_unit)
    unique_times = uniqrows(id_time)
    n_units = rows(unique_units)
    n_times = rows(unique_times)
    
    // Build mapping: original ID -> normalized index
    transmorphic scalar unit_map, time_map
    unit_map = asarray_create("real")
    time_map = asarray_create("real")
    
    for (j = 1; j <= n_units; j++) {
        asarray(unit_map, unique_units[j], j)
    }
    for (j = 1; j <= n_times; j++) {
        asarray(time_map, unique_times[j], j)
    }
    
    // Create normalized IDs
    id_unit_norm = J(N, 1, .)
    id_time_norm = J(N, 1, .)
    
    for (i = 1; i <= N; i++) {
        id_unit_norm[i] = asarray(unit_map, id_unit[i])
        id_time_norm[i] = asarray(time_map, id_time[i])
    }
    
    // Step 3: Populate global did_dat structure
    did_dat.outcome = outcome
    did_dat.treatment = treatment
    did_dat.id_unit = id_unit_norm
    did_dat.id_time = id_time_norm
    did_dat.cluster_var = cluster
    did_dat.covariates = covars
    
    // Metadata
    did_dat.N = N
    did_dat.n_units = n_units
    did_dat.n_periods = n_times
    did_dat.is_panel = 1
    
    // Derived variables will be computed during estimation
    did_dat.Gi = J(0, 1, .)
    did_dat.It = J(0, 1, .)
    did_dat.id_time_std = J(0, 1, .)
    did_dat.outcome_delta = J(0, 1, .)
    did_dat.treat_year = .
    
    return(0)  // Success
}


// ----------------------------------------------------------------------------
// MODULE VERIFICATION FUNCTION
// ----------------------------------------------------------------------------

/*---------------------------------------------------------------------------
 * _did_sa_loaded() - Module Load Verification
 *
 * A confirmation message is displayed when did_sa.mata is loaded successfully.
 *---------------------------------------------------------------------------*/
void _did_sa_loaded()
{
    printf("{txt}did_sa.mata loaded successfully\n")
    printf("{txt}  - sa_estimate(): SA design main entry point\n")
    printf("{txt}  - sa_double_did(): SA point estimation\n")
}

end
