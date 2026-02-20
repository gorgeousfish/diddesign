*! _diddesign_prep.ado - Data preparation for DID estimation
*!
*! Prepares panel data or repeated cross-sectional data for difference-in-
*! differences estimation. Implements time normalization, group indicator
*! construction, and lagged outcome transformation for the double DID estimator.

version 16.0

program define _diddesign_prep, rclass
    
    // =========================================================================
    // SECTION 1: INPUT PARAMETERS
    // =========================================================================
    
    syntax, ///
        OUTCOME(varname)        /// Outcome variable
        TREATment(varname)      /// Treatment indicator
        TIME(varname)           /// Time variable
        [ID(varname)]           /// Unit ID (panel only)
        [POST(varname)]         /// Post indicator (RCS only)
        [CLuster(varname)]      /// Cluster variable
        [COVariates(varlist)]   /// Covariates (already expanded by caller)
        [TOUSE(varname)]        /// Sample marker
        [PANEL]                 /// Panel data flag
        [RCS]                   /// RCS data flag
    
    // Determine data type
    local is_panel = ("`panel'" != "")
    local is_rcs = ("`rcs'" != "")
    
    // Default to panel if neither specified but id is provided
    if !`is_panel' & !`is_rcs' {
        if "`id'" != "" {
            local is_panel = 1
        }
        else {
            display as error "E003: Must specify panel or rcs data type"
            exit 198
        }
    }
    
    // =========================================================================
    // SECTION 2: DATA VALIDATION
    // =========================================================================
    // Validates data structure requirements for DID estimation:
    // - Binary treatment indicator (0/1)
    // - Presence of both treatment and control observations
    // - For RCS: binary post-treatment indicator
    
    // Initialize sample marker if not provided
    if "`touse'" == "" {
        tempvar touse
        gen byte `touse' = 1
    }
    
    // Count observations
    quietly count if `touse'
    local N = r(N)
    if `N' == 0 {
        display as error "E003: No observations in sample"
        exit 2000
    }
    
    // Validate treatment is strictly binary (0/1)
    quietly tabulate `treatment' if `touse'
    if r(r) > 2 {
        display as error "E003: Treatment variable must be binary (0/1)"
        display as error "      Found `r(r)' distinct values"
        exit 198
    }
    quietly summarize `treatment' if `touse'
    // Handle all-missing treatment variable
    if r(N) == 0 | missing(r(min)) {
        display as error "E003: Treatment variable contains only missing values"
        exit 198
    }
    // Use floating-point tolerance for robust binary detection
    if round(r(min), 1e-6) != 0 | round(r(max), 1e-6) != 1 {
        display as error "E003: Treatment variable must contain both 0 and 1 values"
        display as error "      Found min=`r(min)', max=`r(max)'"
        exit 198
    }
    
    // Check for treated observations (preliminary)
    // Definitive check is performed after Gi computation for panel data
    quietly count if abs(`treatment' - 1) < 1e-6 & `touse'
    if r(N) == 0 {
        display as error "E003: No treated units found in data"
        exit 198
    }
    
    // Check for control observations (preliminary)
    // Definitive check is performed after Gi computation for panel data
    quietly count if abs(`treatment') < 1e-6 & `touse'
    if r(N) == 0 {
        display as error "E003: No control units found in data"
        exit 198
    }
    
    // For RCS: validate post indicator is binary
    if `is_rcs' {
        quietly summarize `post' if `touse'
        
        // Handle all-missing post variable
        if r(N) == 0 | missing(r(min)) | missing(r(max)) {
            display as error "E003: Post-treatment indicator contains only missing values"
            exit 198
        }
        
        // Check post indicator contains at most 2 distinct values
        quietly tab `post' if `touse'
        if r(r) > 2 {
            display as error "E003: Post-treatment indicator must be binary (0/1)"
            display as error "       Found " r(r) " distinct values (expected 2)"
            exit 198
        }
        
        // Verify values are exactly 0 and 1 (use floating-point tolerance)
        quietly summarize `post' if `touse'
        if (round(r(min), 1e-6) != 0 & round(r(min), 1e-6) != 1) | (round(r(max), 1e-6) != 0 & round(r(max), 1e-6) != 1) {
            display as error "E003: Post-treatment indicator must contain exactly 0 and 1"
            display as error "       Found min=" r(min) ", max=" r(max)
            exit 198
        }
        
        // Ensure both pre- and post-treatment observations exist
        if r(min) == r(max) {
            if r(min) == 0 {
                display as error "E003: Post-treatment indicator is all 0 (no post-treatment observations)"
            }
            else {
                display as error "E003: Post-treatment indicator is all 1 (no pre-treatment observations)"
            }
            exit 198
        }
    }

    // =========================================================================
    // SECTION 3: PANEL DATA PREPARATION
    // =========================================================================
    // Prepares panel data by:
    // 1. Converting time variable to consecutive integers
    // 2. Identifying treatment time (first period where treatment == 1)
    // 3. Creating group indicator Gi (1 if unit ever treated)
    // 4. Creating post-treatment indicator It
    // 5. Creating standardized time index (centered at treatment time)
    
    if `is_panel' {
        
        // -----------------------------------------------------------------
        // Step 1: Time index conversion
        // Convert time variable to consecutive integers for standardization
        // -----------------------------------------------------------------
        tempvar id_time_n
        egen `id_time_n' = group(`time') if `touse'
        
        // -----------------------------------------------------------------
        // Step 2: Identify treatment time
        // Treatment time is the minimum time period where treatment == 1
        // -----------------------------------------------------------------
        quietly summarize `id_time_n' if abs(`treatment' - 1) < 1e-6 & `touse', meanonly
        if r(N) == 0 {
            display as error "E003: Cannot identify treatment time"
            exit 198
        }
        local treat_year = r(min)
        
        // -----------------------------------------------------------------
        // Step 3: Create group indicator (Gi)
        // Gi = 1 if unit is ever treated, 0 otherwise.
        // Rounding ensures Gi is exactly 0 or 1 for downstream comparisons.
        // -----------------------------------------------------------------
        tempvar Gi_temp
        egen `Gi_temp' = max(`treatment') if `touse', by(`id')
        replace `Gi_temp' = round(`Gi_temp', 1) if `touse'
        
        // -----------------------------------------------------------------
        // Step 4: Create post-treatment indicator (It)
        // It = 1 if time period >= treatment time, 0 otherwise
        // -----------------------------------------------------------------
        tempvar It_temp
        gen `It_temp' = (`id_time_n' >= `treat_year') if `touse'
        
        // -----------------------------------------------------------------
        // Step 5: Create standardized time index
        // Centered at treatment time (treatment time = 0)
        // -----------------------------------------------------------------
        tempvar id_time_std_temp
        gen `id_time_std_temp' = `id_time_n' - `treat_year' if `touse'
        
        // -----------------------------------------------------------------
        // Step 6: Validate control units exist
        // After Gi is computed, units with Gi = 0 must exist
        // -----------------------------------------------------------------
        quietly count if `Gi_temp' == 0 & `touse'
        if r(N) == 0 {
            display as error "E003: No control units found in data (all units eventually treated)"
            exit 198
        }
        
        // -----------------------------------------------------------------
        // Step 7: Count units and periods
        // -----------------------------------------------------------------
        tempvar unit_tag
        egen `unit_tag' = tag(`id') if `touse'
        quietly count if `unit_tag' == 1 & `touse'
        local n_units = r(N)
        
        quietly summarize `id_time_n' if `touse'
        local n_periods = r(max) - r(min) + 1
        
        // Store in permanent variables
        capture drop _did_id_time
        capture drop _did_id_time_std
        capture drop _did_Gi
        capture drop _did_It
        
        gen _did_id_time = `id_time_n' if `touse'
        gen _did_id_time_std = `id_time_std_temp' if `touse'
        gen _did_Gi = `Gi_temp' if `touse'
        gen _did_It = `It_temp' if `touse'
        
        local id_var = "`id'"
    }

    // =========================================================================
    // SECTION 4: REPEATED CROSS-SECTIONAL DATA PREPARATION
    // =========================================================================
    // Prepares repeated cross-sectional (RCS) data by:
    // 1. Assigning Gi directly from treatment variable
    // 2. Assigning It directly from post indicator variable
    // 3. Converting time variable to consecutive integers
    // 4. Creating standardized time index (centered at treatment time)
    //
    // Key difference from panel: Gi and It are directly observed, not derived
    
    else if `is_rcs' {
        
        // -----------------------------------------------------------------
        // Step 1-2: Assign Gi and It directly
        // For RCS data, Gi (treatment group) and It (post period) are
        // directly observed rather than derived from panel structure.
        // -----------------------------------------------------------------
        tempvar Gi_temp It_temp
        gen `Gi_temp' = `treatment' if `touse'
        gen `It_temp' = `post' if `touse'
        
        // -----------------------------------------------------------------
        // Step 3: Time index conversion
        // Convert time variable to consecutive integers for standardization
        // -----------------------------------------------------------------
        tempvar id_time_n
        egen `id_time_n' = group(`time') if `touse'
        
        // -----------------------------------------------------------------
        // Step 4: Identify treatment time
        // For RCS, treatment time is the minimum time period where It == 1
        // -----------------------------------------------------------------
        quietly summarize `id_time_n' if `It_temp' == 1 & `touse', meanonly
        if r(N) == 0 {
            display as error "E003: Cannot identify treatment time from post indicator"
            exit 198
        }
        local treat_year = r(min)
        
        // -----------------------------------------------------------------
        // Step 5: Create standardized time index
        // Centered at treatment time (treatment time = 0)
        // -----------------------------------------------------------------
        tempvar id_time_std_temp
        gen `id_time_std_temp' = `id_time_n' - `treat_year' if `touse'
        
        // -----------------------------------------------------------------
        // Step 6: Count periods
        // Unit count is not applicable for RCS data
        // -----------------------------------------------------------------
        local n_units = .  // Not applicable for RCS
        
        quietly summarize `id_time_n' if `touse'
        local n_periods = r(max) - r(min) + 1
        
        // Store in permanent variables
        capture drop _did_id_time
        capture drop _did_id_time_std
        capture drop _did_Gi
        capture drop _did_It
        
        gen _did_id_time = `id_time_n' if `touse'
        gen _did_id_time_std = `id_time_std_temp' if `touse'
        gen _did_Gi = `Gi_temp' if `touse'
        gen _did_It = `It_temp' if `touse'
        
        local id_var = ""  // No unit ID for RCS
    }

    // =========================================================================
    // SECTION 5: LAGGED OUTCOME TRANSFORMATION
    // =========================================================================
    // Computes the lagged group mean transformation for the sequential DID
    // estimator, which is consistent under the parallel trends-in-trends
    // assumption:
    //
    //   ΔY_{it} = Y_{it} - Ȳ_{Gi,t-1}
    //
    // where Ȳ_{Gi,t-1} is the mean outcome for group Gi at standardized time t-1.
    // This transformation enables the double DID estimator to combine standard
    // DID and sequential DID via GMM for optimal efficiency.
    //
    // Algorithm:
    // 1. Compute mean outcome for each (Gi, id_time_std) combination
    // 2. Shift time index by +1 to create the lag structure
    // 3. Merge lagged means back to original data
    // 4. Compute outcome_delta = outcome - lagged_group_mean
    //
    // The earliest period has missing outcome_delta (no prior period to lag from)
    // =========================================================================
    
    // -----------------------------------------------------------------
    // Step 1: Compute group-period means
    // Mean outcome is computed for each (Gi, id_time_std) combination.
    // Groups containing any missing values are assigned missing mean.
    // -----------------------------------------------------------------
    
    // Use preserve/restore to create the collapsed dataset
    preserve
    
    // Keep only needed variables and sample (quietly to suppress output)
    quietly keep if `touse'
    quietly keep `outcome' _did_Gi _did_id_time_std
    
    // Track groups with missing values
    quietly gen byte _has_na = missing(`outcome')
    
    // Compute means by (Gi, id_time_std)
    // (max) _has_na captures whether any observation in the group has missing value
    quietly collapse (mean) _did_Ymean = `outcome' (max) _has_na, by(_did_Gi _did_id_time_std)
    
    // Propagate missing: if group has any missing, set mean to missing
    quietly replace _did_Ymean = . if _has_na == 1
    quietly drop _has_na
    
    // -----------------------------------------------------------------
    // Step 2: Time-shift by one period
    // Shifting id_time_std by +1 creates the lagged structure.
    // After merge, each observation at time t receives the mean from t-1.
    // -----------------------------------------------------------------
    quietly replace _did_id_time_std = _did_id_time_std + 1
    
    // Save lagged means to tempfile
    tempfile lagged_means
    quietly save `lagged_means', replace
    
    restore
    
    // -----------------------------------------------------------------
    // Step 3: Merge lagged means
    // Left join preserves all observations. Unmatched observations
    // (earliest period) receive missing lagged mean.
    // -----------------------------------------------------------------
    quietly merge m:1 _did_Gi _did_id_time_std using `lagged_means', ///
        keep(master match)
    
    // -----------------------------------------------------------------
    // Validate merge results
    // -----------------------------------------------------------------
    quietly count if _merge == 1  // master only (no match)
    local n_nomatch = r(N)
    
    quietly count if _merge == 3  // matched
    local n_matched = r(N)
    
    if `n_matched' == 0 {
        display as error "E011: Merge failed - no observations matched lagged means"
        display as error "      This may indicate data structure issues"
        quietly drop _merge
        exit 198
    }
    
    // Warning for high non-match rate (expected for earliest period, but warn if excessive)
    quietly count if `touse'
    local n_total_merge = r(N)
    local pct_nomatch = 100 * `n_nomatch' / `n_total_merge'
    if `pct_nomatch' > 50 {
        display as text "Warning: `pct_nomatch'% of observations had no matching lagged mean"
    }
    
    quietly drop _merge
    
    // -----------------------------------------------------------------
    // Step 4: Compute outcome delta
    // ΔY_{it} = Y_{it} - Ȳ_{Gi,t-1}
    // Missing lagged mean produces missing outcome_delta (earliest period).
    // -----------------------------------------------------------------
    capture drop _did_outcome_delta
    quietly gen _did_outcome_delta = `outcome' - _did_Ymean if `touse'
    
    // Clean up temporary Ymean variable (not needed after delta calculation)
    capture drop _did_Ymean
    
    // -----------------------------------------------------------------
    // Step 5: Validate outcome delta
    // -----------------------------------------------------------------
    
    // Count missing outcome_delta
    quietly count if missing(_did_outcome_delta) & `touse'
    local n_missing_delta = r(N)
    
    // Count observations in earliest period (expected to have missing delta)
    quietly summarize _did_id_time_std if `touse', meanonly
    if r(N) == 0 | missing(r(min)) {
        display as text "Warning: Cannot determine earliest time period (all id_time_std missing)"
        local min_time_std = .
        local n_earliest = 0
    }
    else {
        local min_time_std = r(min)
        quietly count if _did_id_time_std == `min_time_std' & `touse'
        local n_earliest = r(N)
    }
    
    // Warning if more observations have missing delta than expected
    quietly count if `touse'
    local n_total = r(N)
    
    if `n_missing_delta' > `n_earliest' {
        local extra_missing = `n_missing_delta' - `n_earliest'
        display as text "Warning: `extra_missing' additional observations have missing outcome_delta"
        display as text "         (beyond the `n_earliest' expected for earliest time period)"
        display as text "         This may indicate missing outcome values in the data"
    }
    
    // General warning if high percentage missing
    // For RCS with few time periods (2-3), expect ~33-50% missing as normal
    local pct_missing = 100 * `n_missing_delta' / `n_total'
    if `pct_missing' > 30 & `is_rcs' & `n_periods' <= 3 {
        // This is expected for RCS with few time periods - suppress warning
        // The earliest period has no lagged data to compute outcome_delta
    }
    else if `pct_missing' > 30 {
        display as text "Warning: `pct_missing'% of observations have missing outcome_delta"
    }

    // =========================================================================
    // SECTION 6: MATA STRUCTURE POPULATION
    // =========================================================================
    // Transfer prepared data to Mata structure for estimation
    
    // Store scalar values for Mata
    local treat_year_std = 0  // Standardized treatment time is defined as zero
    
    // Call Mata function to populate structure
    mata: st_local("mata_rc", strofreal(_diddesign_populate_data( ///
        "`outcome'",           /* outcome variable name     */ ///
        "`treatment'",         /* treatment variable name   */ ///
        "`id_var'",            /* id variable name (or "")  */ ///
        "_did_id_time",        /* id_time variable name     */ ///
        "`covariates'",        /* covariate variable names  */ ///
        "`cluster'",           /* cluster variable name     */ ///
        "_did_Gi",             /* Gi variable name          */ ///
        "_did_It",             /* It variable name          */ ///
        "_did_id_time_std",    /* id_time_std variable name */ ///
        "_did_outcome_delta",  /* outcome_delta var name    */ ///
        `N',                   /* N observations            */ ///
        `n_units',             /* n_units                   */ ///
        `n_periods',           /* n_periods                 */ ///
        `treat_year_std',      /* treat_year (standardized) */ ///
        `is_panel',            /* is_panel flag             */ ///
        "`touse'"              /* touse variable name       */ ///
    )))
    
    if `mata_rc' != 0 {
        display as error "Error populating Mata did_data structure"
        exit 498
    }
    
    // =========================================================================
    // SECTION 7: RETURN VALUES
    // =========================================================================
    // Return prepared variable names and computed scalars to caller
    
    // Variable names (created variables)
    return local id_time = "_did_id_time"
    return local id_time_std = "_did_id_time_std"
    return local Gi = "_did_Gi"
    return local It = "_did_It"
    return local outcome_delta = "_did_outcome_delta"
    
    // Scalar values
    return scalar N = `N'
    return scalar n_units = `n_units'
    return scalar n_periods = `n_periods'
    return scalar treat_year = `treat_year'
    return scalar treat_year_std = `treat_year_std'
    return scalar is_panel = `is_panel'
    return scalar n_missing_delta = `n_missing_delta'

end
