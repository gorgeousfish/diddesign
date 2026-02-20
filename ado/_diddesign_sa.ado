*! _diddesign_sa.ado - Staggered adoption design estimation
*!
*! Implements the staggered adoption (SA) extension of the double DID estimator
*! for settings where treatment timing varies across units. The SA design
*! estimates period-specific ATT at each adoption time t using the double DID
*! framework, then aggregates via time-weighted average:
*!
*!   tau_bar^SA = Sum_t pi_t * tau^SA(t)
*!
*! where pi_t is the proportion of newly treated units at period t. This module
*! serves as the Stata interface for SA estimation, delegating numerical
*! computation to the Mata functions in did_sa.mata.

program define _diddesign_sa, eclass
    version 16.0
    
    // =========================================================================
    // SECTION 1: SYNTAX PARSING
    // =========================================================================
    syntax varlist(min=1 fv) [if] [in], ///
        TREATment(varname)              /// Required: treatment indicator
        ID(varname)                     /// Required for SA: unit identifier
        TIME(varname)                   /// Required for SA: time identifier
        [CLuster(varname)]              /// Cluster variable for SEs
        [COVariates(string)]            /// Additional covariates (supports factor variables)
        [NBoot(integer 30)]             /// Bootstrap iterations (default: 30)
        [LEAD(numlist >=0 integer)]     /// Lead values for SA design
        [THRes(integer 2)]              /// SA threshold (default: 2)
        [LEVEL(cilevel)]                /// Confidence level (default: 95)
        [SEED(integer -1)]              /// Random seed (-1 = not specified)
        [QUIET]                         /// Suppress progress display
        TOUSE(varname)                  /// Sample marker from main program
        [IDORIG(string)]                /// Original id variable name for ereturn
        [TIMEORIG(string)]              /// Original time variable name for ereturn
        [CLUSTERORIG(string)]           /// Original cluster variable name for ereturn
    
    // Read command line from global macro
    local cmdline "$DIDDESIGN_CMDLINE"
    
    // =========================================================================
    // SECTION 2: SET DEFAULTS
    // =========================================================================
    // Parse outcome from varlist
    gettoken outcome rest : varlist
    local covariates_inline = "`rest'"
    local all_covariates = "`covariates_inline' `covariates'"
    local all_covariates = strtrim("`all_covariates'")
    
    // -------------------------------------------------------------------------
    // Duplicate Covariate Check
    // -------------------------------------------------------------------------
    // Remove duplicate covariates when combining inline and covariates() option
    
    if "`all_covariates'" != "" {
        local unique_covars : list uniq all_covariates
        local n_all : word count `all_covariates'
        local n_unique : word count `unique_covars'
        
        if `n_unique' < `n_all' {
            // Find duplicate variables by comparing original and unique lists
            local dups ""
            local seen ""
            foreach v of local all_covariates {
                local is_seen : list v in seen
                if `is_seen' {
                    local is_dup : list v in dups
                    if !`is_dup' {
                        local dups "`dups' `v'"
                    }
                }
                else {
                    local seen "`seen' `v'"
                }
            }
            local dups = strtrim("`dups'")
            display as text "Warning: Duplicate covariates detected and removed: `dups'"
            local all_covariates "`unique_covars'"
        }
    }
    
    // -------------------------------------------------------------------------
    // Factor Variable Expansion
    // -------------------------------------------------------------------------
    // Expand factor variables (i.var, ibn.var) into dummy variables
    // Base category is excluded to avoid collinearity with the intercept
    
    if "`all_covariates'" != "" {
        // Check if covariates contain factor variable notation
        fvexpand `all_covariates' if `touse'
        local has_fvops "`r(fvops)'"
        
        if "`has_fvops'" == "true" {
            // Process each covariate term separately
            // For factor variables: expand and exclude base category (first dummy)
            local expanded_vars ""
            local n_fv_expanded = 0
            
            foreach covar_term of local all_covariates {
                // Check if this term is a factor variable expression
                local is_factor = 0
                if regexm("`covar_term'", "^i\.") | regexm("`covar_term'", "^i\(") {
                    local is_factor = 1
                }
                // Match factor variable notations including multi-digit base categories
                if regexm("`covar_term'", "^ibn\.") | regexm("`covar_term'", "^ib[0-9]+") {
                    local is_factor = 1
                }
                if regexm("`covar_term'", "^c\.") {
                    local is_factor = 0
                }
                
                if `is_factor' {
                    fvrevar `covar_term' if `touse'
                    local fv_vars "`r(varlist)'"
                    local n_fv : word count `fv_vars'
                    
                    if `n_fv' > 1 {
                        // Exclude first variable (base category)
                        forvalues i = 2/`n_fv' {
                            local v : word `i' of `fv_vars'
                            local expanded_vars "`expanded_vars' `v'"
                        }
                        local n_fv_expanded = `n_fv_expanded' + `n_fv' - 1
                    }
                    else if `n_fv' == 1 {
                        local expanded_vars "`expanded_vars' `fv_vars'"
                        local n_fv_expanded = `n_fv_expanded' + 1
                    }
                }
                else {
                    capture fvrevar `covar_term' if `touse'
                    if _rc == 0 {
                        local expanded_vars "`expanded_vars' `r(varlist)'"
                    }
                    else {
                        local expanded_vars "`expanded_vars' `covar_term'"
                    }
                }
            }
            
            local all_covariates "`expanded_vars'"
            local all_covariates : list retokenize all_covariates
            
            if `n_fv_expanded' > 0 {
                display as text "Note: Factor variables expanded to `n_fv_expanded' dummy variables (base categories excluded)"
            }
        }
    }
    
    // Set default values
    if "`lead'" == "" local lead = "0"
    local nboot_val = `nboot'
    local thres_val = `thres'
    local level_val = `level'
    
    // Set default original variable names if not provided
    if "`idorig'" == "" local idorig "`id'"
    if "`timeorig'" == "" local timeorig "`time'"
    if "`clusterorig'" == "" local clusterorig "`cluster'"
    
    // Cluster defaults to id if not specified
    if "`cluster'" == "" {
        local cluster_var = "`id'"
        // Also set clusterorig to idorig if cluster was not specified
        if "`clusterorig'" == "" local clusterorig "`idorig'"
    }
    else {
        local cluster_var = "`cluster'"
    }
    
    // Seed handling
    if `seed' == -1 {
        local seed_val = .
    }
    else {
        local seed_val = `seed'
    }
    
    // Quiet option
    local quiet_val = 0
    if "`quiet'" != "" {
        local quiet_val = 1
    }
    
    // =========================================================================
    // SECTION 3: DATA PREPARATION
    // =========================================================================
    // SA design requires balanced panel structure where each unit is observed
    // across all time periods. The treatment timing matrix G_{it} classifies
    // each unit-period as: newly treated (1), not-yet-treated control (0), or
    // previously treated (-1). Valid periods must have at least 'thres' newly
    // treated units to ensure reliable estimation.
    // -------------------------------------------------------------------------
    
    // Count observations and units
    quietly count if `touse'
    local N = r(N)
    
    // Count unique units
    tempvar unit_tag
    quietly egen `unit_tag' = tag(`id') if `touse'
    quietly count if `unit_tag' == 1 & `touse'
    local n_units = r(N)
    
    // Count unique time periods
    tempvar time_tag
    quietly egen `time_tag' = tag(`time') if `touse'
    quietly count if `time_tag' == 1 & `touse'
    local n_periods = r(N)
    
    // SA design requires at least 2 time periods for valid comparison
    if `n_periods' < 2 {
        display as error "E008: SA design requires at least 2 time periods"
        display as error "      Found only `n_periods' time period(s)"
        display as error "      SA design needs pre-treatment and post-treatment periods for comparison"
        exit 198
    }
    
    // =========================================================================
    // SECTION 4: VALIDATE TREATMENT VARIABLE
    // =========================================================================
    // SA design requires an absorbing (cumulative) binary treatment indicator:
    //   - Binary: D_{it} in {0, 1} for all observations
    //   - Absorbing: once treated, units remain treated (D_{it} = 1 => D_{is} = 1
    //     for all s > t)
    // This structure enables identification of treatment adoption timing A_i,
    // defined as the first period where D_{it} = 1. The treatment timing matrix
    // G_{it} is then constructed based on A_i to classify unit-period cells.
    // -------------------------------------------------------------------------
    
    // Check at most 2 distinct values
    quietly tabulate `treatment' if `touse'
    if r(r) > 2 {
        display as error "E007: Treatment variable must be binary (0/1)"
        display as error "      Found `r(r)' distinct values"
        display as error "      SA design requires cumulative binary treatment indicator"
        exit 459
    }
    
    // Check for missing values and valid range
    quietly summarize `treatment' if `touse'
    if r(N) == 0 | missing(r(min)) {
        display as error "E007: Treatment variable contains only missing values"
        exit 459
    }
    
    // Tolerance 1e-6 for floating-point comparison to handle numeric imprecision
    if round(r(min), 1e-6) != 0 | round(r(max), 1e-6) != 1 {
        display as error "E007: Treatment variable must contain both 0 and 1 values"
        display as error "      Found min=`r(min)', max=`r(max)'"
        display as error "      SA design requires cumulative binary treatment indicator"
        exit 459
    }
    
    // Check cumulative treatment (absorbing treatment)
    // Treatment must only transition from 0 to 1, never decrease
    tempvar treat_lag treat_diff
    quietly {
        bysort `id' (`time'): gen `treat_lag' = `treatment'[_n-1] if `touse'
        gen `treat_diff' = `treatment' - `treat_lag' if `touse' & `treat_lag' != .
        count if `treat_diff' < 0 & `touse'
    }
    if r(N) > 0 {
        display as error "E003: Treatment variable must be cumulative (absorbing)"
        display as error "      Found `r(N)' observations with treatment decreasing over time"
        display as error "      SA design requires treatment to only transition from 0 to 1"
        exit 459
    }
    
    // =========================================================================
    // SECTION 5: CALL MATA SA ESTIMATION
    // =========================================================================
    // Time-average SA-ATT (Staggered Adoption ATT):
    //
    //   tau_bar^SA = Sum_{t in T} pi_t * tau^SA(t)
    //
    // where:
    //   - pi_t = n_{1t} / Sum_{t'} n_{1t'} is the time weight (proportion of
    //     newly treated units at period t)
    //   - tau^SA(t) is the period-specific double DID estimate combining
    //     tau^SA_DID(t) and tau^SA_sDID(t) via GMM optimal weighting
    //   - Estimation uses three consecutive periods {t-2, t-1, t} per period t
    //
    // GMM weight matrix W = Omega^{-1} minimizes variance under heteroskedasticity.
    // Bootstrap variance estimation resamples units (not observations) with
    // replacement, recomputing all period-specific estimates in each iteration.
    // -------------------------------------------------------------------------
    
    // Set random seed if specified for bootstrap reproducibility
    if `seed_val' != . {
        set seed `seed_val'
    }
    
    // Convert lead numlist to Mata format
    local lead_mata = subinstr("`lead'", " ", ", ", .)
    local n_lead : word count `lead'
    
    // Prepare data in Mata
    mata: st_local("mata_rc", strofreal(_did_sa_prepare_data("`outcome'", "`treatment'", "`id'", "`time'", ///
                               "`cluster_var'", "`all_covariates'", "`touse'")))
    
    if `mata_rc' != 0 {
        display as error "E011: SA data preparation failed"
        display as error "      No valid observations selected for analysis"
        exit 498
    }
    
    // Call SA estimation
    mata: st_local("mata_rc", strofreal(_did_sa_main((`lead_mata'), `nboot_val', `thres_val', `level_val', `quiet_val')))
    
    if `mata_rc' != 0 {
        // Provide specific error messages based on error code
        if `mata_rc' == 1 {
            display as error "E011: SA estimation failed - could not create treatment timing matrix (Gmat)"
        }
        else if `mata_rc' == 2 {
            display as error "E011: SA estimation failed - no valid periods found"
            display as error "      Try reducing the threshold value (thres option)"
        }
        else if `mata_rc' == 3 {
            display as error "E011: SA estimation failed - point estimation returned missing values"
        }
        else if `mata_rc' == 4 {
            display as error "E011: SA estimation failed - insufficient valid bootstrap iterations"
            display as error "      Try increasing the number of bootstrap iterations (nboot option)"
        }
        else {
            display as error "E011: SA estimation failed in Mata (error code: `mata_rc')"
        }
        exit 498
    }
    
    // =========================================================================
    // SECTION 6: RETRIEVE RESULTS FROM MATA
    // =========================================================================
    // Transfer estimation metadata from Mata global scalars to Stata locals
    mata: st_local("n_periods_valid", strofreal(_sa_n_periods_valid))
    mata: st_local("n_boot_success", strofreal(_sa_n_boot_success))
    
    // =========================================================================
    // SECTION 7: STORE e() RETURNS
    // =========================================================================
    // Transfer estimation results from Mata to Stata e() class for post-estimation
    // commands. Results include: coefficient vector (b), variance matrix (V),
    // detailed estimates table, GMM weight matrix (W), and time weights (pi_t).
    // -------------------------------------------------------------------------
    
    tempname b_mat V_mat estimates_mat lead_mat weights_mat W_mat vcov_gmm_mat time_weights_mat
    
    mata: st_matrix("`b_mat'", _sa_b)
    mata: st_matrix("`V_mat'", _sa_V)
    mata: st_matrix("`estimates_mat'", _sa_estimates)
    mata: st_matrix("`lead_mat'", _sa_lead_values)
    mata: st_matrix("`weights_mat'", _sa_weights)
    mata: st_matrix("`W_mat'", _sa_W)
    mata: st_matrix("`vcov_gmm_mat'", _sa_vcov_gmm)
    mata: st_matrix("`time_weights_mat'", _sa_time_weights)
    
    // Validate result matrices exist and are non-empty
    capture confirm matrix `lead_mat'
    if _rc != 0 {
        display as error "Error: SA estimation produced no valid results (lead_mat not found)"
        exit 498
    }
    if colsof(`lead_mat') == 0 {
        display as error "Error: SA estimation produced no valid results (lead_mat is empty)"
        exit 498
    }
    
    // Reshape flattened W and VCOV matrices to proper 2x2 form for single lead case
    // Mata vec() uses column-major order: [W11, W21, W12, W22] for a 2x2 matrix
    // For multiple leads, matrices remain as n_lead x 4 (each row is one flattened 2x2)
    local n_lead = colsof(`lead_mat')
    if `n_lead' == 1 {
        // Reconstruct 2x2 GMM weight matrix W = Omega^{-1}
        matrix `W_mat' = (`W_mat'[1,1], `W_mat'[1,3] \ `W_mat'[1,2], `W_mat'[1,4])
        
        // Reconstruct 2x2 variance-covariance matrix Omega of moment conditions
        matrix `vcov_gmm_mat' = (`vcov_gmm_mat'[1,1], `vcov_gmm_mat'[1,3] \ `vcov_gmm_mat'[1,2], `vcov_gmm_mat'[1,4])
    }
    
    // Set row and column names for e(b)
    local b_names ""
    foreach l of numlist `lead' {
        local b_names "`b_names' SA_dDID:lead_`l' SA_DID:lead_`l' SA_sDID:lead_`l'"
    }
    local b_names = trim("`b_names'")
    matrix colnames `b_mat' = `b_names'
    
    // Set row and column names for e(V)
    matrix rownames `V_mat' = `b_names'
    matrix colnames `V_mat' = `b_names'
    
    // Post b and V matrices with sample marker
    ereturn post `b_mat' `V_mat', esample(`touse') obs(`N') depname("`outcome'")
    
    // --- Scalars ---
    ereturn scalar n_units = `n_units'
    ereturn scalar n_periods = `n_periods'
    ereturn scalar n_periods_valid = `n_periods_valid'
    ereturn scalar n_boot = `nboot_val'
    ereturn scalar level = `level_val'
    ereturn scalar n_lead = `n_lead'
    ereturn scalar thres = `thres_val'
    ereturn scalar is_panel = 1
    
    // n_boot_success only if some iterations failed
    if "`n_boot_success'" != "" & "`n_boot_success'" != "." {
        if `n_boot_success' < `nboot_val' {
            ereturn scalar n_boot_success = `n_boot_success'
        }
    }
    
    // --- Macros ---
    ereturn local cmd "diddesign"
    ereturn local cmdline "`cmdline'"
    ereturn local design "sa"
    ereturn local depvar "`outcome'"
    ereturn local treatment "`treatment'"
    ereturn local covariates "`all_covariates'"
    ereturn local id "`idorig'"
    ereturn local time "`timeorig'"
    ereturn local clustvar "`clusterorig'"
    ereturn local ci_method "bootstrap"
    ereturn local lead "`lead'"
    ereturn local properties "b V"
    
    // --- Additional Matrices ---
    // Set row and column names for e(estimates)
    local est_rownames ""
    foreach l of numlist `lead' {
        local est_rownames "`est_rownames' SA_dDID:lead_`l' SA_DID:lead_`l' SA_sDID:lead_`l'"
    }
    local est_rownames = trim("`est_rownames'")
    matrix rownames `estimates_mat' = `est_rownames'
    matrix colnames `estimates_mat' = lead estimate std_error ci_lo ci_hi weight
    
    // Set names for e(lead_values)
    matrix colnames `lead_mat' = `lead'
    
    // Set names for e(weights)
    local wt_rownames ""
    foreach l of numlist `lead' {
        local wt_rownames "`wt_rownames' lead_`l'"
    }
    local wt_rownames = trim("`wt_rownames'")
    matrix rownames `weights_mat' = `wt_rownames'
    matrix colnames `weights_mat' = w_did w_sdid
    
    // Make a copy of estimates_mat for display
    tempname display_mat
    matrix `display_mat' = `estimates_mat'
    
    // Store additional matrices
    ereturn matrix estimates = `estimates_mat'
    ereturn matrix lead_values = `lead_mat'
    ereturn matrix weights = `weights_mat'
    ereturn matrix W = `W_mat'
    ereturn matrix vcov_gmm = `vcov_gmm_mat'
    
    // Set row/column names for time_weights matrix
    local n_tw = rowsof(`time_weights_mat')
    if `n_tw' > 0 {
        local tw_rownames ""
        forvalues i = 1/`n_tw' {
            local tw_rownames "`tw_rownames' period_`i'"
        }
        local tw_rownames = trim("`tw_rownames'")
        matrix rownames `time_weights_mat' = `tw_rownames'
        matrix colnames `time_weights_mat' = weight
    }
    ereturn matrix time_weights = `time_weights_mat'
    
    // =========================================================================
    // SECTION 8: DISPLAY RESULTS
    // =========================================================================
    _diddesign_display_header, cmd("diddesign") design("sa") ///
        datatype("Panel") n(`N') n_units(`n_units') ///
        n_periods(`n_periods') n_boot(`nboot_val') ///
        cluster("`clusterorig'") thres(`thres_val')
    
    // Abbreviate depvar for table header
    local depvar_abbrev = abbrev("`outcome'", 13)
    
    // Display results table for each lead
    local row = 1
    foreach l of numlist `lead' {
        // Table header with depvar~lead label
        if `n_lead' > 1 {
            local tbl_label "`depvar_abbrev'~`l'"
        }
        else {
            local tbl_label "`depvar_abbrev'"
        }
        
        _diddesign_display_table_header, depvar("`tbl_label'") level(`level_val') weight
        
        // SA-Double-DID
        local est = `display_mat'[`row', 2]
        local se = `display_mat'[`row', 3]
        local ci_lo = `display_mat'[`row', 4]
        local ci_hi = `display_mat'[`row', 5]
        local wt = `display_mat'[`row', 6]
        _diddesign_display_result, label("SA-Double DID") ///
            estimate(`est') se(`se') ci_low(`ci_lo') ci_high(`ci_hi') weight(`wt')
        local row = `row' + 1
        
        // SA-DID
        local est = `display_mat'[`row', 2]
        local se = `display_mat'[`row', 3]
        local ci_lo = `display_mat'[`row', 4]
        local ci_hi = `display_mat'[`row', 5]
        local wt = `display_mat'[`row', 6]
        _diddesign_display_result, label("SA-DID") ///
            estimate(`est') se(`se') ci_low(`ci_lo') ci_high(`ci_hi') weight(`wt')
        local row = `row' + 1
        
        // SA-sDID
        local est = `display_mat'[`row', 2]
        local se = `display_mat'[`row', 3]
        local ci_lo = `display_mat'[`row', 4]
        local ci_hi = `display_mat'[`row', 5]
        local wt = `display_mat'[`row', 6]
        _diddesign_display_result, label("SA-sDID") ///
            estimate(`est') se(`se') ci_low(`ci_lo') ci_high(`ci_hi') weight(`wt')
        local row = `row' + 1
    }
    
    display as text "{hline 78}"
    
    // Display notes
    display as text ""
    display as text "Note: Double DID combines DID and sDID via optimal GMM weights."
    
end
