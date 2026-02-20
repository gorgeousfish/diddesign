*! diddesign.ado - Main estimation command for Double DID
*!
*! Implements the Double Difference-in-Differences estimator for standard
*! DID designs with multiple pre-treatment periods. Combines standard DID
*! and sequential DID estimators via GMM for optimal efficiency.

program define diddesign, eclass
    version 16.0
    
    // =========================================================================
    // SECTION 0: Initialize Mata Library
    // =========================================================================
    // Mata functions are loaded if not already available in memory
    
    capture mata: mata describe _did_std_main()
    if _rc != 0 {
        // Mata functions are not loaded; attempt to locate and load them
        local mata_loaded = 0
        
        // Method 1: Direct findfile for diddesign_mata.do (works after net install)
        qui capture findfile diddesign_mata.do
        if _rc == 0 {
            quietly do "`r(fn)'"
            local mata_loaded = 1
        }
        
        // Method 2: Relative path from ado file (works in development environment)
        if !`mata_loaded' {
            qui capture findfile diddesign.ado
            if _rc == 0 {
                local ado_path = subinstr("`r(fn)'", char(92), "/", .)
                local ado_dir = reverse(substr(reverse("`ado_path'"), strpos(reverse("`ado_path'"), "/") + 1, .))
                local mata_path "`ado_dir'/../mata/diddesign_mata.do"
                capture confirm file "`mata_path'"
                if _rc == 0 {
                    quietly do "`mata_path'"
                    local mata_loaded = 1
                }
            }
        }
        
        // Verify loading succeeded
        if !`mata_loaded' {
            capture mata: mata describe _did_std_main()
            if _rc != 0 {
                display as error "E001: DIDdesign Mata functions not found"
                display as error "Mata files could not be located in adopath or relative paths."
                display as error "Solutions:"
                display as error "  1. Reinstall: net install diddesign, from(...) replace"
                display as error "  2. Or manually: do {path}/diddesign_mata.do"
                exit 198
            }
        }
    }
    
    // =========================================================================
    // SECTION 1: Parse Command Syntax
    // =========================================================================
    
    syntax varlist(min=1 fv) [if] [in], ///
        TREATment(varname)              /// Required: treatment indicator
        [ID(varname)]                   /// Unit identifier (required for panel)
        [TIME(varname)]                 /// Time identifier (required)
        [POST(varname)]                 /// Post-treatment indicator (RCS only)
        [CLuster(varname)]              /// Cluster variable for SEs
        [COVariates(string)]            /// Additional covariates (supports factor variables via string)
        [NBoot(integer 30)]             /// Bootstrap iterations (default: 30)
        [LEAD(numlist >=0 integer)]     /// Lead values for SA design
        [THRes(integer 2)]              /// SA threshold (default: 2)
        [LEVEL(cilevel)]                /// Confidence level (Stata default: 95)
        [SEED(integer -1)]              /// Random seed (-1 = not specified)
        [DESIGN(string)]                /// Design type: "did" (default) or "sa"
        [PARALlel]                      /// Use parallel computing
        [SEBoot]                        /// Use bootstrap SE/CI
        [PANEL]                         /// Panel data format
        [RCS]                           /// Repeated cross-section format
        [QUIET]                         /// Suppress progress display
    
    // Full command line is stored for e(cmdline)
    local cmdline "diddesign `0'"
    
    // =========================================================================
    // SECTION 2: Validate Parameters
    // =========================================================================
    // Detailed validation is delegated to _diddesign_parse
    
    local parse_opts "treatment(`treatment')"
    if "`id'" != "" {
        local parse_opts "`parse_opts' id(`id')"
    }
    if "`time'" != "" {
        local parse_opts "`parse_opts' time(`time')"
    }
    if "`post'" != "" {
        local parse_opts "`parse_opts' post(`post')"
    }
    if "`cluster'" != "" {
        local parse_opts "`parse_opts' cluster(`cluster')"
    }
    if "`covariates'" != "" {
        local parse_opts "`parse_opts' covariates(`covariates')"
    }
    local parse_opts "`parse_opts' nboot(`nboot')"
    if "`lead'" != "" {
        local parse_opts "`parse_opts' lead(`lead')"
    }
    local parse_opts "`parse_opts' thres(`thres')"
    local parse_opts "`parse_opts' level(`level')"
    local parse_opts "`parse_opts' seed(`seed')"
    if "`design'" != "" {
        local parse_opts "`parse_opts' design(`design')"
    }
    if "`parallel'" != "" {
        // Parallel bootstrap is reserved for future implementation
        display as text "{p 0 4 2}"
        display as text "Note: The {bf:parallel} option is currently not implemented. "
        display as text "Bootstrap iterations will run sequentially.{p_end}"
        local parse_opts "`parse_opts' parallel"
    }
    if "`seboot'" != "" {
        local parse_opts "`parse_opts' seboot"
    }
    if "`panel'" != "" {
        local parse_opts "`parse_opts' panel"
    }
    if "`rcs'" != "" {
        local parse_opts "`parse_opts' rcs"
    }
    
    _diddesign_parse `varlist' `if' `in', `parse_opts'
    
    // Parsed values are retrieved from r() before subsequent commands overwrite them
    local outcome "`r(outcome)'"
    local treatment_var "`r(treatment)'"
    local id_var "`r(id)'"
    local time_var "`r(time)'"
    local post_var "`r(post)'"
    local cluster_var "`r(cluster)'"
    local covariates_list "`r(covariates)'"
    
    // Scalar returns
    local nboot_val = r(nboot)
    local thres_val = r(thres)
    local level_val = r(level)
    local seed_val = r(seed)
    
    // String returns
    local lead_val "`r(lead)'"
    
    // Scalar returns
    local parallel_val = r(parallel)
    local seboot_val = r(seboot)
    local is_panel = r(is_panel)
    
    // String variable indicators for automatic encoding
    local id_is_string = r(id_is_string)
    local time_is_string = r(time_is_string)
    local cluster_is_string = r(cluster_is_string)
    
    // String return
    local design_val "`r(design)'"
    
    // -------------------------------------------------------------------------
    // String Variable Encoding
    // -------------------------------------------------------------------------
    // Encode string variables to numeric; tempvars must persist through estimation
    // Original variable names are preserved for reporting in e() returns
    
    local id_var_orig "`id_var'"
    local time_var_orig "`time_var'"
    local cluster_var_orig "`cluster_var'"
    
    if `id_is_string' == 1 & "`id_var'" != "" {
        tempvar id_encoded
        quietly egen `id_encoded' = group(`id_var')
        display as text "Note: String variable `id_var' automatically encoded to numeric"
        local id_var "`id_encoded'"
    }
    
    if `time_is_string' == 1 & "`time_var'" != "" {
        tempvar time_encoded
        quietly egen `time_encoded' = group(`time_var')
        display as text "Note: String variable `time_var' automatically encoded to numeric"
        local time_var "`time_encoded'"
    }
    
    if `cluster_is_string' == 1 & "`cluster_var'" != "" {
        tempvar cluster_encoded
        quietly egen `cluster_encoded' = group(`cluster_var')
        display as text "Note: String variable `cluster_var' automatically encoded to numeric"
        local cluster_var "`cluster_encoded'"
    }
    
    // Sample marker is created based on [if] [in] conditions
    marksample touse, novarlist
    
    // Observations with missing values in required variables are excluded
    markout `touse' `outcome'
    markout `touse' `treatment_var'
    if "`id_var'" != "" {
        markout `touse' `id_var'
    }
    if "`time_var'" != "" {
        markout `touse' `time_var'
    }
    if "`post_var'" != "" {
        markout `touse' `post_var'
    }
    if "`cluster_var'" != "" {
        markout `touse' `cluster_var'
    }
    
    // -------------------------------------------------------------------------
    // Factor Variable Expansion
    // -------------------------------------------------------------------------
    // Factor variables (i.var, ibn.var, etc.) are expanded into dummy variables
    // Base category is excluded to avoid collinearity with the intercept term
    
    if "`covariates_list'" != "" {
        // Covariates are checked for factor variable notation
        // The r(fvops) macro returns "true" or "false" as string
        fvexpand `covariates_list' if `touse'
        local has_fvops "`r(fvops)'"
        
        if "`has_fvops'" == "true" {
            // Each covariate term is processed: factors are expanded, base category excluded
            local expanded_vars ""
            local n_fv_expanded = 0
            
            foreach covar_term of local covariates_list {
                // Term is checked for factor variable expression syntax
                local is_factor = 0
                if regexm("`covar_term'", "^i\.") | regexm("`covar_term'", "^i\(") {
                    local is_factor = 1
                }
                // The ibn. prefix or ib# prefix is matched (including multi-digit base categories)
                if regexm("`covar_term'", "^ibn\.") | regexm("`covar_term'", "^ib[0-9]+") {
                    local is_factor = 1
                }
                if regexm("`covar_term'", "^c\.") {
                    // The c. prefix indicates a continuous variable
                    local is_factor = 0
                }
                
                if `is_factor' {
                    // Factor variable is expanded to indicator variables
                    fvrevar `covar_term' if `touse'
                    local fv_vars "`r(varlist)'"
                    local n_fv : word count `fv_vars'
                    
                    if `n_fv' > 1 {
                        // Exclude base category to avoid collinearity with the intercept term
                        forvalues i = 2/`n_fv' {
                            local v : word `i' of `fv_vars'
                            local expanded_vars "`expanded_vars' `v'"
                        }
                        local n_fv_expanded = `n_fv_expanded' + `n_fv' - 1
                    }
                    else if `n_fv' == 1 {
                        // Single-level factor is included (results in a constant)
                        local expanded_vars "`expanded_vars' `fv_vars'"
                        local n_fv_expanded = `n_fv_expanded' + 1
                    }
                }
                else {
                    // Numeric variables are retained without transformation
                    // The fvrevar command handles any remaining factor variable syntax
                    capture fvrevar `covar_term' if `touse'
                    if _rc == 0 {
                        local expanded_vars "`expanded_vars' `r(varlist)'"
                    }
                    else {
                        local expanded_vars "`expanded_vars' `covar_term'"
                    }
                }
            }
            
            local covariates_list "`expanded_vars'"
            local covariates_list : list retokenize covariates_list
            
            if `n_fv_expanded' > 0 {
                display as text "Note: Factor variables expanded to `n_fv_expanded' dummy variables (base categories excluded)"
            }
        }
        
        // Observations with missing covariate values are excluded
        markout `touse' `covariates_list'
    }
    
    // =========================================================================
    // SECTION 3: Design Routing
    // =========================================================================
    // Estimation is routed to staggered adoption (SA) or standard DID design
    
    if "`design_val'" == "sa" {
        // SA design is handled by the _diddesign_sa subprogram
        // Option string is constructed for SA estimation
        local sa_opts "treatment(`treatment_var') id(`id_var') time(`time_var')"
        local sa_opts "`sa_opts' nboot(`nboot_val') thres(`thres_val') level(`level_val')"
        
        if "`cluster_var'" != "" {
            local sa_opts "`sa_opts' cluster(`cluster_var')"
        }
        if "`covariates_list'" != "" {
            local sa_opts "`sa_opts' covariates(`covariates_list')"
        }
        if "`lead_val'" != "" {
            local sa_opts "`sa_opts' lead(`lead_val')"
        }
        if `seed_val' != . {
            local sa_opts "`sa_opts' seed(`seed_val')"
        }
        if "`quiet'" != "" {
            local sa_opts "`sa_opts' quiet"
        }
        
        // Original variable names are passed for e() reporting
        local sa_opts "`sa_opts' idorig(`id_var_orig') timeorig(`time_var_orig')"
        if "`cluster_var_orig'" != "" {
            local sa_opts "`sa_opts' clusterorig(`cluster_var_orig')"
        }
        
        // Global macro is used to pass cmdline (avoids parsing issues with special characters)
        global DIDDESIGN_CMDLINE `"`cmdline'"'
        local sa_opts "`sa_opts' touse(`touse')"
        
        // SA estimation subprogram is invoked
        _diddesign_sa `outcome' `covariates_list', `sa_opts'
        
        // Global macro is cleaned up after SA estimation completes
        capture macro drop DIDDESIGN_CMDLINE
        
        // Execution ends here; _diddesign_sa handles all e() returns and display
        exit
    }
    
    // Continue with standard DID design
    
    // =========================================================================
    // SECTION 4: Data Preparation
    // =========================================================================
    // Data structures are prepared for GMM estimation
    
    local prep_opts "outcome(`outcome') treatment(`treatment_var') time(`time_var')"
    
    if `is_panel' {
        local prep_opts "`prep_opts' id(`id_var') panel"
    }
    else {
        local prep_opts "`prep_opts' post(`post_var') rcs"
    }
    
    if "`cluster_var'" != "" {
        local prep_opts "`prep_opts' cluster(`cluster_var')"
    }
    if "`covariates_list'" != "" {
        local prep_opts "`prep_opts' covariates(`covariates_list')"
    }
    local prep_opts "`prep_opts' touse(`touse')"
    
    _diddesign_prep, `prep_opts'
    
    // Data preparation results are retrieved
    local N = r(N)
    local n_units = r(n_units)
    local n_periods = r(n_periods)
    local n_missing_delta = r(n_missing_delta)
    
    // =========================================================================
    // SECTION 5: GMM Estimation
    // =========================================================================
    // Double DID estimator via Generalized Method of Moments (GMM):
    //
    //   tau_ddid = argmin (m - tau)' W (m - tau)
    //
    // where m = (tau_DID, tau_sDID)' contains the standard DID and sequential
    // DID estimators, and W is the optimal GMM weight matrix (inverse of the
    // variance-covariance matrix of m). The Double DID achieves efficiency
    // under the parallel trends assumption and remains consistent under the
    // weaker parallel trends-in-trends assumption.
    
    // Random seed is set if specified
    if `seed_val' != . {
        set seed `seed_val'
    }
    
    // Lead numlist is converted to Mata format
    local lead_mata = subinstr("`lead_val'", " ", ", ", .)
    local n_lead : word count `lead_val'
    
    // Double DID estimation is executed via GMM
    mata: st_local("mata_rc", strofreal(_did_std_main((`lead_mata'), `nboot_val', `seboot_val', `level_val')))
    
    if `mata_rc' != 0 {
        // Specific error messages are provided based on error code
        if `mata_rc' == 1 {
            display as error "E011: Estimation failed - insufficient valid bootstrap iterations"
            display as error "      Try increasing the number of bootstrap iterations (nboot option)"
        }
        else if `mata_rc' == 2 {
            display as error "E011: Estimation failed - bootstrap VCOV computation failed"
            display as error "      This may be caused by insufficient valid bootstrap samples"
        }
        else {
            display as error "E011: Estimation failed in Mata (error code: `mata_rc')"
        }
        exit 498
    }
    
    // =========================================================================
    // SECTION 6: Store Estimation Results
    // =========================================================================
    // Results are stored in e() for post-estimation commands
    
    // --- Matrices ---
    // Matrices are retrieved from Mata first (before ereturn post clears them)
    tempname b_mat V_mat estimates_mat lead_mat weights_mat W_mat vcov_gmm_mat
    
    mata: st_matrix("`b_mat'", _did_b)
    mata: st_matrix("`V_mat'", _did_V)
    mata: st_matrix("`estimates_mat'", _did_estimates)
    mata: st_matrix("`lead_mat'", _did_lead_values)
    mata: st_matrix("`weights_mat'", _did_weights)
    
    // GMM weight matrix W and variance-covariance matrix of moment conditions
    mata: st_matrix("`W_mat'", _did_W)
    mata: st_matrix("`vcov_gmm_mat'", _did_vcov_gmm)
    
    // Estimation results are validated
    capture confirm matrix `lead_mat'
    if _rc != 0 {
        display as error "Error: Estimation produced no valid results (lead_mat not found)"
        exit 498
    }
    if colsof(`lead_mat') == 0 {
        display as error "Error: Estimation produced no valid results (lead_mat is empty)"
        exit 498
    }
    
    // For single lead, matrices are reshaped to 2x2
    local n_lead = colsof(`lead_mat')
    if `n_lead' == 1 {
        // W matrix is reshaped: [W11, W21, W12, W22] to 2x2 (vec() uses column-major order)
        matrix `W_mat' = (`W_mat'[1,1], `W_mat'[1,3] \ `W_mat'[1,2], `W_mat'[1,4])
        
        // VCOV matrix is reshaped: [V11, V21, V12, V22] to 2x2 (vec() uses column-major order)
        matrix `vcov_gmm_mat' = (`vcov_gmm_mat'[1,1], `vcov_gmm_mat'[1,3] \ `vcov_gmm_mat'[1,2], `vcov_gmm_mat'[1,4])
    }
    
    // Row and column names are set for e(b)
    local b_names ""
    foreach l of numlist `lead_val' {
        local b_names "`b_names' dDID:lead_`l' DID:lead_`l' sDID:lead_`l'"
    }
    // Leading space is trimmed
    local b_names = trim("`b_names'")
    matrix colnames `b_mat' = `b_names'
    
    // Row and column names are set for e(V)
    matrix rownames `V_mat' = `b_names'
    matrix colnames `V_mat' = `b_names'
    
    // Coefficient vector and variance-covariance matrix are posted
    ereturn post `b_mat' `V_mat', esample(`touse') obs(`N') depname("`outcome'")
    
    ereturn local properties "b V"
    
    // --- Scalars ---
    if `is_panel' {
        ereturn scalar n_units = `n_units'
    }
    ereturn scalar n_periods = `n_periods'
    ereturn scalar n_boot = `nboot_val'
    ereturn scalar level = `level_val'
    ereturn scalar n_lead = `n_lead'
    ereturn scalar is_panel = `is_panel'
    ereturn scalar seboot = `seboot_val'
    
    // n_boot_success only if some iterations failed
    mata: st_local("n_boot_success", strofreal(_did_n_boot_success))
    if "`n_boot_success'" != "" & "`n_boot_success'" != "." {
        if `n_boot_success' < `nboot_val' {
            ereturn scalar n_boot_success = `n_boot_success'
        }
    }
    
    // --- Macros ---
    ereturn local cmd "diddesign"
    ereturn local cmdline "`cmdline'"
    ereturn local design "`design_val'"
    ereturn local depvar "`outcome'"
    ereturn local treatment "`treatment_var'"
    ereturn local covariates "`covariates_list'"
    ereturn local id "`id_var_orig'"
    ereturn local time "`time_var_orig'"
    ereturn local clustvar "`cluster_var_orig'"
    ereturn local lead "`lead_val'"
    
    if `seboot_val' {
        ereturn local ci_method "bootstrap"
    }
    else {
        ereturn local ci_method "asymptotic"
    }
    
    // --- Additional Matrices (stored using ereturn matrix after ereturn post) ---
    // Row and column names are set for e(estimates)
    local est_rownames ""
    foreach l of numlist `lead_val' {
        local est_rownames "`est_rownames' dDID:lead_`l' DID:lead_`l' sDID:lead_`l'"
    }
    local est_rownames = trim("`est_rownames'")
    matrix rownames `estimates_mat' = `est_rownames'
    matrix colnames `estimates_mat' = lead estimate std_error ci_lo ci_hi weight
    
    // Names are set for e(lead_values)
    matrix colnames `lead_mat' = `lead_val'
    
    // Names are set for e(weights)
    local wt_rownames ""
    foreach l of numlist `lead_val' {
        local wt_rownames "`wt_rownames' lead_`l'"
    }
    local wt_rownames = trim("`wt_rownames'")
    matrix rownames `weights_mat' = `wt_rownames'
    matrix colnames `weights_mat' = w_did w_sdid
    
    // A copy of estimates_mat is made for display (before ereturn matrix moves it)
    tempname display_mat
    matrix `display_mat' = `estimates_mat'
    
    // Additional matrices are stored
    ereturn matrix estimates = `estimates_mat'
    ereturn matrix lead_values = `lead_mat'
    ereturn matrix weights = `weights_mat'
    ereturn matrix W = `W_mat'
    ereturn matrix vcov_gmm = `vcov_gmm_mat'
    
    // =========================================================================
    // SECTION 7: Display Results
    // =========================================================================
    
    // Header is displayed in Stata official two-column style
    if `is_panel' {
        local datatype "Panel"
        _diddesign_display_header, cmd("diddesign") design("std") ///
            datatype("`datatype'") n(`N') n_units(`n_units') ///
            n_periods(`n_periods') n_boot(`nboot_val') cluster("`cluster_var_orig'")
    }
    else {
        local datatype "Repeated Cross-Section"
        _diddesign_display_header, cmd("diddesign") design("std") ///
            datatype("`datatype'") n(`N') ///
            n_periods(`n_periods') n_boot(`nboot_val') cluster("`cluster_var_orig'")
    }
    
    // Abbreviate depvar for table header
    local depvar_abbrev = abbrev("`outcome'", 13)
    
    // Results table is displayed for each lead value
    local row = 1
    foreach l of numlist `lead_val' {
        // Table header with depvar~lead label
        if `n_lead' > 1 {
            local tbl_label "`depvar_abbrev'~`l'"
        }
        else {
            local tbl_label "`depvar_abbrev'"
        }
        
        _diddesign_display_table_header, depvar("`tbl_label'") level(`level_val') weight
        
        // Double DID
        local est = `display_mat'[`row', 2]
        local se = `display_mat'[`row', 3]
        local ci_lo = `display_mat'[`row', 4]
        local ci_hi = `display_mat'[`row', 5]
        local wt = `display_mat'[`row', 6]
        _diddesign_display_result, label("Double DID") ///
            estimate(`est') se(`se') ci_low(`ci_lo') ci_high(`ci_hi') weight(`wt')
        local row = `row' + 1
        
        // DID
        local est = `display_mat'[`row', 2]
        local se = `display_mat'[`row', 3]
        local ci_lo = `display_mat'[`row', 4]
        local ci_hi = `display_mat'[`row', 5]
        local wt = `display_mat'[`row', 6]
        _diddesign_display_result, label("DID") ///
            estimate(`est') se(`se') ci_low(`ci_lo') ci_high(`ci_hi') weight(`wt')
        local row = `row' + 1
        
        // sDID
        local est = `display_mat'[`row', 2]
        local se = `display_mat'[`row', 3]
        local ci_lo = `display_mat'[`row', 4]
        local ci_hi = `display_mat'[`row', 5]
        local wt = `display_mat'[`row', 6]
        _diddesign_display_result, label("sDID") ///
            estimate(`est') se(`se') ci_low(`ci_lo') ci_high(`ci_hi') weight(`wt')
        local row = `row' + 1
    }
    
    display as text "{hline 78}"
    
    // Footer note
    display as text ""
    display as text "Note: Double DID combines DID and sDID via optimal GMM weights."
    
    // =========================================================================
    // SECTION 8: Cleanup
    // =========================================================================
    
    capture drop _did_id_time
    capture drop _did_id_time_std
    capture drop _did_Gi
    capture drop _did_It
    capture drop _did_outcome_delta
    capture drop _did_Ymean
    
end

