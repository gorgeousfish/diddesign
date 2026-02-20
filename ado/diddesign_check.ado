*! diddesign_check.ado - Diagnostic tests for parallel trends assumption
*!
*! Implements placebo tests for assessing the parallel trends assumption in
*! difference-in-differences designs. Computes standardized pre-treatment DID
*! estimates and equivalence confidence intervals for both standard DID and
*! staggered adoption designs.

program define diddesign_check, eclass
    version 16.0
    
    // -------------------------------------------------------------------------
    // Load Mata Functions
    // -------------------------------------------------------------------------
    capture mata: mata describe did_placebo()
    if _rc != 0 {
        local mata_loaded = 0
        
        // Method 1: Direct findfile for diddesign_mata.do (works after net install)
        qui capture findfile diddesign_mata.do
        if _rc == 0 {
            quietly do "`r(fn)'"
            local mata_loaded = 1
        }
        
        // Method 2: Relative path from ado file (works in development environment)
        if !`mata_loaded' {
            qui capture findfile diddesign_check.ado
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
            capture mata: mata describe did_placebo()
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
    
    // -------------------------------------------------------------------------
    // Syntax Parsing
    // -------------------------------------------------------------------------
    // Supports both panel data (with id()) and repeated cross-section (with post())
    // Standardization is always performed, returning both standardized and raw estimates
    syntax varlist(min=1 fv) [if] [in], ///
        TREATment(varname)              /// Required: treatment indicator
        TIME(varname)                   /// Required: time identifier
        [ID(varname)]                   /// Unit identifier (required for panel)
        [POST(varname)]                 /// Post-treatment indicator (required for RCS)
        [DESIGN(string)]                /// Design type: "did" (default) or "sa"
        [PANEL]                         /// Panel data format
        [RCS]                           /// Repeated cross-section format
        [CLuster(varname)]              /// Cluster variable for SEs
        [NBoot(integer 30)]             /// Bootstrap iterations; default is 30
        [LAG(numlist >0 integer)]       /// Lag values for placebo tests
        [THRes(integer 2)]              /// SA threshold; default is 2
        [PARALlel]                      /// Use parallel computing
        [SEED(integer -1)]              /// Random seed (-1 = not specified)
        [QUIET]                         /// Suppress progress display
    
    // Store full command line for e(cmdline)
    local cmdline "diddesign_check `0'"
    
    // -------------------------------------------------------------------------
    // Parse Variable List
    // -------------------------------------------------------------------------
    // The varlist contains the outcome variable followed by optional covariates
    gettoken depvar covariates : varlist
    
    // Handle covariates (may be empty)
    if "`covariates'" == "" {
        local covars_str ""
    }
    else {
        local covars_str "`covariates'"
    }
    
    // -------------------------------------------------------------------------
    // Duplicate Covariate Check
    // -------------------------------------------------------------------------
    // Remove duplicate covariates to ensure proper model specification
    if "`covars_str'" != "" {
        local unique_covars : list uniq covars_str
        local n_all : word count `covars_str'
        local n_unique : word count `unique_covars'
        
        if `n_unique' < `n_all' {
            // Find duplicate variables by comparing original and unique lists
            local dups ""
            local seen ""
            foreach v of local covars_str {
                local is_seen : list v in seen
                if `is_seen' {
                    local is_dup : list v in dups
                    if !`is_dup' {
                        local dups "`dups' `v'"
                    }
                }
                local seen "`seen' `v'"
            }
            local dups = strtrim("`dups'")
            display as text "Warning: Duplicate covariates detected and removed: `dups'"
            local covars_str "`unique_covars'"
        }
    }
    
    // -------------------------------------------------------------------------
    // Early Validation: Treatment Variable
    // -------------------------------------------------------------------------
    // Verify treatment has both 0 and 1 values before expensive operations
    quietly count if `treatment' == 1
    if r(N) == 0 {
        display as error "E003: No treated observations found in data (treatment all 0)"
        display as error "       Placebo tests require at least one treated unit"
        exit 198
    }
    quietly count if `treatment' == 0
    if r(N) == 0 {
        display as error "E003: No control observations found in data (treatment all 1)"
        display as error "       Placebo tests require at least one control unit"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Factor Variable Expansion
    // -------------------------------------------------------------------------
    // Expand factor variables (i.var, ibn.var) into dummy variables
    // Base category is excluded to avoid collinearity with intercept
    if "`covars_str'" != "" {
        // Create temporary sample marker for fvrevar
        marksample touse_temp, novarlist
        markout `touse_temp' `treatment' `time'
        if "`id'" != "" {
            markout `touse_temp' `id'
        }
        
        // Check if covariates contain factor variable notation
        // Note: r(fvops) returns "true" or "false" as string, not numeric
        fvexpand `covars_str' if `touse_temp'
        local has_fvops "`r(fvops)'"
        
        if "`has_fvops'" == "true" {
            // Process each covariate term separately
            // For factor variables: expand and exclude base category (first dummy)
            local expanded_vars ""
            local n_fv_expanded = 0
            
            foreach covar_term of local covars_str {
                // Check if this term is a factor variable expression
                local is_factor = 0
                if regexm("`covar_term'", "^i\.") | regexm("`covar_term'", "^i\(") {
                    local is_factor = 1
                }
                // Check for ibn. or ib# prefix (supports multi-digit base categories)
                if regexm("`covar_term'", "^ibn\.") | regexm("`covar_term'", "^ib[0-9]+") {
                    local is_factor = 1
                }
                if regexm("`covar_term'", "^c\.") {
                    local is_factor = 0
                }
                
                if `is_factor' {
                    fvrevar `covar_term' if `touse_temp'
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
                    capture fvrevar `covar_term' if `touse_temp'
                    if _rc == 0 {
                        local expanded_vars "`expanded_vars' `r(varlist)'"
                    }
                    else {
                        local expanded_vars "`expanded_vars' `covar_term'"
                    }
                }
            }
            
            local covars_str "`expanded_vars'"
            local covars_str : list retokenize covars_str
            
            if `n_fv_expanded' > 0 {
                display as text "Note: Factor variables expanded to `n_fv_expanded' dummy variables (base categories excluded)"
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Set Default Values
    // -------------------------------------------------------------------------
    // Default design: standard DID
    if "`design'" == "" {
        local design "did"
    }
    else {
        local design = lower("`design'")
    }
    
    // Default lag: 1 period
    if "`lag'" == "" {
        local lag "1"
    }
    
    // Bootstrap iterations: defaults to 30
    local nboot_val = `nboot'
    
    // Staggered adoption threshold: defaults to 2
    local thres_val = `thres'
    
    // Random seed: -1 indicates user did not specify
    if `seed' == -1 {
        local seed_val = .
    }
    else {
        local seed_val = `seed'
    }
    
    // Quiet option - suppress bootstrap progress display
    local quiet_val = 0
    if "`quiet'" != "" {
        local quiet_val = 1
    }
    
    // Parallel option is not yet implemented; bootstrap runs sequentially
    if "`parallel'" != "" {
        display as text "{p 0 4 2}"
        display as text "Note: The {bf:parallel} option is currently not implemented. "
        display as text "Bootstrap iterations will run sequentially.{p_end}"
    }
    
    // -------------------------------------------------------------------------
    // Validate Data Type
    // -------------------------------------------------------------------------
    // Determine data type: panel vs repeated cross-section (RCS)
    // Panel and rcs options are mutually exclusive
    local is_panel_opt = ("`panel'" != "")
    local is_rcs_opt = ("`rcs'" != "")
    
    if `is_panel_opt' & `is_rcs_opt' {
        display as error "E016: Options panel and rcs are mutually exclusive"
        exit 198
    }
    
    // Determine data type (auto-detect if neither specified)
    local is_panel = 0
    if `is_panel_opt' {
        local is_panel = 1
    }
    else if `is_rcs_opt' {
        local is_panel = 0
    }
    else {
        // Auto-detect: panel if id() is specified, RCS if post() is specified
        if "`id'" != "" {
            local is_panel = 1
        }
        else if "`post'" != "" {
            local is_panel = 0
        }
        else {
            // Neither id() nor post() specified - require explicit choice
            display as error "E016: Must specify id() for panel data or post() with rcs option for RCS data"
            exit 198
        }
    }
    
    // Validate required options based on data type
    if `is_panel' {
        // Panel data requires id()
        if "`id'" == "" {
            display as error "E001: Option id() is required for panel data"
            exit 198
        }
    }
    else {
        // RCS data requires post() to identify treatment timing
        if "`post'" == "" {
            display as error "E001: Option post() is required for RCS data"
            display as error "       Specify the post-treatment indicator variable"
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Validate Parameters
    // -------------------------------------------------------------------------
    // nboot >= 2 required for variance estimation (denominator is n_boot - 1)
    if `nboot_val' < 2 {
        display as error "E002: Option nboot() must be at least 2 for variance estimation"
        exit 198
    }
    
    // Validate thres >= 1 (only used in SA design)
    if `thres_val' < 1 {
        display as error "E002: Option thres() must be a positive integer >= 1"
        exit 198
    }
    
    // Validate seed range if specified (must be in [0, 2147483647])
    if `seed' != -1 & (`seed' < 0 | `seed' > 2147483647) {
        display as error "E002: Option seed() must be a valid integer (0 to 2147483647)"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Validate Design Type
    // -------------------------------------------------------------------------
    if !inlist("`design'", "did", "sa") {
        display as error "E002: design() must be 'did' or 'sa'"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Handle Cluster Variable
    // -------------------------------------------------------------------------
    // For panel data, default cluster to unit identifier if not specified
    local clustvar "`cluster'"
    if "`clustvar'" == "" & `is_panel' & "`id'" != "" {
        local clustvar "`id'"
    }
    
    // String cluster variables must be encoded to numeric before markout
    if "`clustvar'" != "" {
        capture confirm string variable `clustvar'
        if _rc == 0 {
            tempvar cluster_encoded
            quietly egen `cluster_encoded' = group(`clustvar')
            display as text "Note: String variable `clustvar' automatically encoded to numeric"
            local clustvar "`cluster_encoded'"
        }
    }
    
    // -------------------------------------------------------------------------
    // Validate SA Design Requires Panel Data
    // -------------------------------------------------------------------------
    // Staggered adoption design only supports panel data structure
    if "`design'" == "sa" & !`is_panel' {
        display as error "E014: SA design requires panel data"
        display as error "       Only the standard DID design supports RCS data"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Mark Sample
    // -------------------------------------------------------------------------
    marksample touse, novarlist
    markout `touse' `depvar' `treatment' `time'
    // Panel data: mark out missing id
    if `is_panel' & "`id'" != "" {
        markout `touse' `id'
    }
    // RCS data: mark out missing post indicator
    if !`is_panel' & "`post'" != "" {
        markout `touse' `post'
    }
    if "`clustvar'" != "" {
        markout `touse' `clustvar'
    }
    if "`covars_str'" != "" {
        markout `touse' `covars_str'
    }
    
    // Count observations
    quietly count if `touse'
    local N = r(N)
    
    if `N' == 0 {
        display as error "E003: No observations"
        exit 2000
    }
    
    // -------------------------------------------------------------------------
    // Compute Number of Clusters
    // -------------------------------------------------------------------------
    if "`clustvar'" != "" {
        tempvar cluster_tag
        quietly egen `cluster_tag' = tag(`clustvar') if `touse'
        quietly count if `cluster_tag' == 1 & `touse'
        local n_clusters = r(N)
    }
    else {
        local n_clusters = `N'
    }
    
    // -------------------------------------------------------------------------
    // Set Random Seed
    // -------------------------------------------------------------------------
    if `seed_val' != . {
        set seed `seed_val'
    }
    
    // -------------------------------------------------------------------------
    // Validate Variables
    // -------------------------------------------------------------------------
    // Validate outcome variable
    capture confirm numeric variable `depvar'
    if _rc {
        display as error "E017: Variable `depvar' must be numeric"
        exit _rc
    }
    
    // Validate treatment variable
    capture confirm numeric variable `treatment'
    if _rc {
        display as error "E017: Variable `treatment' must be numeric"
        exit _rc
    }
    
    // Treatment must be strictly binary (0/1) for valid group matrix construction
    quietly tabulate `treatment' if `touse'
    if r(r) > 2 {
        display as error "E003: Treatment variable must be binary (0/1)"
        display as error "      Found `r(r)' distinct values"
        exit 198
    }
    quietly summarize `treatment' if `touse'
    // Check for all-missing values before numeric comparison
    if r(N) == 0 | missing(r(min)) {
        display as error "E003: Treatment variable contains only missing values"
        exit 198
    }
    // Use floating-point tolerance for binary validation (handles CSV import precision)
    if round(r(min), 1e-6) != 0 | round(r(max), 1e-6) != 1 {
        display as error "E003: Treatment variable must contain both 0 and 1 values"
        display as error "      Found min=`r(min)', max=`r(max)'"
        exit 198
    }
    
    // Validate id variable (only for panel data)
    // String variables are automatically encoded to numeric
    if `is_panel' & "`id'" != "" {
        capture confirm variable `id'
        if _rc {
            display as error "E001: Variable `id' not found"
            exit 111
        }
        // Check if string variable - auto convert using egen group()
        capture confirm string variable `id'
        if _rc == 0 {
            // String variable detected - create numeric encoding
            tempvar id_encoded
            quietly egen `id_encoded' = group(`id')
            display as text "Note: String variable `id' automatically encoded to numeric"
            local id "`id_encoded'"
        }
        else {
            capture confirm numeric variable `id'
            if _rc {
                display as error "E017: Variable `id' must be numeric or string"
                exit _rc
            }
        }
    }
    
    // Validate time variable
    // String variables are automatically encoded to numeric
    capture confirm variable `time'
    if _rc {
        display as error "E001: Variable `time' not found"
        exit 111
    }
    // Check if string variable - auto convert using egen group()
    capture confirm string variable `time'
    if _rc == 0 {
        // String variable detected - create numeric encoding
        tempvar time_encoded
        quietly egen `time_encoded' = group(`time')
        display as text "Note: String variable `time' automatically encoded to numeric"
        local time "`time_encoded'"
    }
    else {
        capture confirm numeric variable `time'
        if _rc {
            display as error "E017: Variable `time' must be numeric or string"
            exit _rc
        }
    }
    
    // Validate post variable (only for RCS data)
    if !`is_panel' & "`post'" != "" {
        capture confirm numeric variable `post'
        if _rc {
            display as error "E017: Variable `post' must be numeric"
            exit _rc
        }
        
        // Post indicator must be strictly binary (0/1)
        // Check only contains 0 and 1 values (not continuous)
        quietly tab `post' if `touse'
        if r(r) > 2 {
            display as error "E003: Post-treatment indicator must be binary (0/1)"
            display as error "       Found " r(r) " distinct values (expected 2)"
            exit 198
        }
        
        // Verify values are exactly 0 and 1 (with floating-point tolerance)
        quietly summarize `post' if `touse'
        if (round(r(min), 1e-6) != 0 & round(r(min), 1e-6) != 1) | (round(r(max), 1e-6) != 0 & round(r(max), 1e-6) != 1) {
            display as error "E003: Post-treatment indicator must contain exactly 0 and 1"
            display as error "       Found min=" r(min) ", max=" r(max)
            exit 198
        }
        
        // Ensure both 0 and 1 values exist for valid DID estimation
        if r(min) == r(max) {
            if r(min) == 0 {
                display as error "E003: Post-treatment indicator is all 0 (no post-treatment observations)"
                display as error "       Placebo tests require at least one post-treatment period"
            }
            else {
                display as error "E003: Post-treatment indicator is all 1 (no pre-treatment observations)"
                display as error "       Placebo tests require at least one pre-treatment period"
            }
            exit 198
        }
    }
    
    // Validate covariates
    if "`covars_str'" != "" {
        foreach var of local covars_str {
            capture confirm numeric variable `var'
            if _rc {
                display as error "E017: Variable `var' must be numeric"
                exit _rc
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Call Mata for Computation
    // -------------------------------------------------------------------------
    // Convert lag numlist to Mata format
    local lag_mata = subinstr("`lag'", " ", ", ", .)
    local n_lags : word count `lag'
    
    // Set id_var for Mata (empty string for RCS)
    local id_var = ""
    if `is_panel' & "`id'" != "" {
        local id_var "`id'"
    }
    
    // Set post_var for Mata (empty string for panel)
    local post_var = ""
    if !`is_panel' & "`post'" != "" {
        local post_var "`post'"
    }
    
    // Call Mata main function for placebo test computation
    mata: _diddesign_check_main( ///
        "`depvar'",           ///
        "`treatment'",        ///
        "`id_var'",           ///
        "`time'",             ///
        "`post_var'",         ///
        "`covars_str'",       ///
        "`clustvar'",         ///
        "`touse'",            ///
        "`design'",           ///
        (`lag_mata'),         ///
        `nboot_val',          ///
        `thres_val',          ///
        `is_panel',           ///
        `quiet_val'           ///
    )
    
    // -------------------------------------------------------------------------
    // Store e() Returns
    // -------------------------------------------------------------------------
    // Retrieve matrices from Mata
    tempname placebo_mat trends_mat Gmat_mat
    
    mata: st_matrix("`placebo_mat'", _check_placebo)
    mata: st_matrix("`trends_mat'", _check_trends)
    
    // Get scalar results
    mata: st_local("n_lags_valid", strofreal(_check_n_lags))
    mata: st_local("n_boot_valid", strofreal(_check_n_boot_valid))
    mata: st_local("filtered_lags", _check_filtered_lags)
    
    // Set matrix row and column names for e(placebo)
    local placebo_rownames ""
    forvalues i = 1/`n_lags_valid' {
        local lag_i = `placebo_mat'[`i', 1]
        local lag_i_int = int(`lag_i')
        local placebo_rownames "`placebo_rownames' `lag_i_int'"
    }
    local placebo_rownames = trim("`placebo_rownames'")
    if "`placebo_rownames'" != "" {
        matrix rownames `placebo_mat' = `placebo_rownames'
    }
    matrix colnames `placebo_mat' = lag estimate std_error estimate_orig std_error_orig EqCI95_LB EqCI95_UB
    
    // Set matrix column names for e(trends)
    matrix colnames `trends_mat' = id_time_std Gi outcome_mean outcome_sd n_obs
    
    // Post results (no b/V matrices for check command)
    ereturn clear
    
    // --- Scalars ---
    ereturn scalar N = `N'
    ereturn scalar n_lags = `n_lags_valid'
    ereturn scalar n_boot = `nboot_val'
    ereturn scalar n_boot_valid = `n_boot_valid'
    ereturn scalar n_clusters = `n_clusters'
    
    // --- Macros ---
    ereturn local cmd "diddesign_check"
    ereturn local cmdline "`cmdline'"
    ereturn local design "`design'"
    ereturn local depvar "`depvar'"
    ereturn local treatment "`treatment'"
    ereturn local clustvar "`clustvar'"
    ereturn local covars "`covars_str'"
    // This command stores e(placebo) and e(trends), not e(b) and e(V)
    
    // --- Matrices ---
    ereturn matrix placebo = `placebo_mat'
    ereturn matrix trends = `trends_mat'
    
    // SA design: store treatment timing matrix (Gmat) only if valid
    // Invalid or placeholder matrices are not stored to prevent misleading plots
    if "`design'" == "sa" {
        capture mata: st_matrix("`Gmat_mat'", _check_Gmat)
        if _rc == 0 {
            capture confirm matrix `Gmat_mat'
            if _rc == 0 {
                local gmat_rows = rowsof(`Gmat_mat')
                local gmat_cols = colsof(`Gmat_mat')
                // Only store Gmat if it represents valid data (not 1x1 placeholder)
                if (`gmat_rows' > 1 | `gmat_cols' > 1) & `gmat_rows' > 0 & `gmat_cols' > 0 {
                    ereturn matrix Gmat = `Gmat_mat'
                }
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Store Data Type Info
    // -------------------------------------------------------------------------
    ereturn scalar is_panel = `is_panel'
    if `is_panel' {
        ereturn local datatype "panel"
    }
    else {
        ereturn local datatype "rcs"
    }
    
    // -------------------------------------------------------------------------
    // Display Results
    // -------------------------------------------------------------------------
    _diddesign_check_display, design("`design'") filtered_lags("`filtered_lags'") ///
        is_panel(`is_panel')
    
end

// =============================================================================
// _diddesign_check_display
// Formats and displays placebo test results
//
// Displays the parallel trends assessment output including design information,
// sample statistics, and a formatted table of placebo estimates with bootstrap
// standard errors and equivalence confidence intervals.
// =============================================================================
program define _diddesign_check_display
    syntax, design(string) [filtered_lags(string) is_panel(integer 1)]
    
    // Header in Stata official two-column style
    if `is_panel' {
        local datatype "Panel"
    }
    else {
        local datatype "Repeated Cross-Section"
    }
    
    _diddesign_display_header, cmd("diddesign_check") design("`design'") ///
        datatype("`datatype'") n(`=e(N)') n_boot(`=e(n_boot)')
    
    // Warning for filtered lags
    if "`filtered_lags'" != "" {
        display ""
        display as text "{p 0 4 2}"
        display as error "Warning: " as text "The following lag(s) were filtered out (exceed max available): `filtered_lags'"
        display as text "{p_end}"
    }
    
    // Table header
    display ""
    display as text "{hline 78}"
    display as text "    Lag|  Estimate  Std. err.   Est.(raw)   SE(raw)       [95% eq. CI]"
    display as text "{hline 8}+{hline 69}"
    
    // Table content
    tempname placebo
    matrix `placebo' = e(placebo)
    local nrows = rowsof(`placebo')
    
    if `nrows' == 0 {
        display as text "    (no valid lags)"
    }
    else {
        forvalues i = 1/`nrows' {
            local lag_val = `placebo'[`i', 1]
            local est = `placebo'[`i', 2]
            local se = `placebo'[`i', 3]
            local est_orig = `placebo'[`i', 4]
            local se_orig = `placebo'[`i', 5]
            local ci_lb = `placebo'[`i', 6]
            local ci_ub = `placebo'[`i', 7]
            
            display as text %7.0f `lag_val' "|" ///
                as result %10.4f `est' %10.4f `se' ///
                %12.4f `est_orig' %10.4f `se_orig' ///
                "   [" %7.4f `ci_lb' ", " %7.4f `ci_ub' "]"
        }
    }
    
    display as text "{hline 78}"
    
end
