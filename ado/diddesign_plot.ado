*! diddesign_plot.ado - Visualization command for DIDdesign estimation results
*!
*! Generates diagnostic plots for Double DID analysis including trend plots
*! for standard DID designs and treatment pattern visualizations for
*! staggered adoption designs.


// =============================================================================
// diddesign_plot
// Main entry point for DIDdesign visualization
//
// Creates five types of plots for Double DID analysis:
//   - estimates : Double-DID estimates across lead values with 90% CI
//   - trends    : Outcome trajectories for treated and control groups
//   - placebo   : Pre-treatment placebo test with equivalence CI
//   - pattern   : Treatment timing heatmap (staggered adoption only)
//   - both      : Combined diagnostic plot (placebo + trends/pattern)
//
// Options:
//   type()     : string - Plot type (default: estimates after diddesign,
//                         both after diddesign_check)
//   saving()   : string - Output file path
//   replace    : flag   - Overwrite existing file
//   scheme()   : string - Graph scheme name
//   title()    : string - Graph title
//   xtitle()   : string - X-axis title
//   ytitle()   : string - Y-axis title
//   ci         : flag   - Display confidence bands on trends plot
//   band       : flag   - Display CI as ribbon on estimates plot
//   name()     : string - Graph name in memory
//   use_check(): string - Name of stored diddesign_check results to overlay
//
// Data source:
//   Reads from e() results produced by diddesign or diddesign_check
// =============================================================================

program define diddesign_plot
    version 16.0
    
    syntax [, TYPE(string) SAVing(string) REPlace SCHeme(string) ///
              TItle(string) XTItle(string) YTItle(string) ///
              CI BAND name(string) USE_CHECK(string)]
    
    // -------------------------------------------------------------------------
    // Validate estimation results
    // -------------------------------------------------------------------------
    
    if "`e(cmd)'" == "" {
        display as error "no estimation results found"
        display as error "run diddesign or diddesign_check first"
        exit 301
    }
    
    if !inlist("`e(cmd)'", "diddesign", "diddesign_check") {
        display as error "diddesign_plot requires diddesign or diddesign_check results"
        display as error "run diddesign or diddesign_check first"
        exit 301
    }
    
    // Source command type is stored for subsequent routing
    local current_cmd "`e(cmd)'"
    
    // -------------------------------------------------------------------------
    // Get design type
    // -------------------------------------------------------------------------
    local design = "`e(design)'"
    if "`design'" == "" {
        local design "did"
    }
    
    // -------------------------------------------------------------------------
    // Set default plot type based on source command
    // -------------------------------------------------------------------------
    // Default type is "estimates" for diddesign results and "both" for
    // diddesign_check results
    
    if "`type'" == "" {
        if "`current_cmd'" == "diddesign" {
            local type "estimates"
        }
        else {
            local type "both"
        }
    }
    
    // Validate type option
    if !inlist("`type'", "trends", "placebo", "pattern", "both", "estimates") {
        display as error "type() must be one of: trends, placebo, pattern, both, estimates"
        exit 198
    }
    
    // -------------------------------------------------------------------------
    // Validate type-design compatibility
    // -------------------------------------------------------------------------
    // Pattern plot is only available for staggered adoption design
    // Trends plot is only available for standard DID design
    
    if "`type'" == "pattern" & "`design'" != "sa" {
        display as error "pattern plot is only available for SA design"
        exit 198
    }
    
    if "`type'" == "trends" & "`design'" == "sa" {
        display as error "trends plot is only available for standard DID design"
        exit 198
    }
    
    // Diagnostic plots require diddesign_check results
    if inlist("`type'", "trends", "placebo", "pattern", "both") & "`current_cmd'" != "diddesign_check" {
        display as error "type(`type') requires diddesign_check results"
        display as error "run diddesign_check first, or use type(estimates)"
        exit 198
    }
    
    // Estimates plot requires e(estimates) matrix
    if "`type'" == "estimates" {
        capture confirm matrix e(estimates)
        if _rc {
            display as error "e(estimates) not found; type(estimates) requires diddesign results"
            exit 301
        }
    }
    
    // -------------------------------------------------------------------------
    // Build common options
    // -------------------------------------------------------------------------
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt `"scheme("`scheme'")"'
    }
    
    // -------------------------------------------------------------------------
    // Dispatch to appropriate plot function
    // -------------------------------------------------------------------------
    
    if "`type'" == "trends" {
        _plot_trends, `scheme_opt' title("`title'") xtitle("`xtitle'") ///
                      ytitle("`ytitle'") `ci' name("`name'")
        if "`saving'" != "" {
            _export_graph, saving("`saving'") `replace'
        }
    }
    else if "`type'" == "placebo" {
        _plot_placebo, `scheme_opt' title("`title'") xtitle("`xtitle'") ///
                       ytitle("`ytitle'") name("`name'")
        if "`saving'" != "" {
            _export_graph, saving("`saving'") `replace'
        }
    }
    else if "`type'" == "pattern" {
        _plot_pattern, `scheme_opt' title("`title'") xtitle("`xtitle'") ///
                       ytitle("`ytitle'") name("`name'")
        if "`saving'" != "" {
            _export_graph, saving("`saving'") `replace'
        }
    }
    else if "`type'" == "both" {
        // Combined diagnostic plot is generated
        _plot_combined, design("`design'") `scheme_opt' ///
                        title("`title'") saving("`saving'") `ci' `replace'
    }
    else if "`type'" == "estimates" {
        // Double-DID estimates plot across lead values
        local use_check_opt ""
        if "`use_check'" != "" {
            local use_check_opt `"use_check("`use_check'")"'
        }
        _plot_estimates, design("`design'") `scheme_opt' ///
                         title("`title'") xtitle("`xtitle'") ytitle("`ytitle'") ///
                         name("`name'") `use_check_opt' `band'
        if "`saving'" != "" {
            _export_graph, saving("`saving'") `replace'
        }
    }
    
end


// =============================================================================
// _plot_combined()
// Generate combined diagnostic plot with two subplots arranged horizontally
//
// Layout:
//   - Standard DID: placebo (left) + trends (right)
//   - Staggered adoption: placebo (left) + pattern (right)
//
// Arguments:
//   design  : string - Design type ("did" or "sa")
//   scheme  : string - Graph scheme name
//   title   : string - Combined graph title
//   saving  : string - Output file path
//   name    : string - Graph name in memory
//   ci      : flag   - Display confidence intervals in trends subplot
//   replace : flag   - Overwrite existing file
// =============================================================================

program define _plot_combined
    version 16.0
    
    syntax , DESIGN(string) [SCHeme(string) TItle(string) SAVing(string) NAME(string) CI REPlace]
    
    // -------------------------------------------------------------------------
    // Build scheme option
    // -------------------------------------------------------------------------
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt `"scheme(`scheme')"'
    }
    
    // -------------------------------------------------------------------------
    // Generate subplots based on design type
    // -------------------------------------------------------------------------
    // Placebo plot is positioned on the left for both designs
    
    if "`design'" == "sa" {
        // Staggered adoption: placebo (left) + pattern (right)
        capture graph drop _tmp_placebo
        capture graph drop _tmp_pattern
        
        quietly _plot_placebo, `scheme_opt' name(_tmp_placebo)
        quietly _plot_pattern, `scheme_opt' name(_tmp_pattern)
        
        local title_opt ""
        if `"`title'"' != "" {
            local title_opt `"title(`"`title'"')"'
        }
        
        local graph_name "combined_plot"
        if "`name'" != "" {
            local graph_name "`name'"
        }
        
        graph combine _tmp_placebo _tmp_pattern, ///
            rows(1) ///
            `title_opt' ///
            `scheme_opt' ///
            name(`graph_name', replace)
        
        capture graph drop _tmp_placebo
        capture graph drop _tmp_pattern
    }
    else {
        // Standard DID: placebo (left) + trends (right)
        capture graph drop _tmp_placebo
        capture graph drop _tmp_trends
        
        quietly _plot_placebo, `scheme_opt' name(_tmp_placebo)
        quietly _plot_trends, `scheme_opt' name(_tmp_trends) `ci'
        
        local title_opt ""
        if `"`title'"' != "" {
            local title_opt `"title(`"`title'"')"'
        }
        
        local graph_name "combined_plot"
        if "`name'" != "" {
            local graph_name "`name'"
        }
        
        graph combine _tmp_placebo _tmp_trends, ///
            rows(1) ///
            `title_opt' ///
            `scheme_opt' ///
            name(`graph_name', replace)
        
        capture graph drop _tmp_placebo
        capture graph drop _tmp_trends
    }
    
    // -------------------------------------------------------------------------
    // Save graph if requested
    // -------------------------------------------------------------------------
    if "`saving'" != "" {
        _export_graph, saving("`saving'") `replace'
    }
    
end


// =============================================================================
// _export_graph()
// Export graph to file with automatic format detection
//
// Arguments:
//   saving  : string - Output file path
//   replace : flag   - Overwrite existing file
//
// Supported formats: .png, .pdf, .eps, .svg, .tif
// Default format is PNG if no recognized extension is provided.
// =============================================================================

program define _export_graph
    version 16.0
    
    syntax , SAVing(string) [REPlace]
    
    // Determine file extension and export format
    local ext = substr("`saving'", -4, .)
    if inlist("`ext'", ".png", ".pdf", ".eps", ".svg", ".tif") {
        graph export "`saving'", `replace'
    }
    else {
        // Default to PNG if no recognized extension
        graph export "`saving'.png", `replace'
    }
    
end



// =============================================================================
// _plot_trends()
// Generate trend plot showing outcome means for treated and control groups
//
// Displays outcome trajectories over time relative to treatment assignment.
// Control and treated groups are shown as connected lines with markers.
// A vertical dashed line at x=0 indicates the treatment timing.
//
// Arguments:
//   saving  : string - Output file path
//   scheme  : string - Graph scheme name
//   title   : string - Graph title
//   xtitle  : string - X-axis title (default: "Time relative to treatment assignment")
//   ytitle  : string - Y-axis title (default: "Mean Outcome")
//   ci      : flag   - Display 90% confidence bands
//   name    : string - Graph name in memory (default: "trends_plot")
//
// Data source:
//   e(trends) matrix with columns: id_time_std, Gi, outcome_mean, outcome_sd, n_obs
// =============================================================================

program define _plot_trends
    version 16.0
    
    syntax [, SAVing(string) SCHeme(string) ///
              TItle(string) XTItle(string) YTItle(string) ///
              CI name(string)]
    
    // -------------------------------------------------------------------------
    // Validate e(trends) exists
    // -------------------------------------------------------------------------
    if "`e(cmd)'" == "" {
        display as error "no estimation results found"
        exit 301
    }
    
    capture confirm matrix e(trends)
    if _rc {
        display as error "e(trends) not found; run diddesign or diddesign_check first"
        exit 301
    }
    
    // -------------------------------------------------------------------------
    // Extract data to temporary dataset
    // -------------------------------------------------------------------------
    preserve
    clear
    
    tempname trends
    matrix `trends' = e(trends)
    local nrows = rowsof(`trends')
    
    if `nrows' == 0 {
        display as error "e(trends) matrix is empty"
        restore
        exit 198
    }
    
    quietly {
        set obs `nrows'
        
        gen double id_time_std = .
        gen byte Gi = .
        gen double outcome_mean = .
        gen double outcome_sd = .
        gen long n_obs = .
        
        forvalues i = 1/`nrows' {
            replace id_time_std = `trends'[`i', 1] in `i'
            replace Gi = `trends'[`i', 2] in `i'
            replace outcome_mean = `trends'[`i', 3] in `i'
            replace outcome_sd = `trends'[`i', 4] in `i'
            replace n_obs = `trends'[`i', 5] in `i'
        }
        
        gen str8 group = cond(Gi == 1, "Treated", "Control")
    }
    
    // -------------------------------------------------------------------------
    // Calculate confidence intervals (optional)
    // -------------------------------------------------------------------------
    // 90% CI: outcome_mean +/- invnormal(0.95) * outcome_sd
    
    if "`ci'" != "" {
        quietly {
            gen double outcome_ci_lb = outcome_mean - invnormal(0.95) * outcome_sd
            gen double outcome_ci_ub = outcome_mean + invnormal(0.95) * outcome_sd
        }
    }
    
    // -------------------------------------------------------------------------
    // Set default values
    // -------------------------------------------------------------------------
    if "`title'" == "" local title ""
    if "`xtitle'" == "" local xtitle "Time relative to treatment assignment"
    if "`ytitle'" == "" local ytitle "Mean Outcome"
    if "`name'" == "" local name "trends_plot"
    
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt "scheme(`scheme')"
    }
    
    // -------------------------------------------------------------------------
    // Build graph command
    // -------------------------------------------------------------------------
    // Colors: Control = gray (gs8), Treated = teal (#1E88A8 = "30 136 168")
    
    local graph_cmd ""
    
    if "`ci'" != "" {
        // CI bands with 30% transparency
        local graph_cmd `graph_cmd' (rarea outcome_ci_lb outcome_ci_ub id_time_std if Gi==0, color(gs8%30) lwidth(none))
        local graph_cmd `graph_cmd' (rarea outcome_ci_lb outcome_ci_ub id_time_std if Gi==1, color("30 136 168"%30) lwidth(none))
    }
    
    // Main trend lines with point markers
    local graph_cmd `graph_cmd' (connected outcome_mean id_time_std if Gi==0, lcolor(gs8) mcolor(gs8) lpattern(solid) msymbol(O))
    local graph_cmd `graph_cmd' (connected outcome_mean id_time_std if Gi==1, lcolor("30 136 168") mcolor("30 136 168") lpattern(solid) msymbol(O))
    
    // -------------------------------------------------------------------------
    // Configure legend
    // -------------------------------------------------------------------------
    local legend_order ""
    if "`ci'" != "" {
        // With CI bands, legend references line elements (3 and 4)
        local legend_order `"order(3 "Control" 4 "Treated") title("Group")"'
    }
    else {
        local legend_order `"order(1 "Control" 2 "Treated") title("Group")"'
    }
    
    // -------------------------------------------------------------------------
    // Generate graph
    // -------------------------------------------------------------------------
    local title_opt ""
    if "`title'" != "" {
        local title_opt `"title("`title'")"'
    }
    
    twoway `graph_cmd', ///
        xline(0, lpattern(dash) lcolor(black)) ///
        legend(`legend_order' rows(1)) ///
        `title_opt' ///
        xtitle("`xtitle'") ///
        ytitle("`ytitle'") ///
        `scheme_opt' ///
        name(`name', replace)
    
    restore
    
end



// =============================================================================
// _plot_placebo()
// Generate placebo plot showing 95% standardized equivalence confidence intervals
//
// Displays error bars for pre-treatment period estimates to assess the
// parallel trends assumption. The 95% equivalence CI is symmetric around zero.
//
// Equivalence CI calculation:
//   EqCI95_UB = max(|estimate + z_{0.95} * se|, |estimate - z_{0.95} * se|)
//   EqCI95_LB = -EqCI95_UB
//
// Arguments:
//   saving  : string - Output file path
//   scheme  : string - Graph scheme name
//   title   : string - Graph title
//   xtitle  : string - X-axis title (default: "Time relative to treatment assignment")
//   ytitle  : string - Y-axis title (default: "95% Standardized Equivalence CI")
//   name    : string - Graph name in memory (default: "placebo_plot")
//
// Data source:
//   e(placebo) matrix with columns: lag, estimate, std_error, estimate_orig,
//   std_error_orig, EqCI95_LB, EqCI95_UB
// =============================================================================

program define _plot_placebo
    version 16.0
    
    syntax [, SAVing(string) SCHeme(string) ///
              TItle(string) XTItle(string) YTItle(string) ///
              name(string)]
    
    // -------------------------------------------------------------------------
    // Validate e(placebo) exists
    // -------------------------------------------------------------------------
    if "`e(cmd)'" == "" {
        display as error "no estimation results found"
        exit 301
    }
    
    capture confirm matrix e(placebo)
    if _rc {
        display as error "e(placebo) not found; run diddesign_check first"
        exit 301
    }
    
    // -------------------------------------------------------------------------
    // Extract data to temporary dataset
    // -------------------------------------------------------------------------
    // e(placebo) columns: lag, estimate, std_error, estimate_orig,
    // std_error_orig, EqCI95_LB, EqCI95_UB
    preserve
    clear
    
    tempname placebo
    matrix `placebo' = e(placebo)
    local nrows = rowsof(`placebo')
    local ncols = colsof(`placebo')
    
    if `nrows' == 0 {
        display as error "e(placebo) is empty"
        restore
        exit 198
    }
    
    // Validate expected matrix structure
    if `ncols' < 7 {
        display as error "e(placebo) has unexpected structure (expected 7 columns, found `ncols')"
        display as error "Please re-run diddesign_check"
        restore
        exit 198
    }
    
    quietly {
        set obs `nrows'
        gen double lag = .
        gen double estimate = .
        gen double std_error = .
        gen double EqCI95_LB = .
        gen double EqCI95_UB = .
        
        forvalues i = 1/`nrows' {
            replace lag = `placebo'[`i', 1] in `i'
            replace estimate = `placebo'[`i', 2] in `i'
            replace std_error = `placebo'[`i', 3] in `i'
            // Pre-computed EqCI95 values (columns 6-7)
            replace EqCI95_LB = `placebo'[`i', 6] in `i'
            replace EqCI95_UB = `placebo'[`i', 7] in `i'
        }
    }
    
    // -------------------------------------------------------------------------
    // Calculate time relative to treatment (time_to_treat = -lag)
    // -------------------------------------------------------------------------
    gen double time_to_treat = -lag
    
    // -------------------------------------------------------------------------
    // Set default values
    // -------------------------------------------------------------------------
    if `"`title'"' == "" local title ""
    if `"`xtitle'"' == "" local xtitle "Time relative to treatment assignment"
    if `"`ytitle'"' == "" local ytitle "95% Standardized Equivalence CI"
    if "`name'" == "" local name "placebo_plot"
    
    // -------------------------------------------------------------------------
    // Handle single lag case (expand x-axis range)
    // -------------------------------------------------------------------------
    local xlim_opt ""
    if `nrows' == 1 {
        quietly sum lag
        local lag_val = r(mean)
        local xlim_lb = -`lag_val' - 1
        local xlim_ub = -`lag_val' + 1
        local xlim_opt "xscale(range(`xlim_lb' `xlim_ub'))"
    }
    
    // -------------------------------------------------------------------------
    // Build options
    // -------------------------------------------------------------------------
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt "scheme(`scheme')"
    }
    
    local title_opt ""
    if `"`title'"' != "" {
        local title_opt `"title(`"`title'"')"'
    }
    
    // -------------------------------------------------------------------------
    // Generate placebo plot
    // -------------------------------------------------------------------------
    // Error bars only (no point markers)
    // Color: teal (#1E88A8 = "30 136 168")
    // Reference line at y=0 (dotted, gray)
    
    twoway (rcap EqCI95_LB EqCI95_UB time_to_treat, ///
                lcolor("30 136 168") msize(vsmall)), ///
           yline(0, lpattern(dot) lcolor(gs8)) ///
           `title_opt' ///
           xtitle(`"`xtitle'"') ///
           ytitle(`"`ytitle'"') ///
           legend(off) ///
           `xlim_opt' ///
           `scheme_opt' ///
           name(`name', replace)
    
    // -------------------------------------------------------------------------
    // Save graph if requested
    // -------------------------------------------------------------------------
    if "`saving'" != "" {
        _export_graph, saving("`saving'") `replace'
    }
    
    restore
    
end



// =============================================================================
// _plot_pattern()
// Generate treatment pattern heatmap for staggered adoption design
//
// Displays treatment status over time for each unit as a tile plot.
// Units are sorted by first treatment timing (earliest at top, never-treated
// at bottom).
//
// Arguments:
//   saving  : string - Output file path
//   scheme  : string - Graph scheme name
//   title   : string - Graph title
//   xtitle  : string - X-axis title (default: "Time")
//   ytitle  : string - Y-axis title (default: "Unit")
//   name    : string - Graph name in memory (default: "pattern_plot")
//
// Data source:
//   e(Gmat) matrix from diddesign_check with design(sa)
//   Rows: units, Columns: time periods
//   Values: 0=not treated, non-zero=treated
// =============================================================================

program define _plot_pattern
    version 16.0
    
    syntax [, SAVing(string) SCHeme(string) ///
              TItle(string) XTItle(string) YTItle(string) ///
              name(string)]
    
    // -------------------------------------------------------------------------
    // Validate e(Gmat) exists and design is SA
    // -------------------------------------------------------------------------
    if "`e(cmd)'" == "" {
        display as error "no estimation results found"
        exit 301
    }
    
    if "`e(design)'" != "sa" {
        display as error "pattern plot only available for SA design"
        exit 198
    }
    
    capture confirm matrix e(Gmat)
    if _rc {
        display as error "e(Gmat) not found; run diddesign_check with design(sa) first"
        exit 301
    }
    
    // -------------------------------------------------------------------------
    // Extract Gmat matrix
    // -------------------------------------------------------------------------
    preserve
    clear
    
    tempname Gmat
    matrix `Gmat' = e(Gmat)
    local n_units = rowsof(`Gmat')
    local n_times = colsof(`Gmat')
    
    if `n_units' == 0 | `n_times' == 0 {
        display as error "e(Gmat) is empty"
        restore
        exit 198
    }
    
    // Check for placeholder matrix indicating failed analysis
    if `n_units' == 1 & `n_times' == 1 {
        if missing(`Gmat'[1,1]) {
            display as error "e(Gmat) is a placeholder matrix (1x1 with missing value)"
            display as error "SA analysis may have failed. Please check diddesign_check output."
            display as error "Tip: Verify that your data has sufficient treated units meeting the threshold."
            restore
            exit 198
        }
    }
    
    // -------------------------------------------------------------------------
    // Compute treatment timing and sort order
    // -------------------------------------------------------------------------
    // Units are sorted by first treatment period:
    //   - Never-treated units: Y = 1 (bottom)
    //   - Latest-treated units: middle positions
    //   - Earliest-treated units: Y = n_units (top)
    
    tempname sorted_pos
    mata: _compute_treat_timing("`Gmat'")
    matrix `sorted_pos' = r(sort_order)
    
    // -------------------------------------------------------------------------
    // Convert to long format dataset
    // -------------------------------------------------------------------------
    quietly {
        local n_obs = `n_units' * `n_times'
        set obs `n_obs'
        
        gen long id_subject = .
        gen long id_time = .
        gen byte treated = .
        gen long id_subject_sorted = .
        
        local obs = 1
        forvalues i = 1/`n_units' {
            forvalues t = 1/`n_times' {
                replace id_subject = `i' in `obs'
                replace id_time = `t' in `obs'
                // Binary treatment indicator: 0=control, non-zero=treated
                replace treated = (`Gmat'[`i', `t'] != 0) in `obs'
                replace id_subject_sorted = `sorted_pos'[`i', 1] in `obs'
                local obs = `obs' + 1
            }
        }
    }
    
    // -------------------------------------------------------------------------
    // Set default values
    // -------------------------------------------------------------------------
    if `"`title'"' == "" local title ""
    if `"`xtitle'"' == "" local xtitle "Time"
    if `"`ytitle'"' == "" local ytitle "Unit"
    if "`name'" == "" local name "pattern_plot"
    
    // -------------------------------------------------------------------------
    // Calculate tile size based on data dimensions
    // -------------------------------------------------------------------------
    if `n_units' <= 20 & `n_times' <= 20 {
        local tile_size = "vlarge"
    }
    else if `n_units' <= 50 & `n_times' <= 50 {
        local tile_size = "large"
    }
    else if `n_units' <= 100 & `n_times' <= 100 {
        local tile_size = "medium"
    }
    else {
        local tile_size = "small"
    }
    
    // -------------------------------------------------------------------------
    // Build options
    // -------------------------------------------------------------------------
    local title_opt ""
    if `"`title'"' != "" {
        local title_opt `"title(`"`title'"')"'
    }
    
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt "scheme(`scheme')"
    }
    
    // -------------------------------------------------------------------------
    // Generate pattern plot (heatmap)
    // -------------------------------------------------------------------------
    // Colors: control = light gray, treated = teal (#1E88A8 = "30 136 168")
    // Axis labels hidden for cleaner appearance
    
    twoway (scatter id_subject_sorted id_time if treated==0, ///
                msymbol(S) msize(`tile_size') mcolor(ltgray)) ///
           (scatter id_subject_sorted id_time if treated==1, ///
                msymbol(S) msize(`tile_size') mcolor("30 136 168")), ///
           `title_opt' ///
           xtitle("`xtitle'") ///
           ytitle("`ytitle'") ///
           legend(order(1 "control" 2 "treated") title("Status") rows(1)) ///
           xlabel(, nolabels noticks) ylabel(, nolabels) ///
           `scheme_opt' ///
           name(`name', replace)
    
    // -------------------------------------------------------------------------
    // Save graph if requested
    // -------------------------------------------------------------------------
    if "`saving'" != "" {
        _export_graph, saving("`saving'") `replace'
    }
    
    restore
    
end


// =============================================================================
// _plot_estimates()
// Generate Double-DID estimates plot with confidence intervals
//
// Displays Double-DID point estimates across lead values with 90% CI.
// Optionally overlays placebo estimates from stored diddesign_check results.
//
// Arguments:
//   design    : string - Design type ("did" or "sa")
//   scheme    : string - Graph scheme name
//   title     : string - Graph title
//   xtitle    : string - X-axis title (default: "Time")
//   ytitle    : string - Y-axis title (default: "Estimates (90% CI)")
//   name      : string - Graph name in memory (default: "estimates_plot")
//   use_check : string - Name of stored diddesign_check results to overlay
//   band      : flag   - Display CI as ribbon instead of error bars
//
// Data source:
//   e(estimates) matrix from diddesign
//   Row structure: 3 rows per lead (Double-DID, DID, sequential DID)
//   Columns: lead, estimate, std_error, ci_lo, ci_hi, weight
// =============================================================================

program define _plot_estimates
    version 16.0
    
    syntax [, DESIGN(string) SCHeme(string) ///
              TItle(string) XTItle(string) YTItle(string) ///
              name(string) USE_CHECK(string) BAND]
    
    // -------------------------------------------------------------------------
    // Validate e(estimates) exists
    // -------------------------------------------------------------------------
    if "`e(cmd)'" == "" {
        display as error "no estimation results found"
        exit 301
    }
    
    capture confirm matrix e(estimates)
    if _rc {
        display as error "e(estimates) not found; run diddesign first"
        exit 301
    }
    
    // -------------------------------------------------------------------------
    // Extract Double-DID estimates
    // -------------------------------------------------------------------------
    // e(estimates) contains 3 rows per lead: Double-DID, DID, sequential DID
    // Extract only Double-DID rows (every 3rd row starting from 1)
    
    preserve
    clear
    
    tempname estimates
    matrix `estimates' = e(estimates)
    local nrows = rowsof(`estimates')
    local ncols = colsof(`estimates')
    
    if `nrows' == 0 {
        display as error "e(estimates) is empty"
        restore
        exit 198
    }
    
    local rownames : rownames `estimates'
    local n_ddid = ceil(`nrows' / 3)
    
    quietly {
        set obs `n_ddid'
        gen double time = .
        gen double estimate = .
        gen double std_error = .
        gen double CI90_LB = .
        gen double CI90_UB = .
        gen byte source = 1  // 1 = diddesign estimates
        
        // Extract Double-DID rows
        local obs = 1
        forvalues i = 1(3)`nrows' {
            replace time = `estimates'[`i', 1] in `obs'
            replace estimate = `estimates'[`i', 2] in `obs'
            replace std_error = `estimates'[`i', 3] in `obs'
            // 90% CI: estimate +/- z_{0.95} * std_error
            replace CI90_LB = `estimates'[`i', 2] - invnormal(0.95) * `estimates'[`i', 3] in `obs'
            replace CI90_UB = `estimates'[`i', 2] + invnormal(0.95) * `estimates'[`i', 3] in `obs'
            local obs = `obs' + 1
        }
    }
    
    // -------------------------------------------------------------------------
    // Overlay placebo estimates (if use_check specified)
    // -------------------------------------------------------------------------
    // Placebo estimates are placed at time = -lag (negative lag values)
    
    if "`use_check'" != "" {
        capture estimates describe `use_check'
        if _rc {
            display as error "stored results '`use_check'' not found"
            display as error "use 'estimates store' to save diddesign_check results first"
            restore
            exit 301
        }
        
        // Temporarily restore stored results to access e(placebo)
        tempname current_est
        _estimates hold `current_est', restore
        
        quietly estimates restore `use_check'
        
        if "`e(cmd)'" != "diddesign_check" {
            _estimates unhold `current_est'
            display as error "'`use_check'' is not diddesign_check results"
            restore
            exit 198
        }
        
        capture confirm matrix e(placebo)
        if _rc {
            _estimates unhold `current_est'
            display as error "e(placebo) not found in '`use_check''"
            restore
            exit 301
        }
        
        tempname placebo
        matrix `placebo' = e(placebo)
        local n_placebo = rowsof(`placebo')
        
        _estimates unhold `current_est'
        
        // Add placebo data using original (non-standardized) estimates
        local new_obs = _N + `n_placebo'
        quietly {
            set obs `new_obs'
            
            local start = _N - `n_placebo' + 1
            forvalues i = 1/`n_placebo' {
                local obs = `start' + `i' - 1
                replace time = -`placebo'[`i', 1] in `obs'
                replace estimate = `placebo'[`i', 4] in `obs'
                replace std_error = `placebo'[`i', 5] in `obs'
                replace CI90_LB = `placebo'[`i', 4] - invnormal(0.95) * `placebo'[`i', 5] in `obs'
                replace CI90_UB = `placebo'[`i', 4] + invnormal(0.95) * `placebo'[`i', 5] in `obs'
                replace source = 0 in `obs'  // 0 = placebo estimates
            }
        }
        
        sort time
    }
    
    // -------------------------------------------------------------------------
    // Set default values
    // -------------------------------------------------------------------------
    if `"`title'"' == "" local title ""
    if `"`xtitle'"' == "" local xtitle "Time"
    if `"`ytitle'"' == "" local ytitle "Estimates (90% CI)"
    if "`name'" == "" local name "estimates_plot"
    
    // -------------------------------------------------------------------------
    // Handle single time point case (expand x-axis range)
    // -------------------------------------------------------------------------
    local xlim_opt ""
    quietly count
    if r(N) == 1 {
        quietly sum time
        local time_val = r(mean)
        local xlim_lb = `time_val' - 1
        local xlim_ub = `time_val' + 1
        local xlim_opt "xscale(range(`xlim_lb' `xlim_ub'))"
    }
    
    // -------------------------------------------------------------------------
    // Build options
    // -------------------------------------------------------------------------
    local scheme_opt ""
    if "`scheme'" != "" {
        local scheme_opt "scheme(`scheme')"
    }
    
    local title_opt ""
    if `"`title'"' != "" {
        local title_opt `"title(`"`title'"')"'
    }
    
    // -------------------------------------------------------------------------
    // Generate estimates plot
    // -------------------------------------------------------------------------
    // X-axis labels at unique time values
    quietly levelsof time, local(time_values)
    local xlabel_opt "xlabel(`time_values')"
    
    if "`band'" != "" {
        // Band mode: CI ribbon + connected line with points
        twoway (rarea CI90_LB CI90_UB time, color(gs8%50) lwidth(none)) ///
               (connected estimate time, lcolor(black) mcolor(black) ///
                    lpattern(solid) msymbol(O)), ///
               yline(0, lpattern(dash) lcolor(gs8)) ///
               `title_opt' ///
               xtitle(`"`xtitle'"') ///
               ytitle(`"`ytitle'"') ///
               `xlabel_opt' ///
               `xlim_opt' ///
               legend(off) ///
               `scheme_opt' ///
               name(`name', replace)
    }
    else {
        // Error bar mode (default): error bars + scatter points
        twoway (rcap CI90_LB CI90_UB time, lcolor(black) msize(vsmall)) ///
               (scatter estimate time, mcolor(black) msymbol(O)), ///
               yline(0, lpattern(dash) lcolor(gs8)) ///
               `title_opt' ///
               xtitle(`"`xtitle'"') ///
               ytitle(`"`ytitle'"') ///
               `xlabel_opt' ///
               `xlim_opt' ///
               legend(off) ///
               `scheme_opt' ///
               name(`name', replace)
    }
    
    restore
    
end
