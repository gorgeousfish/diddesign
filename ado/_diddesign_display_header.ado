*! _diddesign_display_header.ado - Display estimation header in Stata official style
*!
*! Renders a two-column header: left side shows title/design info,
*! right side shows sample statistics aligned with "=" signs.

version 16.0

program define _diddesign_display_header
    version 16.0
    
    syntax , CMD(string) DESIGN(string) DATATYPE(string) ///
             N(integer) [N_units(integer 0) N_periods(integer 0) ///
             N_boot(integer 30) CLuster(string) THRes(integer 2)]
    
    // -------------------------------------------------------------------------
    // Determine labels
    // -------------------------------------------------------------------------
    
    if "`cmd'" == "diddesign" {
        if "`design'" == "sa" {
            local title "Staggered Adoption Double DID"
        }
        else {
            local title "Double Difference-in-Differences"
        }
    }
    else if "`cmd'" == "diddesign_check" {
        local title "Parallel Trends Assessment"
    }
    else {
        local title "DIDdesign"
    }
    
    if "`design'" == "sa" {
        local design_label "Staggered Adoption"
    }
    else {
        local design_label "Standard DID"
    }
    
    // -------------------------------------------------------------------------
    // Line 1: Title + Number of obs
    // -------------------------------------------------------------------------
    
    display as text _n "`title'" _col(41) "Number of obs" _col(63) "=" ///
        _col(65) as result %10.0fc `n'
    
    // -------------------------------------------------------------------------
    // Line 2: Design + second stat
    // -------------------------------------------------------------------------
    
    if "`cmd'" == "diddesign" {
        if `n_units' > 0 {
            display as text "Design:" _col(16) as result "`design_label'" ///
                as text _col(41) "Number of units" _col(63) "=" ///
                _col(65) as result %10.0fc `n_units'
        }
        else {
            display as text "Design:" _col(16) as result "`design_label'"
        }
    }
    else {
        local n_clusters = e(n_clusters)
        display as text "Design:" _col(16) as result "`design_label'" ///
            as text _col(41) "Number of clusters" _col(63) "=" ///
            _col(65) as result %10.0fc `n_clusters'
    }
    
    // -------------------------------------------------------------------------
    // Line 3: Data type + third stat
    // -------------------------------------------------------------------------
    
    if "`cmd'" == "diddesign" {
        if `n_periods' > 0 {
            display as text "Data type:" _col(16) as result "`datatype'" ///
                as text _col(41) "Number of periods" _col(63) "=" ///
                _col(65) as result %10.0fc `n_periods'
        }
        else {
            display as text "Data type:" _col(16) as result "`datatype'"
        }
    }
    else {
        display as text "Data type:" _col(16) as result "`datatype'" ///
            as text _col(41) "Bootstrap reps" _col(63) "=" ///
            _col(65) as result %10.0fc `n_boot'
    }
    
    // -------------------------------------------------------------------------
    // Line 4: Bootstrap reps (diddesign) or Standardized + Bootstrap valid (check)
    // -------------------------------------------------------------------------
    
    if "`cmd'" == "diddesign" {
        display as text _col(41) "Bootstrap reps" _col(63) "=" ///
            _col(65) as result %10.0fc `n_boot'
    }
    else {
        local n_boot_valid = e(n_boot_valid)
        display as text "Standardized:" _col(16) as result "Yes" ///
            as text _col(41) "Bootstrap valid" _col(63) "=" ///
            _col(65) as result %10.0fc `n_boot_valid'
    }
    
    // -------------------------------------------------------------------------
    // Line 5: Equiv CI (check only) or blank
    // -------------------------------------------------------------------------
    
    if "`cmd'" == "diddesign_check" {
        display as text "Equiv. CI:" _col(16) as result "95%"
    }
    
    // -------------------------------------------------------------------------
    // Clustering line (if applicable)
    // -------------------------------------------------------------------------
    
    if "`cluster'" != "" {
        display as text ""
        display as text "Clustering:" _col(16) as result "`cluster'"
    }
    
end
