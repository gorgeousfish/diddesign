*! _diddesign_display_result.ado - Display a single row of estimation results
*!
*! Renders one row in Stata official coefficient table style with columns:
*! Coefficient, Std. err., z, P>|z|, [CI], and optional Weight.
*! Table width matches 78 chars (with weight) or 69 chars (without).

version 16.0

program define _diddesign_display_result
    
    syntax , LABEL(string) ESTIMATE(real) SE(real) ///
             CI_low(real) CI_high(real) [WEIGHT(real -1)]
    
    // -------------------------------------------------------------------------
    // Compute z-statistic and p-value
    // -------------------------------------------------------------------------
    
    if `se' > 0 & `se' < . {
        local z_val = `estimate' / `se'
        local p_val = 2 * (1 - normal(abs(`z_val')))
    }
    else {
        local z_val = .
        local p_val = .
    }
    
    // -------------------------------------------------------------------------
    // Display formatted result row
    // -------------------------------------------------------------------------
    
    if `weight' != -1 {
        // Row with Weight column (78 chars total)
        if `weight' == . {
            local wt_str "       ."
        }
        else {
            local wt_str = string(`weight', "%8.4f")
        }
        
        display as text %13s "`label'" " |" ///
                as result ///
                %9.4f `estimate' ///
                %10.4f `se' ///
                %7.2f `z_val' ///
                %7.3f `p_val' ///
                %11.4f `ci_low' ///
                %10.4f `ci_high' ///
                " `wt_str'"
    }
    else {
        // Row without Weight column (69 chars total)
        display as text %13s "`label'" " |" ///
                as result ///
                %9.4f `estimate' ///
                %10.4f `se' ///
                %7.2f `z_val' ///
                %7.3f `p_val' ///
                %11.4f `ci_low' ///
                %10.4f `ci_high'
    }
end
