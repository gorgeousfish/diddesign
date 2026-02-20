*! _diddesign_display_table_header.ado - Display column headers for estimation results
*!
*! Outputs Stata-style coefficient table header with depvar abbreviation.
*! Table width is 78 chars to fit standard Stata linesize.

version 16.0

program define _diddesign_display_table_header
    version 16.0
    
    syntax , [DEPVAR(string) LEVEL(real 95) WEIGHT]
    
    // Abbreviate depvar to fit 13 chars
    if "`depvar'" != "" {
        local depvar_abbrev = abbrev("`depvar'", 13)
    }
    else {
        local depvar_abbrev ""
    }
    
    // Top separator
    display as text ""
    
    if "`weight'" != "" {
        // Table with Weight column (78 chars)
        display as text "{hline 78}"
        display as text %13s "`depvar_abbrev'" " |" ///
                "  Coef.   Std. err." ///
                "     z   P>|z|" ///
                "   [`level'% conf. interval]" ///
                "  Weight"
        display as text "{hline 14}+{hline 63}"
    }
    else {
        // Standard table without Weight column (69 chars)
        display as text "{hline 69}"
        display as text %13s "`depvar_abbrev'" " |" ///
                "  Coef.   Std. err." ///
                "     z   P>|z|" ///
                "   [`level'% conf. interval]"
        display as text "{hline 14}+{hline 54}"
    }
end
