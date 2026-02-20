*! _diddesign_display_footer.ado - Display estimation footer note
*!
*! Displays a single-line note after the results table.

version 16.0

program define _diddesign_display_footer
    version 16.0
    
    syntax , [N_boot(integer 30) LEVEL(real 95) ///
              CI_method(string) NOTES(string)]
    
    display as text ""
    display as text "Note: Double DID combines DID and sDID via optimal GMM weights."
    
end
