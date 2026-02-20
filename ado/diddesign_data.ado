*! version 1.0.0  20feb2026
*! diddesign_data - Load example datasets for the diddesign package

program define diddesign_data
    version 16
    syntax [anything(name=dataname)] [, CLEAR LIST DESCRIBE]

    // -----------------------------------------------------------------------
    // Dataset registry
    // -----------------------------------------------------------------------
    local valid_datasets "anzia2012 malesky2014 paglayan2019"

    local desc_anzia2012    `"Anzia (2012): Panel data, Texas school districts (2003-2009)"'
    local desc_malesky2014  `"Malesky et al. (2014): RCS data, Vietnam communes (2006, 2008, 2010)"'
    local desc_paglayan2019 `"Paglayan (2019): Staggered adoption, US states (1959-2000)"'

    local nobs_anzia2012    "6,377"
    local nobs_malesky2014  "6,269"
    local nobs_paglayan2019 "2,058"

    local nvars_anzia2012    "11"
    local nvars_malesky2014  "42"
    local nvars_paglayan2019 "9"

    local cite_anzia2012    "Anzia (2012), AJPS"
    local cite_malesky2014  "Malesky, Nguyen & Tran (2014), APSR"
    local cite_paglayan2019 "Paglayan (2019), AJPS"

    // -----------------------------------------------------------------------
    // Option: list - display all available datasets
    // -----------------------------------------------------------------------
    if "`list'" != "" {
        display as text ""
        display as text "{hline 78}"
        display as text "  Available datasets for {bf:diddesign_data}"
        display as text "{hline 78}"
        display as text ""
        display as text "  {bf:Dataset}        {bf:Obs}      {bf:Vars}   {bf:Description}"
        display as text "  {hline 72}"
        foreach dname of local valid_datasets {
            display as result "  `dname'" _col(22) as text "`nobs_`dname''" _col(33) as text "`nvars_`dname''" _col(40) as text "`desc_`dname''"
        }
        display as text ""
        display as text "  Source: `cite_anzia2012'"
        display as text "          `cite_malesky2014'"
        display as text "          `cite_paglayan2019'"
        display as text "{hline 78}"
        display as text ""
        display as text `"  Usage: {bf:diddesign_data {it:dataname} [, clear]}"'
        display as text `"         {bf:diddesign_data {it:dataname} , describe}"'
        display as text ""
        exit
    }

    // -----------------------------------------------------------------------
    // Validate: must specify a dataset name if not using list
    // -----------------------------------------------------------------------
    if `"`dataname'"' == "" {
        display as error "must specify a dataset name or the {bf:list} option"
        display as error ""
        display as error "Usage:"
        display as error `"  {bf:diddesign_data {it:dataname} [, clear]}"'
        display as error `"  {bf:diddesign_data , list}"'
        display as error `"  {bf:diddesign_data {it:dataname} , describe}"'
        display as error ""
        display as error "Type {bf:help diddesign_data} for details."
        exit 198
    }

    // -----------------------------------------------------------------------
    // Validate dataset name
    // -----------------------------------------------------------------------
    local name_valid = 0
    foreach dname of local valid_datasets {
        if "`dataname'" == "`dname'" {
            local name_valid = 1
        }
    }

    if `name_valid' == 0 {
        display as error `""`dataname'" is not a valid dataset"'
        display as error ""
        display as error "Available datasets: {bf:`valid_datasets'}"
        display as error ""
        display as error "Type {bf:diddesign_data, list} for details."
        exit 198
    }

    // -----------------------------------------------------------------------
    // Locate the data file using findfile
    // -----------------------------------------------------------------------
    capture findfile `dataname'.dta
    if _rc != 0 {
        display as error `"file "`dataname'.dta" not found"'
        display as error ""
        display as error "Please verify that the {bf:diddesign} package is properly installed."
        display as error "To install: {bf:net install diddesign, from(...)}"
        display as error "To get data files: {bf:net get diddesign}"
        exit 601
    }
    local filepath `"`r(fn)'"'

    // -----------------------------------------------------------------------
    // Option: describe - show dataset metadata without loading
    // -----------------------------------------------------------------------
    if "`describe'" != "" {
        display as text ""
        display as text "{hline 78}"
        display as text "  Dataset: {bf:`dataname'}"
        display as text "  `desc_`dataname''"
        display as text "{hline 78}"
        display as text ""
        describe using `"`filepath'"'
        exit
    }

    // -----------------------------------------------------------------------
    // Load the dataset
    // -----------------------------------------------------------------------
    if "`clear'" != "" {
        use `"`filepath'"', clear
    }
    else {
        use `"`filepath'"'
    }

    // Display confirmation message
    display as text ""
    display as result "  ({bf:`dataname'} loaded successfully)"
    display as text "  `desc_`dataname''"
    display as text "  `nobs_`dataname'' observations, `nvars_`dataname'' variables"
    display as text ""

end
