{smcl}
{* *! version 1.0.0  16jan2026}{...}
{vieweralsosee "diddesign" "help diddesign"}{...}
{vieweralsosee "diddesign_check" "help diddesign_check"}{...}
{vieweralsosee "diddesign_plot" "help diddesign_plot"}{...}
{vieweralsosee "diddesign_intro" "help diddesign_intro"}{...}
{viewerjumpto "Syntax" "diddesign##syntax"}{...}
{viewerjumpto "Description" "diddesign##description"}{...}
{viewerjumpto "Options" "diddesign##options"}{...}
{viewerjumpto "Stored results" "diddesign##results"}{...}
{viewerjumpto "Methods and formulas" "diddesign##methods"}{...}
{viewerjumpto "References" "diddesign##references"}{...}

{title:Title}

{phang}
{bf:diddesign} {hline 2} Double Difference-in-Differences Estimation

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:diddesign}
{depvar}
[{it:covariates}]
{ifin}
{cmd:,}
{opth treatment(varname)}
{opth time(varname)}
[{it:options}]

{pstd}
where {it:covariates} specifies optional control variables.
Alternatively, use the {opt covariates()} option.

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth treatment(varname)}}treatment indicator variable{p_end}
{synopt:{opth time(varname)}}time variable{p_end}

{syntab:Model}
{synopt:{opth covariates(string)}}optional control variables; supports factor variables{p_end}

{syntab:Design}
{synopt:{opt design(did|sa)}}estimation design; default is {cmd:did}{p_end}
{synopt:{opt panel}}panel data format (default){p_end}
{synopt:{opt rcs}}repeated cross-section data format{p_end}

{syntab:Panel Options}
{synopt:{opth id(varname)}}unit identifier variable (required for panel){p_end}

{syntab:RCS Options}
{synopt:{opth post(varname)}}post-treatment indicator (required for RCS){p_end}

{syntab:Inference}
{synopt:{opth cluster(varname)}}cluster variable for standard errors{p_end}
{synopt:{opt nboot(#)}}number of bootstrap replications; default is {cmd:30}{p_end}
{synopt:{opt seboot}}use bootstrap critical values for confidence intervals{p_end}
{synopt:{opt seed(#)}}random number seed for reproducibility{p_end}

{syntab:Staggered Adoption}
{synopt:{opt thres(#)}}minimum treated units threshold; default is {cmd:2}{p_end}
{synopt:{opt lead(numlist)}}lead periods for SA estimation; default is {cmd:0}{p_end}

{syntab:Reporting}
{synopt:{opt level(#)}}confidence level; default is {cmd:95}{p_end}
{synopt:{opt quiet}}suppress bootstrap progress display{p_end}
{synopt:{opt parallel}}enable parallel bootstrap{p_end}
{synoptline}

{p 4 6 2}
{it:covariates} supports factor-variable notation; see {help fvvarlist}.
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:diddesign} implements the Double Difference-in-Differences (Double DID) 
method from Egami and Yamauchi (2023). The method optimally combines the 
standard DID estimator with the sequential DID (s-DID) estimator using GMM 
weights to improve efficiency while maintaining robustness.

{pstd}
For {bf:standard DID designs} ({cmd:design(did)}), the command estimates both 
DID and s-DID, then computes the GMM-optimal weighted combination (Double DID).

{pstd}
For {bf:staggered adoption (SA) designs} ({cmd:design(sa)}), the command handles 
multiple treatment timing groups and computes time-weighted averages across 
cohorts.

{pstd}
The Double DID estimator is:

{p 8 8 2}
tau_dDID = w1 * tau_DID + w2 * tau_sDID

{pstd}
where w1 and w2 are GMM-optimal weights that minimize the asymptotic variance,
and tau_sDID denotes the sequential DID (s-DID) estimator.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth treatment(varname)} specifies the binary treatment indicator variable.
For panel data, the treatment variable indicates whether unit i is treated at 
time t (1 = treated, 0 = control). For RCS data, the treatment variable 
indicates treatment group membership.

{phang}
{opth time(varname)} specifies the time variable (e.g., year). String 
variables are automatically encoded to numeric using {cmd:egen group()}.

{dlgtab:Model}

{phang}
{it:covariates} specified directly after the dependent variable in the varlist
are optional control variables (e.g., {cmd:diddesign y x1 x2 x3, ...}).
Alternatively, the {opt covariates()} option below can be used to specify 
control variables.

{phang}
{opth covariates(fvvarlist)} specifies additional control variables to include 
in the estimation. This option has no limit on the number of variables and 
can be combined with inline covariates. These variables are included as 
covariates in the difference-in-differences regression to improve precision 
and control for confounding factors.

{pmore}
Factor-variable notation is supported; see {help fvvarlist}. For example, 
{cmd:covariates(i.region)} automatically expands a categorical variable 
{cmd:region} into dummy variables. Factor variables are automatically 
expanded using {cmd:fvrevar} before estimation, and the number of expanded 
variables is displayed.

{pmore}
Examples of factor-variable syntax:

{p 8 12 2}{cmd:i.region} - indicator (dummy) variables for categorical region{p_end}
{p 8 12 2}{cmd:ib3.region} - same as above, but with region=3 as base{p_end}
{p 8 12 2}{cmd:ibn.region} - no base category (include all levels){p_end}
{p 8 12 2}{cmd:c.x1#c.x2} - interaction of continuous variables{p_end}
{p 8 12 2}{cmd:i.region#c.gdp} - interaction of factor and continuous{p_end}

{dlgtab:Design}

{phang}
{opt design(did|sa)} specifies the estimation design. {cmd:did} (default) 
uses the standard DID design with a single treatment time. {cmd:sa} uses 
the staggered adoption design for multiple treatment times across units.

{phang}
{opt panel} indicates that the data is in panel format where the same units 
are observed over time. This is the default when {opt id()} is specified.

{phang}
{opt rcs} indicates that the data is in repeated cross-section format where
different units are observed at each time point. When using RCS data, 
the {opt post()} option is required.

{dlgtab:Panel Options}

{phang}
{opth id(varname)} specifies the variable identifying units (individuals, 
firms, states, etc.). This option is required for panel data. String variables 
are automatically encoded to numeric using {cmd:egen group()}.

{dlgtab:RCS Options}

{phang}
{opth post(varname)} specifies the post-treatment period indicator variable
for repeated cross-section data. This variable should equal 1 for post-treatment
periods and 0 for pre-treatment periods. This option is required when using
RCS data.

{dlgtab:Inference}

{phang}
{opth cluster(varname)} specifies the variable for cluster-robust bootstrap 
standard errors. If not specified for panel data, clustering is done at the 
unit level (same as {cmd:id()}). String variables are automatically encoded 
to numeric.

{phang}
{opt nboot(#)} specifies the number of bootstrap replications for variance 
estimation. The default is 30. For publication-quality results, consider 
using 500 or more replications.

{phang}
{opt seboot} requests that bootstrap-based critical values be used to 
construct confidence intervals instead of asymptotic normal critical values.
This option only affects the standard DID design.

{phang}
{opt seed(#)} sets the random number seed for reproducibility of bootstrap 
results.

{dlgtab:Staggered Adoption}

{phang}
{opt thres(#)} specifies the minimum number of treated units required for a 
time period to be included in the SA analysis. The default is 2. Periods with 
fewer treated units are excluded from the estimation.

{phang}
{opt lead(numlist)} specifies the lead periods for SA estimation. The default 
is 0, which estimates the instantaneous treatment effect. Positive values 
estimate lagged effects.

{dlgtab:Reporting}

{phang}
{opt level(#)} specifies the confidence level, as a percentage, for confidence 
intervals. The default is 95.

{phang}
{opt quiet} suppresses the bootstrap progress display. This is useful when 
running the command in batch mode or when progress messages are not needed.

{phang}
{opt parallel} is reserved for future use; currently has no effect.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:diddesign} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_units)}}number of units{p_end}
{synopt:{cmd:e(n_periods)}}number of time periods{p_end}
{synopt:{cmd:e(n_boot)}}number of bootstrap replications{p_end}
{synopt:{cmd:e(n_lead)}}number of lead periods{p_end}
{synopt:{cmd:e(level)}}confidence level{p_end}
{synopt:{cmd:e(is_panel)}}1 if panel data, 0 otherwise{p_end}
{synopt:{cmd:e(n_boot_success)}}number of successful bootstrap iterations{p_end}
{synopt:{cmd:e(seboot)}}1 if bootstrap CI used, 0 otherwise{p_end}

{p2col 5 20 24 2: Scalars (SA design only)}{p_end}
{synopt:{cmd:e(n_periods_valid)}}number of valid time periods meeting threshold{p_end}
{synopt:{cmd:e(thres)}}minimum treated units threshold used{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:diddesign}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(design)}}estimation design ({cmd:did} or {cmd:sa}){p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(treatment)}}name of treatment variable(s){p_end}
{synopt:{cmd:e(covariates)}}names of covariates{p_end}
{synopt:{cmd:e(id)}}name of unit identifier variable{p_end}
{synopt:{cmd:e(time)}}name of time variable{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(ci_method)}}confidence interval method ({cmd:bootstrap} or {cmd:asymptotic}){p_end}
{synopt:{cmd:e(lead)}}lead values as specified{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}

{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector with columns named {cmd:dDID:lead_}{it:k}{cmd:},
{cmd:DID:lead_}{it:k}{cmd:}, {cmd:sDID:lead_}{it:k} for standard DID design;
{cmd:SA_dDID:lead_}{it:k}{cmd:}, {cmd:SA_DID:lead_}{it:k}{cmd:}, 
{cmd:SA_sDID:lead_}{it:k} for SA design{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix{p_end}
{synopt:{cmd:e(estimates)}}estimation results matrix with columns: lead, estimate, 
std_error, ci_lo, ci_hi, weight{p_end}
{synopt:{cmd:e(lead_values)}}lead values vector{p_end}
{synopt:{cmd:e(weights)}}GMM weights matrix with columns: w_did, w_sdid{p_end}
{synopt:{cmd:e(W)}}GMM optimal weight matrix (2 x 2){p_end}
{synopt:{cmd:e(vcov_gmm)}}variance-covariance matrix of moment conditions (2 x 2){p_end}

{p2col 5 20 24 2: Matrices (SA design only)}{p_end}
{synopt:{cmd:e(time_weights)}}time weights vector (proportional to treated units per period){p_end}

{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}

{marker methods}{...}
{title:Methods and formulas}

{pstd}
The Double DID method optimally combines the standard DID estimator with the 
sequential DID (s-DID) estimator using the generalized method of moments (GMM).

{pstd}
{bf:Identification assumptions:}

{p 8 8 2}
The standard DID estimator (tau_DID) requires the {it:parallel trends assumption}: 
the outcome trend of the control group would have been the same as the trend 
of the treatment group in the absence of treatment.

{p 8 8 2}
The sequential DID estimator (tau_sDID) requires only the weaker 
{it:parallel trends-in-trends assumption}: the change in outcome trends is the 
same across treatment and control groups.

{pstd}
{bf:Estimation:}

{p 8 8 2}
tau_DID uses all pre-treatment periods as the baseline.

{p 8 8 2}
tau_sDID uses only the period immediately before treatment as the baseline.

{pstd}
The Double DID estimator combines both estimators via GMM:

{p 8 8 2}
tau_dDID = w1 * tau_DID + w2 * tau_sDID

{pstd}
where w1 + w2 = 1 and the GMM optimal weight matrix W is the inverse of the 
variance-covariance matrix of the two estimators. This yields the smallest 
asymptotic variance among all weighted combinations. See Egami and Yamauchi 
(2023, equations 12-14) for details.

{pstd}
For staggered adoption designs, the SA-Double-DID estimator computes 
time-weighted averages across treatment cohorts, with weights proportional to 
the number of newly treated units in each period.

{marker references}{...}
{title:References}

{phang}
Egami, N. and S. Yamauchi. 2023.
Using Multiple Pretreatment Periods to Improve Difference-in-Differences
and Staggered Adoption Designs.
{it:Political Analysis} 31(2): 195-212.
{browse "https://doi.org/10.1017/pan.2022.8"}
{p_end}

{title:Author}

{pstd}
Xuanyu Cai{break}
City University of Macau{break}
xuanyuCAI@outlook.com

{pstd}
Wenli Xu{break}
City University of Macau{break}
wlxu@cityu.edu.mo

{title:Also see}

{psee}
Online: {helpb diddesign}, {helpb diddesign_check}, {helpb diddesign_plot}, {helpb diddesign_intro}
{p_end}
