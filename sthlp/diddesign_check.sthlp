{smcl}
{* *! version 1.1.0  29jan2026}{...}
{vieweralsosee "diddesign" "help diddesign"}{...}
{vieweralsosee "diddesign_check" "help diddesign_check"}{...}
{vieweralsosee "diddesign_plot" "help diddesign_plot"}{...}
{vieweralsosee "diddesign_intro" "help diddesign_intro"}{...}
{viewerjumpto "Syntax" "diddesign_check##syntax"}{...}
{viewerjumpto "Description" "diddesign_check##description"}{...}
{viewerjumpto "Options" "diddesign_check##options"}{...}
{viewerjumpto "Stored results" "diddesign_check##results"}{...}
{viewerjumpto "Methods and formulas" "diddesign_check##methods"}{...}
{viewerjumpto "References" "diddesign_check##references"}{...}

{title:Title}

{phang}
{bf:diddesign_check} {hline 2} Diagnostic Tests for Parallel Trends Assumption

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:diddesign_check}
{depvar}
[{it:covariates}]
{ifin}
{cmd:,}
{opth treatment(varname)}
{opth time(varname)}
{c -(}{opth id(varname)} | {opth post(varname)} {opt rcs}{c )-}
[{it:options}]

{pstd}
where {it:covariates} specifies optional control variables.

{synoptset 25 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Required}
{synopt:{opth treatment(varname)}}treatment indicator variable{p_end}
{synopt:{opth time(varname)}}time variable{p_end}

{syntab:Data Type (one required)}
{synopt:{opth id(varname)}}unit identifier for panel data{p_end}
{synopt:{opth post(varname)}}post-treatment indicator for RCS data{p_end}

{syntab:Design}
{synopt:{opt design(did|sa)}}design type; default is {cmd:did}{p_end}
{synopt:{opt panel}}panel data format (default if {opt id()} specified){p_end}
{synopt:{opt rcs}}repeated cross-section format (requires {opt post()}){p_end}

{syntab:Test Options}
{synopt:{opt lag(numlist)}}lag periods for placebo tests; default is {cmd:lag(1)}{p_end}
{synopt:{opt nboot(#)}}number of bootstrap iterations; default is {cmd:nboot(30)}{p_end}
{synopt:{opt thres(#)}}SA design threshold; default is {cmd:thres(2)}{p_end}

{syntab:Inference}
{synopt:{opth cluster(varname)}}cluster variable for standard errors{p_end}
{synopt:{opt seed(#)}}random number seed{p_end}

{syntab:Display}
{synopt:{opt quiet}}suppress bootstrap progress display{p_end}
{synopt:{opt parallel}}parallel bootstrap (reserved){p_end}
{synoptline}

{p 4 6 2}
{it:depvar} is the outcome variable. Control variables can be specified 
directly after {it:depvar} in the varlist. {it:covariates} supports factor-variable 
notation; see {help fvvarlist}. For example, {cmd:diddesign_check y i.region, ...} 
automatically expands categorical variable {cmd:region} into dummy variables.
{p_end}

{marker description}{...}
{title:Description}

{pstd}
{cmd:diddesign_check} performs diagnostic tests for the parallel trends assumption
in difference-in-differences designs. The command implements placebo tests that
estimate treatment effects in pre-treatment periods where no effect should exist
if parallel trends holds.

{pstd}
The command computes:

{phang2}
1. {bf:Placebo estimates}: DID estimates using pre-treatment periods only
{p_end}

{phang2}
2. {bf:Standardized estimates}: Placebo estimates standardized by the outcome 
standard deviation for comparability
{p_end}

{phang2}
3. {bf:Bootstrap standard errors}: Standard errors via block bootstrap
{p_end}

{phang2}
4. {bf:Equivalence confidence intervals (EqCI)}: 95% CIs for equivalence testing
{p_end}

{pstd}
{bf:Interpreting results:}

{pstd}
If the parallel trends assumption holds, placebo estimates should be close to zero.
The equivalence CI provides a way to assess how small the pre-trends are:

{phang2}
- Narrower EqCI bounds provide stronger evidence for parallel trends
{p_end}

{phang2}
- Researchers should interpret the width of the equivalence CI in the context
of their specific application, combining it with domain knowledge about what
constitutes a substantively meaningful deviation for the outcome of interest
{p_end}

{pstd}
For {bf:standard DID designs} ({cmd:design(did)}), placebo tests are computed 
directly from pre-treatment periods.

{pstd}
For {bf:staggered adoption designs} ({cmd:design(sa)}), placebo tests account 
for multiple treatment timing groups and compute time-weighted averages.

{marker options}{...}
{title:Options}

{dlgtab:Required}

{phang}
{opth treatment(varname)} specifies the binary treatment indicator variable.
Units with treatment equal to 1 are in the treatment group; units with 
treatment equal to 0 are in the control group.

{phang}
{opth time(varname)} specifies the time period variable. String variables are 
automatically encoded to numeric using {cmd:egen group()}.

{dlgtab:Data Type}

{pstd}
One of {opt id()} or {opt post()} must be specified to indicate the data type:

{phang}
{opth id(varname)} specifies the variable identifying units (individuals, 
firms, states, etc.). Required for {bf:panel data}. When {opt id()} is specified,
the command assumes panel data format unless {opt rcs} is also specified.
String variables are automatically encoded to numeric using {cmd:egen group()}.

{phang}
{opth post(varname)} specifies the binary post-treatment indicator variable (0/1).
Required for {bf:repeated cross-section (RCS) data}. In RCS data, different 
individuals are sampled at each time period, so there is no unit tracking.
The {opt treatment()} variable serves as the group indicator (treatment group 
vs control group), while {opt post()} indicates pre- vs post-treatment periods.

{dlgtab:Design}

{phang}
{opt design(did|sa)} specifies the design type. {cmd:did} (default) is for 
standard DID with a single treatment time. {cmd:sa} is for staggered adoption 
designs with multiple treatment times across units.
{bf:Note:} SA design only supports panel data; specifying {opt design(sa)} with 
{opt rcs} will result in an error.

{phang}
{opt panel} indicates panel data format. This is automatically assumed when 
{opt id()} is specified without {opt rcs}. For SA design, panel is required.

{phang}
{opt rcs} indicates repeated cross-section data format. When specified, 
{opt post()} is required to identify post-treatment periods. The {opt panel} 
and {opt rcs} options are mutually exclusive.

{dlgtab:Test Options}

{phang}
{opt lag(numlist)} specifies the lag periods for placebo tests. Default is 
{cmd:lag(1)}, testing parallel trends between periods t-1 and t-2. Multiple 
lags can be specified (e.g., {cmd:lag(1 2 3)}) to test multiple pre-treatment 
periods. Lag values must be positive integers less than the maximum available 
pre-treatment periods. Lag values that exceed this limit are automatically 
filtered out and a warning is displayed.

{phang}
{opt nboot(#)} specifies the number of bootstrap iterations for standard error 
estimation. Default is {cmd:nboot(30)}.

{phang}
{opt thres(#)} specifies the minimum number of treated units required per 
period for SA design. Default is {cmd:thres(2)}.

{dlgtab:Inference}

{phang}
{opth cluster(varname)} specifies the variable for cluster-robust bootstrap.
If not specified for panel data, clustering defaults to the unit level.

{phang}
{opt seed(#)} sets the random number seed for reproducibility.

{dlgtab:Display}

{phang}
{opt quiet} suppresses the bootstrap progress display. This is useful when 
running the command in batch mode, in loops, or when embedding 
{cmd:diddesign_check} in programs. Note that progress display is only shown 
for SA design ({cmd:design(sa)}); standard DID does not display progress.

{phang}
{opt parallel} enables parallel computation for bootstrap iterations.
This option is reserved for future implementation.

{marker results}{...}
{title:Stored results}

{pstd}
{cmd:diddesign_check} stores the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(n_lags)}}number of valid lag periods{p_end}
{synopt:{cmd:e(n_boot)}}number of bootstrap iterations{p_end}
{synopt:{cmd:e(n_boot_valid)}}number of valid bootstrap iterations{p_end}
{synopt:{cmd:e(n_clusters)}}number of clusters{p_end}
{synopt:{cmd:e(is_panel)}}1 if panel data, 0 if RCS data{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:diddesign_check}{p_end}
{synopt:{cmd:e(cmdline)}}command as typed{p_end}
{synopt:{cmd:e(design)}}design type ({cmd:did} or {cmd:sa}){p_end}
{synopt:{cmd:e(datatype)}}data type ({cmd:panel} or {cmd:rcs}){p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(treatment)}}name of treatment variable{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(covars)}}names of covariate variables{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(placebo)}}placebo test results matrix with columns:{p_end}
{synopt:}{it:lag}: lag value{p_end}
{synopt:}{it:estimate}: standardized placebo estimate{p_end}
{synopt:}{it:std_error}: bootstrap standard error (standardized){p_end}
{synopt:}{it:estimate_orig}: original (unstandardized) estimate{p_end}
{synopt:}{it:std_error_orig}: bootstrap standard error (original){p_end}
{synopt:}{it:EqCI95_LB}: 95% equivalence CI lower bound{p_end}
{synopt:}{it:EqCI95_UB}: 95% equivalence CI upper bound{p_end}
{synopt:{cmd:e(trends)}}trends data matrix for plotting{p_end}
{synopt:{cmd:e(Gmat)}}treatment timing matrix (SA design only){p_end}

{marker methods}{...}
{title:Methods and formulas}

{pstd}
{bf:Placebo estimation}

{pstd}
For lag k, the placebo estimate is computed as:

{p 8 8 2}
tau_placebo(k) = DID(t-k, t-k-1)

{pstd}
where DID(t1, t2) is the difference-in-differences estimator comparing periods 
t1 and t2. Under parallel trends, tau_placebo(k) should equal zero.

{pstd}
{bf:Standardization}

{pstd}
The command always computes both standardized and unstandardized placebo 
estimates.

{pstd}
Estimates are standardized by the control group outcome standard deviation 
at the baseline period (pre-treatment control observations):

{p 8 8 2}
tau_std = tau / sd(Y_control_baseline)

{pstd}
where {it:Y_control_baseline} is the outcome for control units (Gi=0) at the 
baseline time period (It=0). This standardization allows comparison across 
different lag periods and outcome variables, providing effect sizes in 
standard deviation units (similar to Cohen's d).

{pstd}
{bf:Equivalence confidence intervals}

{pstd}
The equivalence CI is computed using the TOST (Two One-Sided Tests) method:

{p 8 8 2}
90% CI = estimate +/- 1.645 * SE

{p 8 8 2}
EqCI95 = (-nu, nu) where nu = max(|90% CI upper|, |90% CI lower|)

{pstd}
Narrower intervals provide stronger evidence that the placebo estimate is 
close to zero, supporting the parallel trends assumption. Researchers should 
interpret the width in the context of their specific application.

{pstd}
{bf:Bootstrap variance}

{pstd}
Standard errors are computed via block bootstrap at the cluster level:

{p 8 12 2}
1. Resample clusters with replacement
{p_end}
{p 8 12 2}
2. Recompute placebo estimates on bootstrap sample
{p_end}
{p 8 12 2}
3. Repeat n_boot times
{p_end}
{p 8 12 2}
4. SE = sd(bootstrap estimates)
{p_end}

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
