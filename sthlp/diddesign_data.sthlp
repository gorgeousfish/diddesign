{smcl}
{* *! version 1.0.0  20feb2026}{...}
{vieweralsosee "diddesign" "help diddesign"}{...}
{vieweralsosee "diddesign_check" "help diddesign_check"}{...}
{vieweralsosee "diddesign_plot" "help diddesign_plot"}{...}
{vieweralsosee "diddesign_intro" "help diddesign_intro"}{...}
{vieweralsosee "diddesign_data" "help diddesign_data"}{...}
{viewerjumpto "Syntax" "diddesign_data##syntax"}{...}
{viewerjumpto "Description" "diddesign_data##description"}{...}
{viewerjumpto "Options" "diddesign_data##options"}{...}
{viewerjumpto "Examples" "diddesign_data##examples"}{...}
{viewerjumpto "Datasets" "diddesign_data##datasets"}{...}
{viewerjumpto "Also see" "diddesign_data##alsosee"}{...}

{title:Title}

{phang}
{bf:diddesign_data} {hline 2} Load example datasets for the diddesign package

{marker syntax}{...}
{title:Syntax}

{pstd}
Load a dataset

{p 8 17 2}
{cmd:diddesign_data}
{it:dataname}
[{cmd:,} {opt clear}]

{pstd}
List available datasets

{p 8 17 2}
{cmd:diddesign_data}
{cmd:,} {opt list}

{pstd}
Describe a dataset without loading

{p 8 17 2}
{cmd:diddesign_data}
{it:dataname}
{cmd:,} {opt describe}

{pstd}
where {it:dataname} is one of {bf:anzia2012}, {bf:malesky2014}, or {bf:paglayan2019}.

{marker description}{...}
{title:Description}

{pstd}
{cmd:diddesign_data} loads example datasets included with the {cmd:diddesign}
package. The command uses Stata's {cmd:findfile} to dynamically locate data
files along the {cmd:adopath}, so no hardcoded paths are needed.

{pstd}
Three example datasets are provided, covering standard DID with panel data,
DID with repeated cross-sectional data, and staggered adoption designs.

{marker options}{...}
{title:Options}

{phang}
{opt clear} replaces the data in memory. If data in memory has been changed
and {opt clear} is not specified, an error is issued.

{phang}
{opt list} displays a table of all available datasets with their names,
observation counts, variable counts, and descriptions.

{phang}
{opt describe} displays the variable metadata of the specified dataset
(using {cmd:describe using}) without loading the data into memory.

{marker datasets}{...}
{title:Available datasets}

{synoptset 18 tabbed}{...}
{synopthdr:Dataset}
{synoptline}
{synopt:{bf:anzia2012}}Panel data, Texas school districts (2003-2009).
Anzia (2012, AJPS). 6,377 obs, 11 vars.{p_end}
{synopt:{bf:malesky2014}}Repeated cross-section, Vietnam communes (2006, 2008, 2010).
Malesky, Nguyen & Tran (2014, APSR). 6,269 obs, 42 vars.{p_end}
{synopt:{bf:paglayan2019}}Staggered adoption, US states (1959-2000).
Paglayan (2019, AJPS). 2,058 obs, 9 vars.{p_end}
{synoptline}

{marker examples}{...}
{title:Examples}

{pstd}
Load a dataset:

{phang2}{cmd:. diddesign_data malesky2014, clear}{p_end}

{pstd}
List all available datasets:

{phang2}{cmd:. diddesign_data, list}{p_end}

{pstd}
Describe a dataset without loading:

{phang2}{cmd:. diddesign_data paglayan2019, describe}{p_end}

{pstd}
Typical workflow:

{phang2}{cmd:. diddesign_data malesky2014, clear}{p_end}
{phang2}{cmd:. diddesign_check pro4, treatment(treatment) time(year) post(post_treat) rcs cluster(id_district) lag(1) nboot(200)}{p_end}
{phang2}{cmd:. diddesign pro4, treatment(treatment) time(year) post(post_treat) rcs cluster(id_district) nboot(200)}{p_end}

{marker alsosee}{...}
{title:Also see}

{psee}
Online: {helpb diddesign}, {helpb diddesign_check}, {helpb diddesign_plot}, {helpb diddesign_intro}
{p_end}
