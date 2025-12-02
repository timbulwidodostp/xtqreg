{smcl}
{* *! version 1.0.0  30mar2015}{...}
{cmd:help xtqreg}
{hline}

{title:Title}

{p2colset 5 17 19 2}{...}
{p2col :{hi:xtqreg} {hline 2}}Linear Quantile Mixed Models{p_end}
{p2colreset}{...}

{title:Syntax}
{phang}

{p 8 13 2}
{cmd:xtqreg} {depvar} [{indepvars}] {ifin} 
	[{cmd:,} {it:{help xtqreg##xtreg_options:xtqreg_options}}]

{synoptset 25 tabbed}{...}
{marker xtqreg_options}{...}
{synopthdr :xtqreg_options}
{synoptline}
{syntab :Model}
{synopt :{cmdab:q:uantiles(}{it:#}[{it:#}[{it:# ...}]]{cmd:)}}estimate {it:#} quantiles; default is {cmd:quantiles(.5)}{p_end}
{synopt :{opt r:eps(#)}}perform {it:#} bootstrap replications; default is {cmd:reps(20)}{p_end}
{synopt :{opt seed(string)}}seed for the bootstrap{p_end}
{synopt :{opt m:ethod(string)}}specifies the optimization algorithm{p_end}
{synopt :{opt r:andom(varlist)}}specifies the random effects{p_end}
{synopt :{opth cov:ariance(xtqreg##vartype:vartype)}}specifies the variance-covariance of the random effects{p_end}
{synopt :{opt noc:onstant}}suppress constant term{p_end} 

{syntab :Paths}
{synopt :{opt pathr(R_pathname)}}specifies a path name for invoking the R command{p_end}
{synopt :{opt ver:sion(R_version)}}specifies the R version to use{p_end}
{synopt :{opt dir:ectory(pathname)}}defines the working directory{p_end}

{syntab :Reporting}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synoptline}

{synoptset 23}{...}
{marker vartype}{...}
{synopthdr :vartype}
{synoptline}
{synopt :{opt ind:ependent}}one variance parameter per random effect, 
all covariances zero; the default unless a factor variable is specified{p_end}
{synopt :{opt ex:changeable}}equal variances for random effects, 
and one common pairwise covariance{p_end}
{synopt :{opt id:entity}}equal variances for random effects, all 
covariances zero; the default for factor variables{p_end}
{synopt :{opt un:structured}}all variances and covariances distinctly 
estimated{p_end}
{synoptline}
{p2colreset}{...}

{p2colreset}{...}
{phang}See {manhelp sqreg_postestimation R:sqreg postestimation} for features
available after estimation.

{title:Description}

{pstd}
{cmd:xtqreg} is used to fit linear quantile mixed models based on the asymmetric Laplace distribution. 
The command {cmd:xtqreg} is a wrapper for the {browse "http://cran.r-project.org/web/packages/lqmm/":lqmm} package developed by Marco Geraci in 
{browse "http://cran.r-project.org/":R}. Therefore {browse "http://cran.r-project.org/":R} needs to be installed together with the package {browse "http://cran.r-project.org/web/packages/lqmm/":lqmm}.

{title:Options for xtqreg}

{dlgtab:Model}

{phang}{cmd:quantiles(}{it:#} [{it:#} [{it:#} {it:...}]]{cmd:)} specifies the quantiles as numbers between 0 and 1;
numbers larger than 1 are interpreted as percentages. The default value is 0.5, which corresponds to the median.

{phang}{opt reps(#)} specifies the number of bootstrap replications for estimating 
variance-covariance matrix and standard errors of the regression coefficients.

{phang}{opt seed(string)} specifies the seed for the bootstrap for estimating 
variance-covariance matrix and standard errors of the regression coefficients.

{phang}
{opt r:andom} specifies the random effects. The default is the intercept.

{phang}
{opt m:ethod} specifies the optimization algorithm.
The optimization algorithm is based on the gradient of the Laplace log-likelihood (gs), the default.
An alternative optimization algorithm is based on a Nelder-Mead algorithm (df).

{phang}
{opt noconstant}; see
{helpb estimation options##noconstant:[R] estimation options}.

{dlgtab:Paths}

{phang}
{opt pathr(R_pathname)} specifies the path name for invoking the R command.
The default path is "C:\Program Files\R\R-4.5.0\bin\x64\R.exe". The R version 
can be changed with option {cmd:version()}.
The R_pathname can also be set by defining a {help macro:global macro} {hi:Rterm_path}
 (See {help rsource:rsource}, {hi:{help rsource##rsource_technote:Technical note}}). 

{phang}
{opt ver:sion(R_version)} defines the working R version. The default is "4.5.0".

{phang}
{opt dir:ectory(pathname)} defines the working directory.
The default directory is the current directory. If the current directory is on the 
server the program will return an error. The alternative directory must
be on the local computer. 

{dlgtab:Reporting}

{phang}{opt level(#)}; see 
{helpb estimation options##level():[R] estimation options}.

{title:Examples}

{pstd}Set the path for Windows{p_end}
{phang2}{stata `"global Rterm_path "C:\Program Files\R\R-4.5.0\bin\x64\Rterm.exe""'}{p_end}

{pstd}Set the path for Mac{p_end}
{phang2}{stata `"global Rterm_path "/usr/bin/r""'}{p_end}
 
{phang2}{stata `"import excel "https://raw.githubusercontent.com/timbulwidodostp/xtqreg/main/xtqreg/xtqreg.xlsx", sheet("Sheet1") firstrow clear"'}{p_end}

{phang2}{stata "xtset id month"}{p_end}

{pstd}Random intercept{p_end}
{phang2}{stata "xtqreg y  prog month  inter"}{p_end}

{pstd}Random intercept and multiple quantiles{p_end}
{phang2}{stata "xtqreg y  prog month  inter, q(25 50 75)"}{p_end}

{pstd}Random intercept and random slope{p_end}
{phang2}{stata "xtqreg y  prog month  inter, random(month)"}{p_end}

{title:Author}

{phang2}Timbul Widodo{p_end}

{phang2}Olah Data Semarang{p_end}

{phang2}www.youtube.com/@amalsedekah{p_end}

{phang2}Whatsapp +6285227746673 (085227746673){p_end}

{title:References}

{phang2}Geraci M and Bottai M (2007). Quantile regression for longitudinal data using the asymmetric Laplace distribution. Biostatistics 8(1), 140-154.{p_end}

{phang2}Bottai M, Orsini N and Geraci M (2014). A Gradient Search Maximization Algorithm for Laplace Likelihood. Journal of Statistical Computation and Simulation.{p_end}

{phang2}Geraci M and Bottai M (2014). Linear Quantile Mixed Models. Statistics and Computing. May 2014, Volume 24, Issue 3, pp 461-479.{p_end}

{phang2}Geraci M (2014). Linear Quantile Mixed Models: The lqmm Package for Laplace Quantile Regression. Journal of Statistical Software. Vol. 57, Issue 13, May 2014.{p_end}

{title:Also see}

{psee}
Manual:  {bf:[R] qreg}

{psee}
Online:  {manhelp qreg_postestimation R:qreg postestimation};{break}
{manhelp bootstrap R}
{p_end}
