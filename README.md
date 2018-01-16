optimus
===========
[![](http://cranlogs.r-pkg.org/badges/optimus)](https://cran.r-project.org/package=optimus)  
[![Build Status](https://travis-ci.org/mitchest/optimus.svg?branch=master)](https://travis-ci.org/mitchest/optimus)[![Code coverage](https://codecov.io/gh/mitchest/optimus/branch/master/graphs/badge.svg?branch=master)](https://codecov.io/gh/mitchest/optimus/)  

## What is optimus??

An R package for assessment and diagnostics of competing
clustering solutions, using predictive models. The main intended
use is for comparing clustering/classification solutions of
ecological data (e.g. presence/absence, counts, ordinal scores) to:

1) find an optimal partitioning solution  
2) identify characteristic species and  
3) refine a classification by merging clusters such that it
increases predictive performance.  

However, in a more general sense, this package can do the above for 
any set of clustering solutions for i observations of j variables. 
More details on the background and theory behind using predictive 
models for classification assessment, in an ecological context, 
can be found in Lyons et al. (2016).

## Installation

In R, simply use

    install.packages("optimus")

See the package page on CRAN for more details:  
https://cran.r-project.org/package=optimus  
  
### Development version
If you want to install the development version of optimus,
for example if I've added something new that you want to use,
but it's not yet up on CRAN, then you can also install directly
from github. It's very easy - simply use Hadley Wickham's 
(excellent) **devtools** package - install **devtools** from
 CRAN within R using

    install.packages("devtools")

then call

    library(devtools)
	devtools::install_github("mitchest/optimus")

### Bugs

There are some probably. If you find them, please let me know
about them - either directly on github, or the contact details below. 

## How to use optimus?
You can find the vignette on the CRAN home page, or you can access it
here too (might be new things here before CRAN occasionally).  
Check out [the tutorial here](https://rawgit.com/mitchest/optimus/master/optimus-workflow.html).  
	
### Contact

* Mitchell Lyons
* mitchell.lyons@gmail.com / mitchell.lyons@unsw.edu.au
	
### References

Lyons et al. 2016. Model-based assessment of ecological community classifications. Journal of Vegetation Science: 27 (4) 704--715. DOI: http://dx.doi.org/10.1111/jvs.12400
