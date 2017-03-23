optimus
===========
[![Build Status](https://travis-ci.org/mitchest/optimus.svg?branch=master)](https://travis-ci.org/mitchest/optimus) [![codecov.io](https://codecov.io/github/mitchest/optimus/coverage.svg?branch=master)](https://codecov.io/github/mitchest/optimus?branch=master)

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

The best way to install optimus is to use Hadley Wickham's 
(excellent) **devtools** package - it's easy as to install optimus 
within the R environment, directly from github. Simply install 
**devtools** from CRAN with

    install.packages("devtools")

then call

    library(devtools)
	devtools::install_github("mitchest/optimus")

Once I stop frequently making changes I will put optimus on CRAN.

I have not yet figured out how to host binary packages for optimus.
If you have the correct devtools and compilers on your system,
then you can compile the package yourself from source (on github). 

### Development

This package has just come out of 'official' development, but 
functionality will continue to be refined and added. This is also 
the first package I have written for public access, so there will
inevitably be bugs and issues. If you find them, please let me know
about them - either directly on github, or the contact details below. 

## How to use optimus?
The vignette probably won't compile when installing from github,
so check out [the tutorial here](https://rawgit.com/mitchest/optimus/master/optimus-workflow.html).  
	
### Contact

* Mitchell Lyons
* mitchell.lyons@gmail.com / mitchell.lyons@unsw.edu.au
	
### References

Lyons et al. 2016. Model-based assessment of ecological community classifications. Journal of Vegetation Science: 27 (4) 704--715. DOI: http://dx.doi.org/10.1111/jvs.12400
