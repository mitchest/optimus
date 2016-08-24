optimus
===========

## What is optimus??

An R package for assessment and diagnostics of competing
clustering solutions, using predictive models. The main intended
use is for comparing clustering/classification solutions of
ecological data (e.g. presence/absence, counts, ordinal scores) to
1) find an optimal partitioning solution, 2) identify
characteristic species and 3) refine a classification by merging
clusters that increase predictive performance. However, in a more
general sense, this package can do the above for any set of
clustering solutions for i observations of j variables. More
details on the background and theory behind using predictive models
for classification assessment, in an ecological context, can be
found in Lyons et al. (2016).

## Development

This is a development version, and this is also the first package
I have written. There are certainly bug and issues, and if you
find them, please let me know about them - either directly on
github, or the contact details below.

## Installation

I have not yet figured out how to host binary packages for optimus.
If you have the correct dev tools and compilers on your system,
then you can compile the package yourself from source (on github).

Honestly, the easiest way to install (until I get optimus on CRAN),
is to use Hadley Wickham's (excellent) **devtools** package. Then
it's easy as to install optimus within the R environment, directly
from github. Simply install **devtools** from CRAN with

    install.packages("devtools")

then call

    library(devtools)
	devtools::install_github("mitchest/optimus")

	
### Contact

* Mitchell Lyons
* mitchell.lyons@gmail.com / mitchell.lyons@unsw.edu.au
	
### References

Lyons et al. 2016. Model-based assessment of ecological community classifications. Journal of Vegetation Science: 27 (4) 704--715. DOI: http://dx.doi.org/10.1111/jvs.12400
