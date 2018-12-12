# Source to reproduce results from Sun and Stephens (2018)

This repository contains source code to reproduce all results and plots of Sun and Stephens (2018), which develops `cashr` methods to solve the Empirical Bayes Normal Means (EBNM) with the correlated noise problem. The methods can be found in the [`cashr`](https://github.com/LSun/cashr) package.

If you find a bug, please create an [issue](https://github.com/LSun/cashr_paper/issues).

Citing this work
----------------

If you find any of the source code in this repository useful for your work, please cite our paper:

> Sun, Lei, and Matthew Stephens. 2018. "Empirical Bayes Normal Means with Correlated Noise."

License
-------

Copyright (c) 2018, Lei Sun.

All source code and software in this repository are made available under the terms of the [GNU General Public License](http://www.gnu.org/licenses/gpl.html). See the [LICENSE](LICENSE) file for the full text of the license.

Instructions
------------

To reproduce the results of Sun and Stephens (2018), you need to install the appropriate R packages.

### Install software and R packages

1.  Install [R](https://cran.r-project.org).
3.  Install the required R packages by running the following commands in the R interactive environment. (Note the order of these commands is important---the Bioconductor packages should be installed before the CRAN packages.)

``` r
source("https://bioconductor.org/biocLite.R")
biocLite(c("qvalue", "limma", "edgeR"), suppressUpdates = TRUE)
install.packages(c("ashr", "locfdr", "deconvolveR", "EbayesThresh", "ggplot2"))
devtools::install_github("LSun/cashr")
```

-------------

This project is created with [workflowr][]. Grateful acknowledgement to [David Gerard][], [Peter Carbonetto][], [John Blischak][].

[David Gerard]: https://dcgerard.github.io/
[Peter Carbonetto]: https://pcarbo.github.io/
[John Blischak]: https://jdblischak.com/
[workflowr]: https://github.com/jdblischak/workflowr
