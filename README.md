# Source to reproduce results from Sun and Stephens (2018)

This repository contains source code to reproduce all results and plots of Sun and Stephens (2018), which develops `cashr` methods to solve the Empirical Bayes Normal Means (EBNM) with the correlated noise problem. The methods can be found in the [`cashr`](https://github.com/LSun/cashr) package.

If you find a bug, please create an [issue](https://github.com/LSun/cashr_paper/issues).

Citing this work
----------------

If you find any of the source code in this repository useful for your work, please cite our paper:

> Sun, Lei, and Matthew Stephens. 2018. "Empirical Bayes Normal Means with Correlated Noise."

License
-------

Copyright (c) 2018, Lei Sun, Matthew Stephens

All source code and software in this repository are made available under the terms of the [GNU General Public License](http://www.gnu.org/licenses/gpl.html). See the [LICENSE](LICENSE) file for the full text of the license.

Instructions
------------

To reproduce the results of Sun and Stephens (2018), you need to install the appropriate R packages.

### Install software and R packages

1.  Install [`R`](https://cran.r-project.org).
2.  Install [`RStudio`](https://www.rstudio.com/).
2.  Install the required `R` packages by running the following commands in the `R` interactive environment. (Note the order of these commands is important---the Bioconductor packages should be installed before the CRAN packages.)
``` r
source("https://bioconductor.org/biocLite.R")
biocLite(c("qvalue", "limma", "edgeR"), suppressUpdates = TRUE)
install.packages(c("ashr", "locfdr", "deconvolveR", "EbayesThresh", "ggplot2", "latex2exp", "devtools"))
```
3. Install the `Rmosek` package according to online instructions such as
- https://docs.mosek.com/8.1/rmosek/install-interface.html
- https://gist.github.com/mikelove/67ea44d5be5a053e599257fe357483dc
- https://rdrr.io/cran/ashr/f/inst/rmosek-mac.md
Once `Rmosek` is intalled, install additional `R` packages by running the following commands in the `R` interactive environment.
``` r
install.packages("REBayes")
devtools::install_github("LSun/cashr")
```

Data
----

The RNA-seq data set used for realistic simulations in this paper was created from the real human liver tissue RNA-seq gene expression read counts. In particular, the sample and gene identifiers have been removed from the data.

The RNA-seq gene expression data were originally generated by the GTEx Project, which was supported by the Common Fund of the Office of the Director of the National Institutes of Health, and by NCI, NHGRI, NHLBI, NIDA, NIMH, and NINDS. The data used for the analyses described in this paper were obtained from the GTEx Portal at [the GTEx Portal](https://www.gtexportal.org).

The leukemia data set used for real data illustrations in this paper is maintained by researchers ([Bradley Efron][] and [Trevor Hastie][]) at Stanford. The data set is available at http://statweb.stanford.edu/~ckirby/brad/LSI/datasets-and-programs/datasets.html.

The mouse data set used for real data illustrations in this paper was generated by Scott Adrian Smemo (RIP) and [Marcelo Nobrega][] at the University of Chicago.

All the data needed for reproducing the results and plots are stored in the `data/` directory.

-------------

This project is created with [workflowr][]. Special thanks to [David Gerard][], [Peter Carbonetto][], [John Blischak][].

[Bradley Efron]: http://statweb.stanford.edu/~ckirby/brad/
[Trevor Hastie]: https://web.stanford.edu/~hastie/
[Marcelo Nobrega]: http://nobregalab.uchicago.edu/
[David Gerard]: https://dcgerard.github.io/
[Peter Carbonetto]: https://pcarbo.github.io/
[John Blischak]: https://jdblischak.com/
[workflowr]: https://github.com/jdblischak/workflowr
