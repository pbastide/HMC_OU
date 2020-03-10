# Analyses of the dataset from Blanquart et al. (2017) PLoS Biology.

* `xml_files_2019-11-19`: `xml` files to fit the various models described in
  the manuscript with BEAST.

* `BEEHIVE_heritability.R`: helper `R` functions to compute the heritability
  (according to the definition of the paper). This script should be run after
  the `xml`, on the log files produced during the analyses.

* `blanquartMSM.tree`: The ML tree used, from Blanquart et al. (2017), trimmed
  to sub-type B and MSM patients.
