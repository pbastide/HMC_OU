# Repository to reproduce the analyses in Bastide et al. 2020.

## Content

The general structure is the following.
Please refer to the `README` files in each sub-directory for more information on the analyses.

* `exampleDataSets`: data and scripts to reproduce the analyses
	* `BEEHIVE`: Dataset from Blanquart et al. (2017) *PLoS Biology*.
		* `fixedTree`: Inference of the model with a fixed tree.
		* `simulations`: Data simulation using the HIV fixed tree.
	* `SCHNITZLER`: Dataset from Schnitzler et al. (2017) *Evolutionary Ecology Research*.
		* `totalEvidence`: Joint inference of the tree and the process.
		* `fixedTree`: Model selection on the process on a fixed tree.

* `R_Utility_Files`: helper R functions.

## Requirements

* `BEAST`: http://beast.community/installing. The latest (development) version
  from `BEAST` is needed, please follow instructions to install the package
  from the `master` branch.

* `R`: https://cran.r-project.org/index.html. Version 3.6.

* `R` package `PhylogeneticEM`: https://CRAN.R-project.org/package=PhylogeneticEM. 
  Version 1.4 (for simulations and heritability computations).
