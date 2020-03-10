# Simulation and analyses according to the scenario described in the manuscript.

* `simulations_BEEHIVE_campain.R`: `R` script to generate the simulated data.
  This script calls `write_xml_hmc.R` and `write_xml_factor_hmc.R` to generate
  the relevant `xml` files for subsequent inference.

* `blanquartMSM_sim_template_hmc.xml`: template `xml` file for the inference.
  This `xml` file is edited by the previous `R` script.

* `simulations_BEEHIVE_campaign_heritability.R`: Heritability computations for
  the simulation results. To be run after the `BEAST` analyses, on the resulting
  `log` files.
