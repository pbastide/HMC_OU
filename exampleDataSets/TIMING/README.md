# Timing comparisons between MH and HMC

* `xml_templates`: `xml` files to fit the various models described in
  the manuscript with BEAST, on both datasets.

* `write_xml.R`: helper `R` script to generate 10 identical copies of the template files. Each file must then be launched with beast, eg with different seeds, to study the variability of the timing of the various approaches.

* `summary.R`: helper `R` script that read the results, and produce a summary of the timings of the various approaches.
