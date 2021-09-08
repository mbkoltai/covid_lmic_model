## Fitting satellite-based excess mortality data from Mogadishu (Somalia) Feb-Sept/2020

### Contact

Mihaly Koltai  
May 2021  
[CMMID, LSHTM](https://www.lshtm.ac.uk/aboutus/people/koltai.mihaly)

### Summary

This project is about fitting a SARS-CoV-2 epidemiological model to excess deaths data from Mogadishu (Somalia) during January-September 2020 to estimate the probable date of introduction and the basic reproduction number (R_0) of the pathogen in the region. The repository contains the accompanying code to [this preprint](https://www.medrxiv.org/node/368103.external-links.html). 
The scripts below can be used to reproduce figures in the publication.  
To run simulations you need to download the age-structured compartmental COVID19 model COVIDM (developed by [Nicholas Davies](https://github.com/nicholasdavies)), available eg. [here](https://github.com/nicholasdavies/newcovid/tree/master/covidm_for_fitting). Once downloaded, update the path to your COVIDM directory in the scripts (only needed once, to compile/launch CovidM).  
Input data files needed to reproduce simulations are in the folder [repo_data](https://github.com/mbkoltai/covid_lmic_model/tree/master/repo_data). The paths in the scripts below mostly point to the input files in *repo_data*, but in a few cases (eg. fitting output files that were not uploaded) to a local folder on the author's computer, in these cases change the path accordingly.

### Annotation for scripts

-  *somalia_data_model.R*: read in satellite data on burials and NPI stringency. Set up COVIDM, run and visualise individual simulations and compare to data. Read in ACLED and visualise ACLED data on deaths due to political violence. 

- *somalia_data_model_fcns.R*: contains custom-made functions for visualisation, data processing etc.

- *run_multiple_mcmc.R*: run Markov Chain Monte Carlo fitting of COVIDM with excess death data, with variable inputs (IFR, seed size, fitting window)

- *extract_multiple_mcmc_scan_IFR_seedsize.R*: extract results of fits where IFR and seed size were varied as input parameters. Used for full fits for period Jan/2020-Sep/2020 (or Febr/2020-Aug/2020)

- *extract_multiple_mcmc_scan_IFR.R*: extract results of fits where only IFR was varied as an input parameter. This was used for the fits with the truncated fitting window Febr/2020-Apr/2020.
