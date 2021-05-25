## Model fitting satellite-based excess mortality data in Mogadishu, Somalia

**repo under construction/cleanup**

This repo contains scripts to run model fitting and to extract and visualise the results of fits. The scripts are accompanying [this preprint](https://www.medrxiv.org/) (to be uploaded).
To run simulations you need to download the age-structured compartmental COVID19 model COVIDM; available in several repos, for instance [here](https://github.com/nicholasdavies/newcovid/tree/master/covidm_for_fitting).

Annotation for the scripts:

-  *somalia_data_model.R*: read in satellite data on burials and NPI stringency. Set up COVIDM, run and visualised individual simulations and compare to data. Read in ACLED and visualise ACLED data on deaths due to political violence. 

- *somalia_data_model_fcns.R*: contains custom-made functions for visualisation, data processing etc.

- *run_multiple_mcmc.R*: run Markov Chain Monte Carlo fitting of COVIDM with excess death data, with variable inputs (IFR, seed size, fitting window)

- *extract_multiple_mcmc_scan_IFR_seedsize.R*: 

- *extract_multiple_mcmc_scan_IFR.R*: 