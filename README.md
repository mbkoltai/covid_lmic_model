## Model fitting satellite-based excess mortality data in Mogadishu, Somalia

### Contact

Mihaly Koltai  
May 2021  
[CMMID, LSHTM](https://www.lshtm.ac.uk/aboutus/people/koltai.mihaly)

### Summary

This project is about fitting a SARS-CoV-2 epidemiological model to excess deaths data from Mogadishu (Somalia) during January-October 2020 to estimate the probable date of introduction and the basic reproduction number (R_0) of the pathogen in the region. The repository contains the accompanying code to [this preprint](https://www.medrxiv.org/) (to be uploaded). 
The scripts below can be used to reproduce figures in the publication.  
To run simulations you need to download the age-structured compartmental COVID19 model COVIDM, available in several repos, for instance [here](https://github.com/nicholasdavies/newcovid/tree/master/covidm_for_fitting). Once COVIDM is downloaded, update the paths to your directory in the scripts below.  
Input data files needed to reproduce simulations can be found in the folder [repo_data](https://github.com/mbkoltai/covid_lmic_model/tree/master/repo_data). The paths in the scripts below point to the local folder structure, you will need to change paths to point to *repo_data* instead when using the scripts.

### Annotation for scripts

-  *somalia_data_model.R*: read in satellite data on burials and NPI stringency. Set up COVIDM, run and visualised individual simulations and compare to data. Read in ACLED and visualise ACLED data on deaths due to political violence. 

- *somalia_data_model_fcns.R*: contains custom-made functions for visualisation, data processing etc.

- *run_multiple_mcmc.R*: run Markov Chain Monte Carlo fitting of COVIDM with excess death data, with variable inputs (IFR, seed size, fitting window)

- *extract_multiple_mcmc_scan_IFR_seedsize.R*: extract results of fits where IFR and seed size were varied as input parameters. Used for full fits for period Jan/2020-Oct/2020 (or Febr/2020-Aug/2020)

- *extract_multiple_mcmc_scan_IFR.R*: extract results of fits where only IFR was varied as an input parameter. This was used for the fits with the truncated fitting window Febr/2020-Apr/2020.