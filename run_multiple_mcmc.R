# LOAD DATA and settings
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
lapply(c("tidyverse","deSolve","qs","gtools","rstudioapi","wpp2019","countrycode","coronavirus","wesanderson","dttr2","RcppRoll",
         "scales","wpp2019","GGally","corrr","ungeviz"), library,character.only=TRUE)
# functions and plotting theme
source("somalia_data_model_fcns.R")
# burial data
burial_data=read_csv("repo_data/out_bdr_daily_estimates.csv")
baseline_daily_burials=mean(subset(burial_data,date>="2019-07-01" & date<="2019-11-01")$new_graves_best_ipol)
# subset for relevant period and columns
out_bdr_daily_estimates=burial_data[!rowSums(is.na(burial_data))==(ncol(burial_data)-1),
                                    !colSums(is.na(burial_data))==nrow(burial_data)] %>% filter(date>"2019-11-01") %>%
  mutate(daily_baseline_subtr=ifelse(new_graves_best_ipol-baseline_daily_burials>0,
                                     new_graves_best_ipol-baseline_daily_burials,0),
         rollmeanweek=roll_mean(daily_baseline_subtr,7,align="center", fill=NA), # rolling mean BASELINE subtracted
         rollsumweek=roll_sum(daily_baseline_subtr,7,align="left",fill=NA),
         rollmeanweek_no_subtr=roll_mean(new_graves_best_ipol,7,align="center", fill=NA),
         rollsumweek_no_subtr=roll_sum(new_graves_best_ipol,7,align="left",fill=NA))
# Oxford Stringency Index
OxCGRT_url="https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"
OxCGRT_input=fcn_get_OxCGRT(OxCGRT_url,"Somalia") %>% mutate(OxCGRT_scaled_smoothed=roll_mean(OxCGRT_scaled,30,align="center",fill=NA))
# separate into 4 phases
NPI_phases=list(first=c("2020-03-19","2020-06-30"),second=c("2020-07-01","2020-08-29"),
                third=c("2020-08-30","2020-10-08"),fourth=c("2020-10-09","2020-11-01"))
NPIvals=sapply(NPI_phases,function(x) mean(OxCGRT_input$OxCGRT_scaled[OxCGRT_input$date>as.Date(x)[1]&OxCGRT_input$date<as.Date(x)[2]]))
npi_df=left_join(data.frame(t(data.frame(on_off=c("on","off"),NPI_phases))) %>% rownames_to_column(var="name") %>% 
                   filter(!name=="on_off") %>% rename(on=X1,off=X2) %>% mutate(on=as.Date(on),off=as.Date(off)),
                 data.frame(NPIvals) %>% rownames_to_column(var="name"),by="name") %>% mutate(name=factor(name,levels=unique(name))) %>%
  rename(contact_level=NPIvals) %>% mutate(contact_reduction=1-contact_level)
# IFR estimates (from: https://www.imperial.ac.uk/mrc-global-infectious-disease-analysis/covid-19/report-34-IFR/)
data(pop)
somalia_agegroups_IFR=fcn_merge_ifr_above_age(left_join(fcn_load_age_str("Somalia",n_year="2015",90),
          fcn_load_ifr("repo_data/IFR_by_age_imperial.csv"),by=c("agegroup","agegroup_min")),75) %>% 
  mutate(ifr_mean=ifelse(ifr_mean==0,min(ifr_mean[ifr_mean>0]),ifr_mean),log_ifr=log(ifr_mean),logit_ifr=log(ifr_mean/(1-ifr_mean)))
# IFR estimates from Sandmann 2021 cmmid paper
IFR_estimates_Sandmann2021 <- read_csv("repo_data/IFR_estimates_Sandmann2021.csv")
if (any(IFR_estimates_Sandmann2021$value_percent>1)) {n_cols<-2:ncol(IFR_estimates_Sandmann2021)
IFR_estimates_Sandmann2021[,n_cols]<-IFR_estimates_Sandmann2021[,n_cols]/1e2
IFR_estimates_Sandmann2021 <- left_join(IFR_estimates_Sandmann2021 %>% rename(agegroup=Age,ifr_mean=value_percent), 
    somalia_agegroups_IFR %>% select(!c(ifr_mean,log_ifr,logit_ifr)),by="agegroup") %>% mutate(logit_ifr=log(ifr_mean/(1-ifr_mean))) }
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# COVIDM
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
cm_force_rebuild=F; cm_build_verbose=T; cm_version=1; source(file.path(cm_path,"R","covidm.R"))
countryval="Somalia"; params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval),
                                                 date_start="2019-11-01",date_end="2020-10-01")
# set population: Somalia --> Mogadishu
params$pop[[1]]$name=countryval; mogadishu_popul=2.2e6
params$pop[[1]]$size=somalia_agegroups_IFR$agegroupsize*(mogadishu_popul/sum(somalia_agegroups_IFR$agegroupsize))
params$pop[[1]]$dist_seed_ages=cm_age_coefficients(30,70,5*(0:length(params$pop[[1]]$size)))
# set clinical fraction values (from Davies 2020 -> "repo_data/suscept_clinfract_posteriors_davies2010.csv")
# set approximated clin fract values
params$pop[[1]]$y=fun_lin_approx_agedep_par(agegroups=somalia_agegroups_IFR,min_val=0.25,max_val=0.7,rep_min=6,rep_max=2)
###
# multiple mcmc fits
fitting_params <- c("R0_fit","introd_date","npi_scale") # "seed_size","compliance"
### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ###
# functions for MCMC fitting
# pf handles fitting parameters
pf <- function(parameters, x){x=as.list(x); n_groups=length(parameters$pop[[1]]$size)
# R0
target_R0=2; parameters$pop[[1]]$u=c(rep(0.0145,4),rep(0.0305,12))*(x$R0_fit/target_R0) #  # cm_calc_R0(params,1)
# seed size and introd date: continuous seeding or one-time event?
seeding_t_window=round(x$introd_date):round(x$introd_date+round(n_seedsize/10))
parameters$pop[[1]]$seed_times=rep(seeding_t_window,each=10) # x new infections/day for n days
# parameters$pop[[1]]$seed_times=rep(round(x$introd_date):round(x$introd_date),each=n_seedsize) # x$seed_size
# IFR
agegroupmeans=c(2.5+(0:14)*5,80.255) # slope_val<-0.0999 # as.numeric(linregr$coefficients["agegroup_mean"])
# estimates used in Sandmann 2021 (and other CMMID publs)
ifr_logit <- c(-10.414283,-11.512915,-11.512915,-10.414283,-9.721106,-8.947846,-8.334632,-7.823646,-7.194687,
               -6.715924,-6.178135,-5.732038,-5.385862,-4.522041,-4.073073,-2.026972)
# x$ifr_logit_intercept+slope_val*agegroupmeans
parameters$processes<-list(cm_multinom_process("Ip",
                                               outcomes=data.table(death=inv.logit(ifr_logit + ifr_logit_increm)/parameters$pop[[1]]$y),
                                               delays=data.table(death=data.table(death=cm_delay_gamma(22,22,60,1/4)$p)),report="o"))
# compliance
t_npi=list(first=c("2020-03-19","2020-06-30"),second=c("2020-07-01","2020-08-29"),
           third=c("2020-08-30","2020-10-08"),fourth=c("2020-10-09","2020-11-01")); npi_vals=c(0.455,0.736,0.593,0.624)
for (k in 1:length(npi_vals)) { if (k==1) {iv=cm_iv_build(parameters)}
  cm_iv_contact(iv, t_npi[[k]][1], t_npi[[k]][2], 1 - (1-as.numeric(rep(npi_vals[k],4)) )*x$npi_scale)
  if (k==length(npi_vals)) {parameters$pop[[1]]$schedule=NULL; parameters=cm_iv_apply(parameters,iv)} }
return (parameters) 
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# calculating posterior likelihood by a Poisson
likelihood = function(parameters, dynamics, data, x){
  inc = data; inc[, t := as.numeric(date - ymd(parameters$date0))]
  # simulations output scaled!!!
  eval = merge(dynamics[compartment == "death_o", .(model_case = sum(value)/parameters$scaling), by = t], inc, by = "t");
  ll = sum(dpois(eval$new_deaths, lambda = pmax(0.1, eval$model_case), log = T)); return (ll) }
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# priors
priors=list(R0_fit="N 3 1 T 1 5", introd_date="N 182 20",npi_scale="U 0 1")
### ### ### ### ### ### ### ###
# PERIOD we are fitting
fitting_date_window=as.Date(c("2020-02-23","2020-08-24")) # c("2020-01-15","2020-04-13") # c("2020-01-15","2020-10-01")
# select fitting data
fitting_incidence <- data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
                                  mutate(daily_baseline_subtr=round(daily_baseline_subtr))) %>%
  filter(date>=fitting_date_window[1] & date<=fitting_date_window[2]) %>% rename(new_deaths=daily_baseline_subtr)
# plot
ggplot(fitting_incidence,aes(x=date,y=new_deaths)) + geom_line() + geom_point() +
  theme_bw() + scale_x_date(date_breaks="2 days",expand=expansion(0,0)) + standard_theme +
  theme(axis.text.x=element_text(vjust=0.5),legend.position="top") # + geom_vline(xintercept=as.Date("2020-04-13"),size=0.4)
# n_compliance,n_seedsize
### ### ### ### ### ### ### ###
parscan_dirname=paste0("simul_output/somalia/scan_seedsize_ifr_",
                       gsub(" ","_",paste0(paste0(gsub("_","",names(priors)[2]),"_",unlist(priors)[2],""),collapse="_")),
                       "_fitperiod_",paste0(gsub("-","",fitting_date_window),collapse = "_"))
if (!dir.exists(parscan_dirname)) {dir.create(parscan_dirname); print("created dir")} else {print("dir exists")}
CDR_vals=c(1e4*baseline_daily_burials/mogadishu_popul,0.1,0.2,0.4)[1] 
# start date for simulations
params$date0="2019-09-01"
# select range of seed sizes
scan_seed_vals=c(20,50,100,200,500)
# select range of IFR values
ifr_increm_vals=c(2.5) # c(0,1,2,2.5,3)
fits_death_scaling=list(); fits_compliance=list()
# fcn_somal_sir_singlesimul<-.GlobalEnv$fcn_somal_sir_singlesimul
# df_output <- foreach(k=1:k_lim,.combine=rbind,.packages=c("tidyr","deSolve","dplyr","RcppRoll")) %dopar% { }
for (n_seedsize in scan_seed_vals) {
  for (ifr_logit_increm in ifr_increm_vals) {
    for (k in 1:length(CDR_vals)) { # introd_date ifr_logit_intercept
      scale_factor=(mogadishu_popul*CDR_vals[k]/1e4)/baseline_daily_burials; params$scaling=scale_factor
      fitting_incidence <- data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
                                        mutate(daily_baseline_subtr=round(daily_baseline_subtr))) %>%
        filter(date>=fitting_date_window[1] & date<=fitting_date_window[2]) %>% rename(new_deaths=daily_baseline_subtr)
      ### fitting -------------------
      print(paste0("#### fitting to CDR=",round(CDR_vals[k],3),"##, ifr_logit increm=",ifr_logit_increm," ### seedsize=",n_seedsize," ##"))
      fits_death_scaling[[k]]=cm_fit(base_parameters=params, priors=priors, parameters_func=pf, likelihood_func=likelihood,
                          data=fitting_incidence, mcmc_burn_in=5e2, mcmc_samples=2e3, mcmc_init_opt=F, opt_maxeval=25 )
      if (k==length(CDR_vals)) {
        saveRDS(fits_death_scaling, paste0(parscan_dirname,"/fits_death_",paste0(fitting_date_window,collapse = "_"),
                                           "_ifr_increm",ifr_logit_increm,"_seedsize",n_seedsize,".rds") )}
    }
    # fits_compliance[[which(scan_compliance_vals %in% n_compliance)]]=fits_death_scaling
  }
  # fits_all[[which(scan_seed_vals %in% n_seedsize)]]=fits_compliance
}
# SAVE
# saveRDS(fits_all,file=paste0("simul_output/somalia/fits_death_",paste0(fitting_date_window,collapse = "_"),".rds"))
