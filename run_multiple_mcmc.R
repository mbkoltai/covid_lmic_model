# multiple mcmc fits
fitting_params <- c("R0_fit","introd_date","ifr_logit_intercept") # "seed_size","compliance"
######
# function to input fitting params
pf <- function(parameters, x){
  x=as.list(x); n_groups=length(parameters$pop[[1]]$size)
# R0
target_R0=2 # cm_calc_R0(params,1) # params$pop[[1]]$u
parameters$pop[[1]]$u=c(rep(0.0145,4),rep(0.0305,12))*(x$R0_fit/target_R0)
# seed size and introd date
parameters$pop[[1]]$seed_times=rep(x$introd_date:x$introd_date,each=n_seedsize) # x$seed_size
# IFR
agegroupmeans=c(2.5+(0:14)*5,80.255); slope_val<-0.0999 # check by: as.numeric(linregr$coefficients["agegroup_mean"])
parameters$processes<-list(cm_multinom_process("Ip",
                  outcomes=data.table(death=inv.logit(x$ifr_logit_intercept+slope_val*agegroupmeans)/params$pop[[1]]$y),
                                          delays=data.table(death=data.table(death=cm_delay_gamma(22,22,60,1/4)$p)),report="o"))
# compliance (scaling of NPI)
t_npi=list(first=c("2020-03-19","2020-06-30"),second=c("2020-07-01","2020-08-29"),
           third=c("2020-08-30","2020-10-08"),fourth=c("2020-10-09","2020-11-01")); npi_vals=c(0.455,0.736,0.593,0.624)
for (k in 1:length(npi_vals)) { if (k==1) {iv=cm_iv_build(parameters)}
  cm_iv_contact(iv, t_npi[[k]][1], t_npi[[k]][2], 1 - (1-as.numeric(rep(npi_vals[k],4)) )*n_compliance) # x$compliance
  if (k==length(npi_vals)) {parameters$pop[[1]]$schedule=NULL; parameters=cm_iv_apply(parameters,iv)} }
return (parameters) 
}
#####
# function to calculate likelihood
likelihood = function(parameters, dynamics, data, x){
  inc = data; inc[, t := as.numeric(date - ymd(parameters$date0))]
  # simulations output scaled!!!
  eval = merge(dynamics[compartment == "death_o", .(model_case = sum(value)/parameters$scaling), by = t], inc, by = "t");
  ll = sum(dpois(eval$new_deaths, lambda = pmax(0.1, eval$model_case), log = T)); return (ll) }
####
# priors
priors=list(R0_fit="N 3 1 T 0.5 5.5", introd_date="N 120 20",ifr_logit_intercept="N -10.8 -2") # ,compliance="U 0 1" #  T 10 120
CDR_vals=c(1e4*baseline_daily_burials/mogadishu_popul,0.1,0.2,0.4)
fitting_date_window=as.Date(c("2020-01-15","2020-10-01")) # "2020-10-01"
# n_compliance,n_seedsize
parscan_dirname=paste0("simul_output/somalia/scan_seedsize_compliance_",
                       gsub(" ","_",paste0(paste0(names(priors)[2],"_",unlist(priors)[2],""),collapse="_")))
scan_seed_vals=round(10^((2:5)/2)); scan_compliance_vals=c(0,0.25,0.5); fits_death_scaling=list(); fits_compliance=list()
for (n_seedsize in scan_seed_vals) {
  for (n_compliance in scan_compliance_vals) {
    for (k in 1:length(CDR_vals)) { # introd_date ifr_logit_intercept
      scale_factor=(mogadishu_popul*CDR_vals[k]/1e4)/baseline_daily_burials; params$scaling=scale_factor
  fitting_incidence <- data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
    mutate(daily_baseline_subtr=round(daily_baseline_subtr))) %>%
    filter(date>=fitting_date_window[1] & date<=fitting_date_window[2]) %>% rename(new_deaths=daily_baseline_subtr)
  ### fitting -------------------
  print(paste0("#### fitting to CDR=",round(CDR_vals[k],3),"###, compliance=",n_compliance," ### seedsize=",n_seedsize," ########"))
  fits_death_scaling[[k]]=cm_fit(base_parameters=params, priors=priors, parameters_func=pf, likelihood_func=likelihood,
                                   data=fitting_incidence, mcmc_burn_in=5e2, mcmc_samples=2e3, mcmc_init_opt=F, opt_maxeval=25 )
if (k==length(CDR_vals)) {saveRDS(fits_death_scaling, paste0(parscan_dirname,"/fits_death_",paste0(fitting_date_window,collapse = "_"),
                                         "_compliance",n_compliance,"_seedsize",n_seedsize,".rds") )}
  }
  fits_compliance[[which(scan_compliance_vals %in% n_compliance)]]=fits_death_scaling
  } # fits_all[[which(scan_seed_vals %in% n_seedsize)]]=fits_compliance
}
# SAVE
# saveRDS(fits_all,file=paste0("simul_output/somalia/fits_death_",paste0(fitting_date_window,collapse = "_"),".rds"))
