#' ---
#' title: Somalia COVID19
#' author: Mihaly Koltai
#' date: February 2021
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#'      fig_width: 20
#'      fig_height: 10
#' ---

#' Load libraries, functions
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
lapply(c("tidyverse","deSolve","qs","gtools","rstudioapi","wpp2019","countrycode","coronavirus","wesanderson",
         "RcppRoll","scales","dttr2","wpp2019","foreach","parallel","doParallel"), library,character.only=TRUE)
# functions and plotting theme
source("covid_LIC_fcns.R")

### JHU global covid19 data ----------------
#' ## JHU global covid19 data
data("coronavirus")
# age structure
N_tot=fun_cntr_agestr("Somalia",i_year="2020",age_groups=data.frame(age_group=c(1:16),age_low=c(seq(0,75,5)),age_high=c(seq(4,74,5),100)))
# reported case and deaths data
covid_somal=coronavirus %>% filter(country %in% "Somalia") %>% mutate(rollingmean=roll_mean(cases,7,align="center",fill=NA)) %>%
  mutate(per_million=rollingmean/(sum(N_tot)/1e6),name=str_replace(type,"confirmed","confirmed cases")) %>% rename(value=cases)

#' ### Plot timecourse of cases/deaths
 ggplot(subset(covid_somal,!type %in% "recovered"),aes(x=date)) + 
  geom_line(aes(y=rollingmean)) + # geom_point(aes(y=cases),size=0.5,fill=NA,shape=1) + #
  scale_x_date(name="", limits=c(as.Date("2020-03-01"), NA_Date_),date_breaks="weeks") + facet_wrap(~name,scales="free",nrow=2) +
  geom_rect(aes(xmin=as.Date("2020-04-20"),xmax=as.Date("2020-04-24"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_rect(aes(xmin=as.Date("2020-05-10"),xmax=as.Date("2020-05-15"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_rect(aes(xmin=as.Date("2020-05-28"),xmax=as.Date("2020-06-01"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + theme_bw() + standard_theme + 
  labs(title="COVID19 in Somalia: confirmed cases and deaths (7-day rolling mean)",caption="source: JHU CCSE") + ylab("")
#' No lag in deaths compared to cases
# save
# tcourse_filename="simul_output/somalia_output/reported_cases_deaths.png" # _dots
# ggsave(tcourse_filename,width=30,height=18,units="cm") # _dots

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### SATELLITE DATA on cemeteries ----------------
#' satellite image data
out_bdr_daily_estimates=read_csv("data/somalia/Mogadishu_data/mogadishu_burial_analysis-main/out_bdr_daily_estimates.csv")
# baseline of daily burials 
baseline_daily_burials=mean(subset(out_bdr_daily_estimates,date>="2019-07-01" & date<="2019-11-01")$new_graves_best_ipol)
# subset for relevant period and columns
out_bdr_daily_estimates=out_bdr_daily_estimates[!rowSums(is.na(out_bdr_daily_estimates))==(ncol(out_bdr_daily_estimates)-1),
      !colSums(is.na(out_bdr_daily_estimates))==nrow(out_bdr_daily_estimates)] %>% filter(date>"2019-11-01") %>%
      mutate(daily_baseline_subtr=ifelse(new_graves_best_ipol-baseline_daily_burials>0,new_graves_best_ipol-baseline_daily_burials,0),
      rollmeanweek=roll_mean(daily_baseline_subtr,7,align="center", fill=NA), # rolling mean BASELINE subtracted
      rollsumweek=roll_sum(daily_baseline_subtr,7,align="left",fill=NA))

# plot number of burials (7-day mean and sum)
ggplot(subset(out_bdr_daily_estimates %>% pivot_longer(col=c(new_graves_best_ipol,rollmeanweek,rollsumweek)),grepl("roll",name)),
              aes(x=date,y=value)) + facet_wrap(~name,scales="free") + geom_line() + geom_point(size=0.3) + theme_bw() + standard_theme + 
  scale_x_date(date_breaks="2 weeks",expand=expansion(0.0)) + theme(axis.text.x=element_text(vjust=0.5))

# compare to reported deaths
weekly_deaths_reported=data.frame(subset(covid_somal,name %in% "death")[,c("date","value")],datasource="reported") %>% 
  mutate(datasource=as.character(datasource)) %>% rename(value_report_daily=value) %>%
  mutate(rollmeanweek=roll_mean(value_report_daily,7,align="center",fill=NA),
         rollsumweek=roll_sum(value_report_daily,7,fill=NA,align="right"))
# plot together
df_compare_report_satell=bind_rows(weekly_deaths_reported,
  out_bdr_daily_estimates %>% select(date,daily_baseline_subtr,rollmeanweek,rollsumweek) %>% mutate(datasource="satellite") ) %>% 
  pivot_longer(cols=!c(date,datasource)) %>% filter(!is.na(value)) %>% #  & name=="rollsumweek"
  mutate(type=ifelse(!grepl("sum|mean",name),"daily",name),week=format(date,"%Y/%W")) %>%
  filter(type=="daily") %>% group_by(week,datasource) %>% 
  summarise(date=min(date),datasource=unique(datasource),name=unique(name),value=sum(value)) %>% mutate(week=gsub("/","/w",week))
# plot
ggplot(subset(df_compare_report_satell,date>"2020-01-15" & date<"2020-10-07"),aes(x=week,y=value,group=datasource)) + 
  # geom_line(aes(color=datasource)) + geom_point(aes(color=datasource,shape=datasource),size=0.9) + 
  geom_bar(aes(fill=datasource),stat="identity",position=position_dodge(width=0.75),color="black",size=0.2) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) + # facet_wrap(~type,scales="free",nrow=3) + 
 # scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0),limits=as.Date(c("2020-01-10","2020-08-17"))) +
 scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*5) + ggtitle(expression("weekly burials - confirmed COVID19 deaths (weekly sum)"))
# SAVE
# ggsave("simul_output/somalia_output/satellite_burials_reported_deaths_weekly.png",units="cm",height=15,width=25)
# ggsave("simul_output/somalia_output/satellite_burials_reported_deaths_weekly_barplot.png",units="cm",height=15,width=25)

# Rt estimate: https://epiforecasts.io/covid/posts/national/somalia/

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Oxford Stringency Index
# truncate until a given timepoint
OxCGRT_url="https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"
OxCGRT_input = fcn_get_OxCGRT(OxCGRT_url,"Somalia") %>% mutate(OxCGRT_scaled_smoothed=roll_mean(OxCGRT_scaled,30,align="center",fill=NA))
# Somalia population
somal_popul_tot=sum(N_tot); mogadishu_popul=2.2e6

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# COVIDM
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
cm_force_rebuild=F; cm_build_verbose=T; cm_version=1; source(file.path(cm_path,"R","covidm.R"))
countryval="Somalia"; params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval),date_start="2019-11-01",date_end="2020-09-01")
params$pop[[1]]$name=countryval
# seeding
params$pop[[1]]$seed_times = rep(31:61,each=1) # 5 new infections each day for first 7 days
# infections start in individuals aged 20-50
params$pop[[1]]$dist_seed_ages=cm_age_coefficients(20,60,5*(0:length(params$pop[[1]]$size)))
### incident deaths -------
### set up age structure + IFR
somalia_agegroups_IFR=fcn_merge_ifr_above_age(left_join(fcn_load_age_str("Somalia",90),fcn_load_ifr("data/IFR_by_age_imperial.csv"),
                                              by=c("agegroup","agegroup_min")),75); somalia_agegroups_IFR$ifr_mean[1]=3e-6

### add death process to model ------
params$processes <- list(cm_multinom_process("Ip",outcomes=data.table(death=somalia_agegroups_IFR$ifr_mean), 
                                             delays=data.table(death=ponset2death), report="o"))
# RUN SIMUL
delta_clinfr=0.3
params$pop[[1]]$y=fun_lin_approx_agedep_par(agegroups=somalia_agegroups_IFR,min_val=0.1,max_val=0.75,rep_min=3,rep_max=5)
target_R0=1.8; scale_r0=target_R0/cm_calc_R0(params,1)
# change susceptibility to get R0
params$pop[[1]]$u=params$pop[[1]]$u*(target_R0/cm_calc_R0(params,1))
# set population to Somalia --> Mogadishu
N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups=data.frame(age_group=c(1:16),age_low=c(seq(0,75,5)),age_high=c(seq(4,74,5),100)))
params$pop[[1]]$size=N_tot*mogadishu_popul/sum(N_tot)
# NPIs
# first="2020-03-19","2020-05-28"
NPI_phases=list(first=c("2020-03-19","2020-05-28"),second=c("2020-05-29","2020-06-30"),third=c("2020-07-01","2020-09-24"))
NPIvals=sapply(NPI_phases,function(x) mean(OxCGRT_input$OxCGRT_scaled[OxCGRT_input$date>as.Date(x)[1]&OxCGRT_input$date<as.Date(x)[2]]))
k_compl=0.1
for (k in 1:length(NPIvals)) { # setup for version 1
  if (cm_version==1) {if (k==1) {params$pop[[1]]$schedule=NULL; iv=cm_iv_build(params)} # sets up data structure for interventions
    cm_iv_contact(iv, NPI_phases[[k]][1], NPI_phases[[k]][2], 1-(1-as.numeric(rep(NPIvals[k],4)))*k_compl) 
    if (k==length(NPIvals)) {params=cm_iv_apply(params,iv)} } else { # sets the "schedule" parameter to follow interventions in iv
  # version 2
  params$schedule[[k]]=list(parameter="contact",pops=numeric(),mode="multiply",
                          values=list(rep(NPIvals[k],4),rep(1,4)),times=NPI_phases[[k]])} }
### ### ### ### ### ### ### ### ### ###
# RUN SIMULATION
ptm<-proc.time(); run=cm_simulate(params,1); proc.time()-ptm 
# selvars=c("cases","subclinical","death_o");covidm_simul_agesep=fcn_covidm_df(run$dynamics,sel_vars=c("cases","subclinical","death_o","R"),params) 
# sum of age groups
compartm_types=list(case_vars=c("cases","subclinical","Ia","Ip","Is","S","R"),death_vars=c("D","death_o"))
dynamics_type=list(cumul=c("D","R","S"),incid=c("cases","subclinical","death_o"),preval=c("Ia","Ip","Is"))
covidm_simul=fcn_covidm_process_output(run$dynamics,filter_vars=c("E","foi","cases_reported"),
      compartm_types,dynamics_type,populval=mogadishu_popul,params)
# PLOT
npi_df=data.frame(on_off=c("on","off"),NPI_phases) %>% pivot_longer(!on_off) %>% mutate(date=as.Date(value)) %>% filter(on_off %in% "on")
seeding_df=data.frame(seed_date=unique(covidm_simul$date)[unique(params$pop[[1]]$seed_times)]) %>% summarise(min=min(seed_date),max=max(seed_date))
ggplot(subset(covidm_simul,!dynam_type %in% "preval")) + geom_area(aes(x=date,y=value,fill=compartment),color="black",size=0.3) +
  facet_wrap(dynam_type~compartm_type,scales="free") + theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) +
  scale_x_date(limits=as.Date(c("2019-11-30","2020-08-15")),date_breaks="2 weeks",expand=expansion(0,0)) + ylab("number") + 
  scale_y_continuous(expand=expansion(0.01,0)) + geom_vline(data=npi_df,aes(xintercept=date,color=name),size=1.1) + # linetype="dashed"
  geom_rect(data=seeding_df,aes(xmin=min,xmax=max,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.4)
# SAVE
# ggsave(paste0("simul_output/somalia_output/covidm",cm_version,"_output.png"),width=30,height=20,units="cm")

# PLOT incident deaths with data
# out_bdr_daily_estimates %>% select(date,new_graves_best_ipol,daily_baseline_subtr,rollmeanweek)
fitting_date_window=as.Date(c("2019-01-15","2020-08-01"))
fcn_covidm_singlesim_error(covidm_simul,out_bdr_daily_estimates,fitting_date_window)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### MCMC --------------
# define fitting parameters
# linear regression of ifr so we can generate new vals
somalia_agegroups_IFR = somalia_agegroups_IFR %>% mutate(log_ifr=log(ifr_mean),logit_ifr=log(ifr_mean/(1-ifr_mean)))
# ggplot(somalia_agegroups_IFR,aes(x=agegroup_mean)) + geom_line(aes(y=log_ifr)) + geom_point(aes(y=log_ifr)) + 
#   geom_line(aes(y=logit_ifr),color="red") + geom_point(aes(y=logit_ifr),color="red")+ scale_x_continuous(breaks=2.5+(0:16)*5)
linregr=lm(logit_ifr~agegroup_mean,data=somalia_agegroups_IFR %>% select(agegroup_mean,logit_ifr) )
ggplot(somalia_agegroups_IFR %>% mutate(pred_ifr=
  inv.logit(linregr$coefficients["(Intercept)"] + linregr$coefficients["agegroup_mean"]*somalia_agegroups_IFR$agegroup_mean),
  shifted_ifr=inv.logit(-8 + linregr$coefficients["agegroup_mean"]*somalia_agegroups_IFR$agegroup_mean)) %>%
  select(agegroup_mean,ifr_mean,pred_ifr,shifted_ifr) %>% pivot_longer(!c(agegroup_mean)),aes(x=agegroup_mean,y=value*1e2,group=name,color=name)) +
  geom_line() + geom_point() + theme_bw() + standard_theme + scale_x_continuous(breaks=2.5+(0:16)*5) + scale_y_log10()
# predicted IFR: exp(-10.8 + 0.1*c(2.5+(0:14)*5,80.255))
# clinical fraction
# d_clinfr=0.3; params$pop[[1]]$y=fun_lin_approx_agedep_par(somalia_agegroups_IFR,min_val=0.35-d_clinfr,max_val=0.35+d_clinfr,rep_min=3,rep_max=5)
# plot(somalia_agegroups_IFR$agegroup_mean,params$pop[[1]]$y,type = "b")
# compliance with NPIs
### ### ### ### ### ### ### ###
# define parameters func, which interprets a proposal for the posterior distribution as a parameter set usable by the underlying model.
fitting_params <- c("R0_target","introd_date","seed_size","ifr_exp_intercept", "compliance")
pf <- function(parameters, x){ x = as.list(x); n_groups = length(parameters$pop[[1]]$size);
    # R0
    parameters$pop[[1]]$u=rep(0.0232,n_groups)*(x$R0_target/2.4)
    # seed size and introd date
    parameters$pop[[1]]$seed_times=rep(x$introd_date:139,each=x$seed_size)
  # IFR
    agegroupmeans=c(2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,80.255)
  parameters$processes<-list(cm_multinom_process("Ip",outcomes=data.table(death=inv.logit(x$ifr_exp_intercept+0.1*agegroupmeans)),
                                                 delays=data.table(death=ponset2death),report="o"))
  # compliance
  t_npi=list(c("2020-03-19","2020-05-28"),c("2020-05-29","2020-06-30"),c("2020-07-01","2020-09-24")); npi_vals=c(0.423,0.524,0.714)
    for (k in 1:length(npi_vals)) { if (k==1) {iv=cm_iv_build(parameters)}
      cm_iv_contact(iv, t_npi[[k]][1], t_npi[[k]][2], 1 - ( 1-as.numeric(rep(npi_vals[k],4)) )*x$compliance ) 
      if (k==length(npi_vals)) {parameters$pop[[1]]$schedule=NULL; parameters=cm_iv_apply(parameters,iv)} }
    return (parameters) }
# priors
# priors = list(u="N 0.1 0.025 T 0 0.2", t_intro="U 0 10")
priors=list(R0_target="N 3 1 T 1 5", introd_date="N 50 10 T 10 90",seed_size="U 1 5",ifr_exp_intercept="U -12 -6",compliance="U 0 1")
# data
# scaling by CDR
CDR_vals=c(0.3,0.4,0.5); scale_factor=(mogadishu_popul*CDR_vals[1]/1e4)/baseline_daily_burials
fitting_incidence <- data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
                          mutate(daily_baseline_subtr=round(daily_baseline_subtr*scale_factor)) %>% 
                    filter(date>=fitting_date_window[1]&date<=fitting_date_window[2]) %>% rename(new_deaths=daily_baseline_subtr))
# define likelihood function
likelihood = function(parameters, dynamics, data, x){
  inc = data; inc[, t := as.numeric(date - ymd(parameters$date0))];
  eval = merge(dynamics[compartment == "death_o", .(model_case = sum(value)), by = t], inc, by = "t");
  ll = sum(dpois(eval$new_deaths, lambda = pmax(0.1, eval$model_case), log = T)); return (ll) }
### ### #
# fitting
fit=cm_fit(base_parameters=params, priors = priors, parameters_func = pf, likelihood_func=likelihood,
  data = fitting_incidence, mcmc_burn_in=500, mcmc_samples=2000, mcmc_init_opt = F, opt_maxeval = 25 )

# show posteriors
# cm_plot_posterior(fit); cm_plot_pairwise(fit)
# histogram of posteriors
df_posteriors=fit$posterior %>% select(all_of(fitting_params)) %>% mutate(ifr_age_weighted_perc=1e2*sapply(ifr_exp_intercept,
  function(x) sum(exp(x+0.1*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc) ) )
# plot distribs
ggplot(subset(df_posteriors %>% mutate(n=1:nrow(fit$posterior)) %>% pivot_longer(!n), !name %in% "ifr_exp_intercept")) + 
  geom_histogram(aes(x=value,fill=name),bins=50,color="black",size=0.2) + facet_wrap(~name,scales="free") + 
  theme_bw() + standard_theme

# use posterior to generate sample dynamics from the model
dyn = cm_sample_fit(fit, 25)

# summarize these runs
sel_compartm="death_o"
summ=dyn[compartment=="death_o",.(death_o=sum(value)),by=.(t, run)]
summ=summ[,cm_mean_hdi(death_o),by=t] %>% mutate(date=as.Date(seq(fitting_date_window[1],fitting_date_window[2]+50,1)[t+1]))
# show model fit
ggplot(summ,aes(x=date)) + geom_ribbon(aes(ymin=lower,ymax=upper), fill="blue") + geom_line(aes(y=mean), colour = "blue") + 
 geom_line(data=out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% mutate(new_deaths=daily_baseline_subtr*scale_factor),
            aes(y=new_deaths),linetype="dashed") + standard_theme + theme_bw() + xlab("") + 
  theme(axis.text.x=element_text(vjust=0.5,angle = 90)) + 
  scale_x_date(date_breaks = "week",expand=expansion(0.0)) + scale_y_continuous(breaks=(0:20)*10,expand=expansion(0,0))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### simple SEIR model for Somalia satellite data ----------------------
# transmission parameters
# d_e (pre-infectious days) | d_p (presympt inf.) | d_c (Duration of symptomatic infectiousness in days) | d_death (timescale of deaths)
# d_e=4; d_p=1.5; d_c=3.5; d_s=5; d_death=22.5 
# TRANSMISSION parameter (same as susceptibility) (without popul normalisation, done in the fcn)
# beta_val=0.5 # ; R0=beta_val*(sympt_share*(d_p+d_c) + (1-sympt_share)*d_s)
# clinical fraction by age group
clinical_fraction=fun_lin_approx_agedep_par(agegroups=data.frame(N_tot),min_val=0.4-0.35,max_val=0.4+0.35,rep_min=5,rep_max=3)
# variables
# var names
var_categ_list=list(name_vars=c("S","E","I_P","I_C","I_S","R","cumul_sympt_inf","cumul_asympt_inf","cumul_presympt_inf","incid_E"),
  case_vars=c("S","I_P","I_C","I_S","R","R_recov","symptom_cases","new_sympt_inf","new_asympt_inf"),
  sel_vars=c("t","S","I_P","I_C","I_S","R_recov","D","new_E","new_presympt_inf","new_sympt_inf","new_asympt_inf","new_recov","new_deaths","new_deaths_poiss"),
  cumul_var=c("S","R","D","R_recov","cumul_sympt_inf","cumul_asympt_inf"), 
  delta_var=c("new_sympt_inf","new_asympt_inf","new_deaths","new_recov"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### single simulation ---------------------------------------------
# timecourse for one param set (R0=beta/(1/d_c)=beta*d_c)
# ifrval=0.14 # sympt_share=sum(clinical_fraction*N_tot/sum(N_tot)); 
ptm<-proc.time()
singlesimul <- fcn_somal_sir_singlesimul(var_categ_list,time_steps=seq(0,nrow(OxCGRT_input),by=0.25),
  num_params=c("beta"=0.69,"d_e"=4,"d_p"=1.5,"d_c"=3.5,"d_s"=5,"d_death"=15,"IFR"=0.14,"sympt_share"=sum(clinical_fraction*N_tot/sum(N_tot)),
  "asympt_infness"=0.45,"seed_size_val"=1,"seeding_duration"=60,"compliance"=0.1,"popul_tot"=mogadishu_popul),
  day0=as.Date("2019-12-15"),OxCGRT_input,fcn_flag="proc")
time_sim=proc.time()-ptm; time_sim
# PLOT different variables (takes some time)
fcn_plot_singlesim_separ_vars(singlesimul,OxCGRT_input,popultot=mogadishu_popul)
# SAVE
# filetag_tcourse=paste0(c(paste0(names(num_params)[c(1,4)],"_",round(num_params,3)[c(1,4)]*c(1,1e2)),"compliance",compliance_val,"day",introd_date),collapse="_")
# ggsave(paste0("simul_output/somalia_output/SIR_timecourse_",filetag,".png"),width=25,height=20,units="cm")
####
# compare to data (simulation scaled)
fcn_plot_singlesim_data_error_deaths(singlesimul,out_bdr_daily_estimates,selvars=c("new_deaths","new_deaths_poiss"),OxCGRT_input,
                                     baselineburialrate=baseline_daily_burials,death_rate_percap=0.5,n_per_persday=1e4,
                                     datelims=c("2019-11-15","2020-08-05"),popul=mogadishu_popul,"plot")
# attack rate, cumul death rate: # fcn_attack_death_rate(singlesimul[[1]],plot_flag="calc",y_breaks_lims=NA)
fcn_attack_death_rate(singlesimul,plot_flag="plot",y_breaks_lims=list(breaks=round(10^((-8:8)/4),2), limits=c(1,100)))

# PEAK cases and deaths
fcn_plot_peakdates(singlesimul,datelims=c("2020-04-05","2020-06-15"),datebreakval = "2 day",
                   c("new_sympt_inf","new_deaths","new_recov"),yminval=5e1,logflag="log") # ,
# c("new_E","new_presympt_inf","new_asympt_inf","new_sympt_inf")
# ggsave(paste0("simul_output/somalia_output/peak_cases_compartms.png"),width=25,height=15,units="cm")
# ggsave(paste0("simul_output/somalia_output/peak_recov_deaths_compartms.png"),width=25,height=15,units="cm")

### MCMC ---------------------------------------------

### parameter scan ---------------------------------------------
# param table of all permutations
scan_params=list(seedsize_scanvals=5,introd_date_scanvals=c(seq(as.Date("2019-09-03"),as.Date("2019-09-29"),7),
                seq(as.Date("2019-10-01"),as.Date("2019-11-10"),4),seq(as.Date("2019-11-10"),as.Date("2020-01-30"),14)),
    betaval=seq(0.24,0.29,0.01),compliance=c(0,5,10,15,20,25,26,28,30,32,35,37.5,40,42,50,70,100)/100, 
    IFR=c(1.5,2,3,4,4.5,4.75,5,5.25,5.5,6,7,8,10,15,20,30)/100) # # seed_dur=60
param_table=expand.grid(scan_params); colnames(param_table)=names(scan_params)
# truncate OxCGRT_input not to include data points after 13/Aug
# OxCGRT_input = subset(OxCGRT_input,date < "2020-08-20")
# RUN PARSCAN with parallel (foreach) or sequentially (for) | write "full" for k_lim to scan whole parameter table, otherwise a number
cores=detectCores(); cl<-makeCluster(cores[1]-1); registerDoParallel(cl); ptm<-proc.time()
# 1 calc ~ 0.07 (with parall)
df_ode_solution_scan=fcn_paramscan_parallel_seq("paral",k_lim="full",sir_varnames,var_categ_list,timesteps=seq(0,nrow(OxCGRT_input),by=0.5),
    param_table,transm_pars=c("d_c"=d_c,"d_e"=d_e,"d_death"=d_death,"seeding_duration"=seed_dur),clinical_fraction,OxCGRT_input,N_tot,
    subpopul=mogadishu_popul,filterdates=out_bdr_daily_estimates$date,
    parnames_number=c("beta","npi_compliance","IFR","seed_size","introd_date"))
time_sim=proc.time()-ptm; time_sim
stopCluster(cl)
# qsave(df_ode_solution_scan,paste0("simul_output/somalia_output/df_ode_solution_scan",prod(sapply(scan_params, length)),".qs"))
# df_ode_solution_scan=qread(paste0("simul_output/somalia_output/df_ode_solution_scan",prod(sapply(scan_params, length)),".qs"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### PARSCAN results (errors) --------------------------------------------
fitting_dates_lims=c(as.Date("2020-01-19"),as.Date("2020-08-14"))
# ESTIMATE of death rate for Mogadishu
# between 0.2-0.6 deaths/10^4 person-days. 0.4 deaths/10K ppl/day -> Mogadishu is 2.2M -> 88 deaths/day -> ~616/week
# CDR_val=0.2; # mogad_daily_death_rate=mogadishu_popul*CDR_val/1e4
CDR_vals=c(0.3,0.4,0.5); error_pred_data_comp=list()
for (k_CDR in 1:length(CDR_vals)) {
scale_factor=baseline_daily_burials/(mogadishu_popul*CDR_vals[k_CDR]/1e4)
# NOT fitting the period before the rise of burials begin -> as.Date("2020-01-19")
error_pred_data_comp[[k_CDR]]=left_join(subset(df_ode_solution_scan,date>fitting_dates_lims[1]&date<fitting_dates_lims[2]) %>% 
  rename(simul_value=value),out_bdr_daily_estimates %>% select(date,new_graves_best_ipol),by="date",suffix=c("_model","_satel_data")) %>%
  mutate(data_baseline_subtr=ifelse(new_graves_best_ipol-baseline_daily_burials>0,new_graves_best_ipol-baseline_daily_burials,0)) %>%
  group_by(parset_ID) %>% # introd_date,IFR,npi_compliance,beta,
  summarise(mae=sqrt(mean(abs(data_baseline_subtr-simul_value*scale_factor),na.rm=T)),
  # rmse=sqrt(mean((data_baseline_subtr-simul_value*scale_factor )^2,na.rm=T)),
  loglh_poiss=-sum(dpois(x=round(data_baseline_subtr)+1,lambda=round(simul_value*scale_factor)+1,log=TRUE) ),
  loglh_negbinom=-sum(dnbinom(x=round(data_baseline_subtr)+1,mu=round(simul_value*scale_factor)+1,log=TRUE,size=10))) %>%
  pivot_longer(cols=!c(parset_ID)) %>% mutate(CDR=CDR_vals[k_CDR]) # ,IFR,npi_compliance,beta,introd_date
}
error_pred_data_comp=bind_rows(error_pred_data_comp)
# compare likelihood with MAE
# err_types=data.frame(sapply(unique(error_pred_data_comp$name), function(x) {subset(error_pred_data_comp,grepl(x,name))$value}))
# ggplot(err_types) + geom_point(aes(x=mae,y=loglh_poiss)) + geom_point(aes(x=mae,y=loglh_negbinom),color="red") + theme_bw() + standard_theme +
#   scale_x_log10() + scale_y_log10() + ylab("neg. log-likelihood") + # limits=c(2e4,1e5) | limits=c(2e2,5e3)
#   xlab("mean abs error") + labs(title="MSE vs negative log-likelihood",subtitle="red=neg. binom., black=poisson")
# ggsave("simul_output/somalia_output/errors/logll_mse_comparison.png",width=32,height=20,units="cm")

### PLOT trajectories for best fits ------------------------
for (err_type in unique(error_pred_data_comp$name)){ # [grepl("negbin",unique(error_pred_data_comp$name))]
for (k_CDR in 1:length(CDR_vals)) {
k_thresh=1.5
df_plot=left_join(left_join(df_ode_solution_scan,data.frame(param_table %>% 
                                                              select(-seedsize_scanvals),parset_ID=1:nrow(param_table)),by="parset_ID"),
  subset(error_pred_data_comp,CDR==CDR_vals[k_CDR] & grepl(err_type,name)),by="parset_ID",suffix=c("","_error")) %>%
  mutate(R0=round(betaval*((1-IFR)*d_c+IFR*d_death),1)) %>% filter(value_error<k_thresh*min(value_error)) %>% 
  group_by(introd_date_scanvals,R0) %>% filter(value_error==min(value_error)) %>% group_by(introd_date_scanvals) %>%
  mutate(best_score=ifelse(value_error==min(value_error),2,1))
# using wesanderson palette
fcn_plot_best_fits(df_plot,out_bdr_daily_estimates,OxCGRT_input,seed_dur,baseline_daily_burials,CDR_vals[k_CDR],CDR_metric=1e4,
    fitting_dates_lims+c(-5,5),popul=mogadishu_popul,geomtext_vals=list(as.Date("2020-01-20"),0.88),ylimval=25,sizelimval=10,d_c,d_death,
      title_str="model output ~ f(introd. date)",colorpal=wes_palette("Zissou1",min(length(unique(df_plot$R0))/2,10),type="continuous"))
# SAVE
errfolder=paste0("simul_output/somalia_output/bestfits/",err_type,"/"); if (!dir.exists(errfolder)) {dir.create(errfolder)}
ggsave(paste0(errfolder,"trajs_err_",(k_thresh-1)*1e2,"pct_abovebest_CDR",CDR_vals[k_CDR],".png"),width=40,height=24,units="cm")
}
}

### PLOT best fits per parameter set -------------------------
g(x,dfplotdistr) %=% fcn_best_paramsets(error_pred_data_comp,param_table,threshold=0.01,partype="numerical") # "numerical" dates
ggplot(subset(dfplotdistr,grepl("loglh",err_type) & !name %in% "betaval")) + 
 geom_point(aes(x=CDR,y=meanval,group=interaction(CDR,err_type),color=factor(err_type)),pch="-",size=16,position=position_dodge2(width=0.03)) +
 geom_linerange(aes(x=CDR,ymin=minval,ymax=maxval,group=interaction(CDR,err_type),color=factor(err_type)),position=position_dodge2(width=0.03),
                 alpha=0.3,size=3) + facet_wrap(~name,scales="free") + theme_bw() + standard_theme + 
  xlab("") + ylab("parameter values") + labs(color="CDR") # scale_y_date(date_breaks="2 weeks") + 
  # theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) 
# save
ggsave(paste0("simul_output/somalia_output/bestfits/best_paramsets_distrib",
              ifelse(class(dfplotdistr$meanval)=="Date","_dates",""),".png"),width=10,height=10,units="cm")

### calculate attack rates for best paramsets -------------------------
g(x,dfplotdistr) %=% fcn_best_paramsets(error_pred_data_comp,param_table,threshold=0.001,partype="numerical"); df_attackrates_list=list()
for (k_best in 1:nrow(param_table[unique(x$parset_ID),])) {
  parset=param_table[unique(x$parset_ID)[k_best],]
  dfonesim=fcn_somal_sir_singlesimul(sir_varnames,var_categ_list,timesteps=seq(0,nrow(OxCGRT_input),by=0.25),
    day0=parset$introd_date_scanvals, num_params=c("beta"=parset$betaval,"gamma"=1/d_c,"sigma"=1/d_e,"d_death"=1/d_death,"IFR"=parset$IFR),
        clin_fract=clinical_fraction,seed_size_val=5,seeding_duration=60,OxCGRT_input,compliance=parset$compliance,
        mogadishu_popul*N_tot/sum(N_tot),xlimvals=c(),plot_flag="")[[1]]
df_attackrates_list[[k_best]]=cbind(fcn_attack_death_rate(dfonesim,plot_flag="calc",y_breaks_lims=list())[2,c("attack_rate","cumul_death_rate")],
      parset_ID=unique(x$parset_ID)[k_best]); print(k_best/length(unique(x$parset_ID))) }
df_attackrates=left_join(bind_rows(df_attackrates_list),
  param_table %>% rownames_to_column(var="parset_ID") %>% mutate(parset_ID=as.numeric(parset_ID))) %>% select(-seedsize_scanvals) %>%
  mutate(R0=betaval*((1-IFR)*d_c + IFR*d_death)) %>% pivot_longer(cols=-c(introd_date_scanvals,R0,betaval,compliance,IFR,parset_ID)) %>%
  mutate(IFR=1e2*IFR,value=1e2*value)
# plot
colorpal=wes_palette("Zissou1",10,type="continuous")# min(subset(df_attackrates,name %in% "attack_rate")$value)
ggplot(subset(df_attackrates,name %in% "attack_rate")) + geom_point(aes(x=introd_date_scanvals,y=R0,size=value,color=IFR),shape=16) + 
  scale_color_gradientn(colours=colorpal) + theme_bw() + standard_theme + scale_x_date(date_breaks="week") + xlab("introduction rate") +
  labs(size="attack rate (%)",color="IFR (%)")
ggsave(paste0("simul_output/somalia_output/bestfits/best_paramsets_attackrates.png"),width=15,height=10,units="cm")

### ### ### ### ### ### ### ### ### ### ### ###
### profile likelihood from grid search
profile_likelihood=left_join(error_pred_data_comp,param_table %>% 
  mutate(parset_ID=1:nrow(param_table),R0=round(betaval*((1-IFR)*d_c + IFR*d_death)/0.05)*0.05),by="parset_ID") %>%
  mutate(introd_date=as.numeric(introd_date_scanvals)) %>% filter(grepl("loglh",name)) %>%
  rename(error_type=name,error=value) %>% pivot_longer(cols=c(introd_date,IFR,compliance,R0)) %>% # seed_size
  group_by(error_type,name,value,CDR) %>% summarise(best_fit=min(error)) %>% 
  mutate(value=ifelse(name %in% c("IFR","npi_compliance"),value*1e2,value),value_date=as.Date(value,origin=as.Date("1970-01-01")))
# confidence intervals on parameters from profile likelihood 
ci_params=profile_likelihood %>% group_by(CDR,name,error_type) %>% # ,grepl("loglh",error_type)
  summarise(best_fit_all=min(best_fit),ci95=min(best_fit)+qchisq(0.95,1)/2,xminval=min(value),xmaxval=max(value))
# left_join with profile likelihood
profile_likelihood=left_join(profile_likelihood,ci_params %>% select(-c(xminval,xmaxval)),by=c("CDR","name","error_type")) %>%  
  mutate(within_CI=ifelse(best_fit<=ci95,T,F))

### ### ### ### ###
for (errtypeval in unique(profile_likelihood$error_type)) { # unique(profile_likelihood$error_type)
  df_plot=subset(profile_likelihood,grepl(errtypeval,error_type) & best_fit<6e2) #  
ggplot(df_plot,aes(x=value,y=best_fit,group=CDR)) + geom_line() + geom_point(aes(color=within_CI)) + scale_color_manual(values=c("black","red")) +
  facet_wrap(CDR~name,labeller=labeller(CDR=label_both),scales="free") + theme_bw() + standard_theme + theme(legend.position="top") + 
  xlab("parameter value") + ylab("sum(negLL)") + scale_y_log10() + scale_x_continuous(expand=expansion(0.01,0)) + 
  labs(title=paste0("Profile likelihood (",errtypeval,")"),subtitle="best fit for given parameter value", caption=paste0("dates:",paste0(gsub(
    "2020","'20",gsub("2019","'19",as.Date(unique(subset(profile_likelihood,name %in% "introd_date")$value),origin="1970-01-01"))),
    collapse=", ")),color="CI95") + geom_rect(data=subset(ci_params,grepl(errtypeval,error_type)),
    aes(xmin=-Inf,xmax=Inf,ymin=best_fit_all,ymax=ci95),inherit.aes=F, size=NA,fill="pink",alpha=0.4)
# SAVE
ggsave(paste0("simul_output/somalia_output/bestfits/",errtypeval,"/profile_likelihood_CDR.png"),width=32,height=26,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### optimisation --------------
# ODE in fcns file

# DATA | # baseline_daily_burials=min(newburials_weekly$new_graves_best_ipol)
fitting_data=data.frame(date=newburials_weekly$date, t=as.numeric(newburials_weekly$date)-as.numeric(newburials_weekly$date)[1],
            new_burials_weekly=(newburials_weekly$new_graves_best_ipol-baseline_daily_burials)*mogad_weekly_death_rate/baseline_daily_burials)
# fcn to calc error
seir_error=function(opt_params,params_list,data,err_type,npi_date){
  ode_sol <- euler(y=params_list$initvals,times=params_list$timesteps,func=seir_deaths_npi_interpol,parms=opt_params) %>%
   as.data.frame() %>% filter(time %% 1 ==0) %>% mutate(new_deaths=round(D-lag(D, default=D[1])),date=npi_date[time+1]) %>% 
    mutate(new_deaths_weekly=roll_sum(new_deaths,7,fill=NA,align="right")) %>% filter(date %in% data$date) 
# error
  data_simul=subset(left_join(ode_sol,fitting_data,by="date"), date>params_list$fitting_start_date)
if (grepl("sse",err_type)){ simul_error=sum((data_simul$new_deaths_weekly - data_simul$new_burials_weekly)^2)}
if (grepl("rmse",err_type)){ simul_error=sqrt(mean((data_simul$new_deaths_weekly - data_simul$new_burials_weekly)^2))}
if (grepl("poiss",err_type)){ simul_error=-sum(dpois(x=round(data_simul$new_burials_weekly),
                                                     lambda=round(data_simul$new_deaths_weekly),log=TRUE)) }
if (grepl("negbin",err_type)){simul_error=-sum(dnbinom(x=round(data_simul$new_burials_weekly), 
                                          mu=round(data_simul$new_deaths_weekly),log=TRUE,size=params_list$negbinom_size)) }
simul_error
}

# calc ERROR/likelihood
params_list=list("initvals"=c("S"=mogadishu_popul,"E"=0,"I_R"=0,"I_D"=0,"D"=0),"popul"=mogadishu_popul,
                 "timesteps"=seq(0,nrow(OxCGRT_input),0.25),"fitting_start_date"="2020-01-19","negbinom_size"=10,
                 "introd_date"=as.Date("2019-11-12"),"seeding_duration"=21)
# parameters to optimize
opt_params=c("beta"=0.3,"IFR_estim"=4/100,"seed_size"=17,"npi_compl"=0.1)
# time-varying inputs
# day0_num=as.Date("2019-11-12"); seeding_duration=21
timevar_signal<-data.frame(t=(1:nrow(OxCGRT_input))-1,date=OxCGRT_input$date,npi_index=(1-OxCGRT_input$OxCGRT_scaled)*OxCGRT_input$NPI_on,seeding=0)
timevar_signal$seeding[timevar_signal$date %in% params_list$introd_date:(params_list$introd_date + params_list$seeding_duration)]=1
input_npi <- approxfun(timevar_signal[,c("t","npi_index")]); input_seeding<-approxfun(timevar_signal[,c("t","seeding")])
# error value
seir_error(opt_params,params_list,fitting_data,err_type="negbin",npi_date=OxCGRT_input$date)

# Optimisation
# Nelder-Mead
optim_proc=fcn_extract_optimresults(capture.output(optim(par=opt_params,fn=seir_error,method='Nelder-Mead',data=fitting_data,
                  params_list=params_list,err_type="rmse",npi_date=OxCGRT_input$date,control=c(trace=2))),parnames=names(opt_params))
optim_output=optim(par=opt_params,fn=seir_error,method='Nelder-Mead',data=fitting_data,
                   params_list=params_list,err_type="rmse",npi_date=OxCGRT_input$date,control=c(trace=2))
# plot convergence
ggplot(optim_proc,aes(x=count,y=sse)) + geom_line()
# pre vs post-optim error
sapply(list(opt_params,optim_output$par),function(x){seir_error(x,params_list,fitting_data,err_type="negbin",npi_date=OxCGRT_input$date)})

# single simul comparison
# gamma=0.29; sigma=0.25; d_death=1/21; inv_timestep=1/0.25
# singlesimul <- fcn_somal_sir_singlesimul(sir_varnames,var_categ_list,timesteps,day0=day0_num,
#  num_params=c("beta"=array(optim_output$par["beta"]),"gamma"=0.29,"sigma"=0.25,"d_death"=1/21,"IFR"=array(optim_output$par["IFR_estim"])),
#  clin_fract=clinical_fraction,seed_size_val=array(optim_output$par["seed_size"]),seeding_duration=seeding_duration,
#  OxCGRT_input,compliance=optim_output$par["npi_compl"],(N_tot/sum(N_tot))*mogadishu_popul,xlimvals=c(1,350),"plot"); singlesimul[[2]]

# concat BEST FIT with data
optim_sim_data=left_join(subset(singlesimul[[1]][,c("date","value","name")], name %in% "new_deaths") %>% rename(new_deaths=value) %>%
          mutate(new_deaths_weekly=roll_sum(new_deaths,7,fill=NA,align="right")), fitting_data,by="date") %>% filter(!is.na(t))
# PLOT BEST FIT
ggplot(optim_sim_data) + geom_line(aes(x=date,y=new_burials_weekly),color="blue") + geom_line(aes(x=date,y=new_deaths_weekly)) + 
  geom_point(aes(x=date,y=new_burials_weekly),color="blue") + theme_bw() + standard_theme + scale_x_date(date_breaks = "2 weeks") + 
  labs(title=paste0("Fit after optimisation (",paste0(paste0(c("introd_date","seeding_duration"),"=",
                  c(as.character(params_list$introd_date),params_list$seeding_duration)),collapse=", "),")"),
       subtitle=paste0(paste0(names(optim_output$par),"=",round(optim_output$par,2)),collapse = ", "))
# BEST FIT
ggsave(paste0("simul_output/somalia_output/errors/best_fit_neldermeadopt.png"),width=32,height=18,units="cm") # _facetgrid

### MCMC fitting --------------
# remotes::install_github('sbfnk/fitR') # package from S Funk for MCMC
# library('fitR')
