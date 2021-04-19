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
OxCGRT_input=fcn_get_OxCGRT(OxCGRT_url,"Somalia") %>% mutate(OxCGRT_scaled_smoothed=roll_mean(OxCGRT_scaled,30,align="center",fill=NA))
ggplot(OxCGRT_input) + geom_line(aes(x=date,y=OxCGRT_scaled)) +theme_bw()+standard_theme + theme(axis.text.x=element_text(vjust=0.5)) +
  scale_x_date(limits=c(min(OxCGRT_input$date[OxCGRT_input$NPI_on>0])-5,max(out_bdr_daily_estimates$date)+30),breaks="week",
               expand=expansion(0.01,0))

### Somalia population, IFR ------------
somalia_agegroups_IFR=fcn_merge_ifr_above_age(left_join(fcn_load_age_str("Somalia",90),fcn_load_ifr("data/IFR_by_age_imperial.csv"),
                                    by=c("agegroup","agegroup_min")),75); somalia_agegroups_IFR$ifr_mean[1]=3e-6
somal_popul_tot=sum(somalia_agegroups_IFR$agegroupsize); mogadishu_popul=2.2e6

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# COVIDM
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
cm_force_rebuild=F; cm_build_verbose=T; cm_version=1; source(file.path(cm_path,"R","covidm.R"))
countryval="Somalia"; params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval),
                                                 date_start="2019-11-01",date_end="2020-11-01")
# set population: Somalia --> Mogadishu
params$pop[[1]]$name=countryval
params$pop[[1]]$size=somalia_agegroups_IFR$agegroupsize*(mogadishu_popul/sum(somalia_agegroups_IFR$agegroupsize))
### seeding ---
npi_on_day=min(OxCGRT_input$date[OxCGRT_input$NPI_on>0]); introd_day=as.Date("2019-12-01")
seeding_t_window=sapply(c(introd_day,introd_day),function(x) as.numeric(x-as.Date(params$date0)))
params$pop[[1]]$seed_times=rep(seeding_t_window[1]:seeding_t_window[2],each=2) # x new infections/day for n days
# infections start in individuals aged 20-50
params$pop[[1]]$dist_seed_ages=cm_age_coefficients(20,60,5*(0:length(params$pop[[1]]$size)))
### add death process to model ------
params$processes <- list(cm_multinom_process("Ip",outcomes=data.table(death=somalia_agegroups_IFR$ifr_mean), 
                                             delays=data.table(death=cm_delay_gamma(22, 22, 60, 1/4)$p), report="o"))
# suscept and clinical fraction age dependent
suscept_clinfract_posteriors<-read_csv("data/suscept_clinfract_posteriors_davies2010.csv") %>% 
  mutate(agegroup=factor(agegroup,levels=unique(agegroup)))
# ggplot(suscept_clinfract_posteriors,aes(x=agegroup,y=value,group=1)) + geom_line() + geom_point() + 
# facet_wrap(~name,scales="free") + theme_bw() + scale_y_continuous(breaks=(2:14)/20)
params$pop[[1]]$y=fun_lin_approx_agedep_par(agegroups=somalia_agegroups_IFR,min_val=0.25,max_val=0.7,rep_min=6,rep_max=2)
# change susceptibility to get R0
target_R0=2; params$pop[[1]]$u=c(rep(0.38,4),rep(0.8,12)); params$pop[[1]]$u=params$pop[[1]]$u*(target_R0/cm_calc_R0(params,1))
# R0: cm_calc_R0(params,1)
### ### ### ### ###
# NPIs
NPI_phases=list(first=c("2020-03-19","2020-06-30"),second=c("2020-07-01","2020-08-20"),
                third=c("2020-08-30","2020-10-08"),fourth=c("2020-10-09","2020-11-01"))
NPIvals=sapply(NPI_phases,function(x) mean(OxCGRT_input$OxCGRT_scaled[OxCGRT_input$date>as.Date(x)[1]&OxCGRT_input$date<as.Date(x)[2]]))
k_compl=0.4
for (k in 1:length(NPIvals)) { # setup for version 1
  if (cm_version==1) {if (k==1) {params$pop[[1]]$schedule=NULL; iv=cm_iv_build(params)} # sets up data structure for interventions
    cm_iv_contact(iv, NPI_phases[[k]][1], NPI_phases[[k]][2], 1-(1-as.numeric(rep(NPIvals[k],4)))*k_compl) 
    if (k==length(NPIvals)) {params=cm_iv_apply(params,iv)} } else { # sets the "schedule" parameter to follow interventions in iv
  # version 2
  params$schedule[[k]]=list(parameter="contact",pops=numeric(),mode="multiply",
                          values=list(rep(NPIvals[k],4),rep(1,4)),times=NPI_phases[[k]])} }
### ### ### ### ### ### ### ### ### ###
# RUN SIMULATION
ptm <- proc.time(); run=cm_simulate(params,1); proc.time()-ptm 
# covidm_simul_agesep=fcn_covidm_df(run$dynamics,sel_vars=c("cases","subclinical","death_o","R"),params) 
# sum of age groups
covidm_simul=fcn_covidm_process_output(run$dynamics,filter_vars=c("E","foi","cases_reported"),
      compartm_types=list(case_vars=c("cases","subclinical","Ia","Ip","Is","S","R"),death_vars=c("D","death_o")),
      dynamics_type=list(cumul=c("D","R","S"),incid=c("cases","subclinical","death_o"),preval=c("Ia","Ip","Is")),
      populval=mogadishu_popul,params)
# PLOT
# df for npi and seeding
npi_df=data.frame(on_off=c("on","off"),NPI_phases) %>% pivot_longer(!on_off) %>% mutate(date=as.Date(value)) %>% filter(on_off %in% "on")
seeding_df=data.frame(seed_date=unique(covidm_simul$date)[unique(params$pop[[1]]$seed_times)]) %>% summarise(min=min(seed_date),max=max(seed_date))
# make the plot
ggplot(subset(covidm_simul,!dynam_type %in% "preval")) + geom_area(aes(x=date,y=value,fill=compartment),color="black",size=0.3) +
  facet_wrap(dynam_type~compartm_type,scales="free") + theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) +
  scale_x_date(limits=as.Date(c(params$date0,params$time1)),date_breaks="2 weeks",expand=expansion(0,0)) + ylab("number") + 
  scale_y_continuous(expand=expansion(0.01,0)) + geom_vline(data=npi_df,aes(xintercept=date),size=0.2,linetype="dashed") + # 
  geom_rect(data=seeding_df,aes(xmin=min,xmax=max,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.4) + labs(fill="")
# SAVE
# ggsave(paste0("simul_output/somalia_output/covidm",cm_version,"_output.png"),width=30,height=20,units="cm")

# PLOT incident deaths with data
# out_bdr_daily_estimates %>% select(date,new_graves_best_ipol,daily_baseline_subtr,rollmeanweek)
fitting_date_window=as.Date(c("2020-01-15","2020-10-01"))
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
  inv.logit(linregr$coefficients["(Intercept)"] + linregr$coefficients["agegroup_mean"]*somalia_agegroups_IFR$agegroup_mean)) %>%
  # shifted_ifr=inv.logit(-8 + linregr$coefficients["agegroup_mean"]*somalia_agegroups_IFR$agegroup_mean)
  select(agegroup_mean,ifr_mean,pred_ifr) %>% rename(estimate=ifr_mean,fit=pred_ifr) %>% pivot_longer(!c(agegroup_mean)),
  aes(x=agegroup_mean,y=value*1e2,group=name,color=name)) + geom_line(size=1.05) + geom_point(size=2) +
  theme_bw() + standard_theme + scale_x_continuous(breaks=2.5+(0:16)*5) + labs(color="") + xlab("age (year)") + ylab("IFR %") +
  scale_y_log10(breaks=scales::trans_breaks("log10",function(x) 10^x),labels=scales::trans_format("log10",scales::math_format(10^.x)),
                limits=c(1e-4,10^1.1)) + theme(axis.text.x=element_text(vjust=0.5,size=12),axis.text.y=element_text(size=12))
### ### ###
# ggsave("simul_output/somalia_output/IFR_consensus_estimate_fit.png",width=15,height=10,units="cm")
# predicted IFR: exp(-10.8 + 0.1*c(2.5+(0:14)*5,80.255))
### ### ### ### ### ### ### ###
# define parameters func, which interprets a proposal for the posterior distribution as a parameter set usable by the underlying model.
fitting_params <- c("R0_target","introd_date","seed_size","ifr_logit_intercept", "compliance")

pf <- function(parameters, x){x=as.list(x); n_groups=length(parameters$pop[[1]]$size);
    # R0
    target_R0=2 # cm_calc_R0(params,1) # params$pop[[1]]$u
    parameters$pop[[1]]$u=c(rep(0.0145,4),rep(0.0305,12))*(x$R0_fit/target_R0)
    # seed size and introd date
    parameters$pop[[1]]$seed_times=rep(x$introd_date:x$introd_date,each=x$seed_size)
  # IFR
  agegroupmeans=c(2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,80.255); slope_val=0.1
  parameters$processes<-list(cm_multinom_process("Ip",outcomes=data.table(death=inv.logit(x$ifr_logit_intercept+slope_val*agegroupmeans)),
                                          delays=data.table(death=data.table(death=cm_delay_gamma(22, 22, 60, 1/4)$p)),report="o"))
  # compliance
  t_npi=list(c("2020-03-19","2020-05-28"),c("2020-05-29","2020-06-30"),c("2020-07-01","2020-09-24")); npi_vals=c(0.423,0.524,0.714)
    for (k in 1:length(npi_vals)) { if (k==1) {iv=cm_iv_build(parameters)}
      cm_iv_contact(iv, t_npi[[k]][1], t_npi[[k]][2], 1 - ( 1-as.numeric(rep(npi_vals[k],4)) )*x$compliance ) 
      if (k==length(npi_vals)) {parameters$pop[[1]]$schedule=NULL; parameters=cm_iv_apply(parameters,iv)} }
    return (parameters) }
# priors
priors=list(R0_fit="N 3 1 T 1 5", introd_date="N 50 10 T 10 90",seed_size="U 1 5",ifr_logit_intercept="U -12 -6",compliance="U 0 1")
# data
# scaling by CDR
# ESTIMATE of death rate for Mogadishu
# between 0.2-0.6 deaths/10^4 person-days. 0.4 deaths/10K ppl/day -> Mogadishu is 2.2M -> 88 deaths/day -> ~616/week
# CDR_val=0.2; # mogad_daily_death_rate=mogadishu_popul*CDR_val/1e4
CDR_vals=c(0.3,0.4,0.5); scale_factor=(mogadishu_popul*CDR_vals[1]/1e4)/baseline_daily_burials
fitting_incidence <- data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
                          mutate(daily_baseline_subtr=round(daily_baseline_subtr*scale_factor)) %>% 
                    filter(date>=fitting_date_window[1]&date<=fitting_date_window[2]) %>% rename(new_deaths=daily_baseline_subtr))
# define likelihood function
likelihood = function(parameters, dynamics, data, x){
  inc = data; inc[, t := as.numeric(date - ymd(parameters$date0))];
  eval = merge(dynamics[compartment == "death_o", .(model_case = sum(value)), by = t], inc, by = "t");
  ll = sum(dpois(eval$new_deaths, lambda = pmax(0.1, eval$model_case), log = T)); return (ll) }
### ### ###
# fitting
fit=cm_fit(base_parameters=params, priors = priors, parameters_func = pf, likelihood_func=likelihood,
  data = fitting_incidence, mcmc_burn_in=500, mcmc_samples=2000, mcmc_init_opt = F, opt_maxeval = 25 )

# show posteriors
# cm_plot_posterior(fit); cm_plot_pairwise(fit)
# histogram of posteriors
# fit$posterior 
df_posteriors = fit$posterior  %>% select(all_of(fitting_params)) %>% mutate(ifr_age_weighted_perc=1e2*sapply(ifr_logit_intercept,
  function(x) sum(inv.logit(x+0.1*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc) ) )
# IFR for ALL infections
df_posteriors = posterior_ifr_logit  %>% select(all_of(fitting_params)) %>% 
  mutate(ifr_sympt_perc=1e2*sapply(ifr_logit_intercept,
  function(x) sum(inv.logit(x+0.1*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc ) ),
  ifr_all_perc=1e2*sapply(ifr_logit_intercept,
        function(x) sum(inv.logit(x+0.1*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc*params$pop[[1]]$y ) ))

# plot distribs
ggplot(subset(df_posteriors %>% mutate(n=1:nrow(df_posteriors)) %>% pivot_longer(!n), !name %in% "ifr_logit_intercept")) + 
  geom_histogram(aes(x=value,fill=name),bins=50,color="black",size=0.2) + facet_wrap(~name,scales="free") + 
  theme_bw() + standard_theme
# SAVE
ggsave("simul_output/somalia_output/posteriors.png",width=30,height=18,units="cm")

# CI95
posterior_CI95 = df_posteriors %>% mutate(n=row_number()) %>% pivot_longer(!n) %>% group_by(name) %>% 
  summarise(ci95_low=quantile(value,probs=c(2.5,97.5)/1e2)[1],ci95_up=quantile(value,probs=c(2.5,97.5)/1e2)[2] )

ifr_estimates=cbind(somalia_agegroups_IFR %>% select(agegroup_mean,ifr_mean),
data.frame(sapply(as.numeric(array(subset(posterior_CI95,name %in% "ifr_logit_intercept")[,2:3])), function(x)
       inv.logit(x + 0.1*somalia_agegroups_IFR$agegroup_mean)*params$pop[[1]]$y )) %>% rename(ifr_lower=X1,ifr_upper=X2) )
# plot
ggplot(ifr_estimates %>% rename(`estimate from literature`=ifr_mean) %>% pivot_longer(!c(agegroup_mean)),aes(x=agegroup_mean,y=value*1e2,group=name,color=name)) +
  geom_line(size=1.05) + geom_point(size=2) + theme_bw() + standard_theme + scale_x_continuous(breaks=2.5+(0:16)*5) + labs(color="") +
  scale_y_log10(breaks=scales::trans_breaks("log10",function(x) 10^x),labels=scales::trans_format("log10",scales::math_format(10^.x)),
  limits=c(1e-4,10^2)) + theme(axis.text.x=element_text(vjust=0.5,size=12),axis.text.y=element_text(size=12)) + 
  xlab("age (year)") + ylab("IFR %")
# SAVE
ggsave("simul_output/somalia_output/ifr_mcmc_estimates.png",width=20,height=12,units="cm")
### ### ### ### ### 
# use posterior to generate sample dynamics from the model
dyn = cm_sample_fit(fit, 25)

# summarize these runs
sel_compartm="death_o"
summ=dyn[compartment=="death_o",.(death_o=sum(value)),by=.(t, run)]
summ=summ[,cm_mean_hdi(death_o),by=t] %>% mutate(date=as.Date(seq(fitting_incidence$date[1],fitting_incidence$date[1]+max(summ$t),1)[t+1]))
# show model fit
ggplot(summ) + geom_ribbon(aes(x=date,ymin=lower,ymax=upper), fill="blue") + geom_line(aes(x=date,y=mean), colour = "blue") + 
 geom_line(data=fitting_incidence %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),linetype="dashed") + 
  standard_theme + theme_bw() + xlab("") + theme(axis.text.x=element_text(vjust=0.5,angle = 90)) + 
  # scale_x_date(date_breaks = "week",expand=expansion(0.0)) + 
  scale_y_continuous(breaks=(0:20)*10,expand=expansion(0,0))
