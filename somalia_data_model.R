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
lapply(c("tidyverse","deSolve","qs","gtools","rstudioapi","wpp2019","countrycode","coronavirus","wesanderson","dttr2","RcppRoll",
         "scales","wpp2019"), library,character.only=TRUE) 
# detach("package:fitdistrplus", unload = TRUE); detach("package:MASS", unload = TRUE) # "foreach","parallel","doParallel"
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
  geom_line(aes(y=rollingmean)) + facet_wrap(~name,scales="free",nrow=2) + # geom_point(aes(y=cases),size=0.5,fill=NA,shape=1) + #
  scale_x_date(name="", limits=c(as.Date("2020-03-01"), NA_Date_),date_breaks="weeks",expand=expansion(0.01,0)) + 
  # geom_rect(aes(xmin=as.Date("2020-04-20"),xmax=as.Date("2020-04-24"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  # geom_rect(aes(xmin=as.Date("2020-05-10"),xmax=as.Date("2020-05-15"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  # geom_rect(aes(xmin=as.Date("2020-05-28"),xmax=as.Date("2020-06-01"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + 
  theme_bw() + standard_theme + ylab("") + scale_y_continuous(expand=expansion(0.02,0)) +
  labs(title="COVID19 in Somalia: confirmed cases and deaths (7-day rolling mean)",caption="source: JHU CCSE")
#' No lag in deaths compared to cases
# save
# tcourse_filename="simul_output/somalia/reported_cases_deaths.png" # _dots
# ggsave(tcourse_filename,width=30,height=18,units="cm") # _dots

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### SATELLITE DATA on cemeteries ----------------
#' satellite image data
burial_data=read_csv("data/somalia/Mogadishu_data/mogadishu_burial_analysis-main/out_bdr_daily_estimates.csv")
# pop_wp2015,pop_wp2020,new_graves_best_ipol,
ggplot(burial_data %>% select(date,br_wp2015,br_wp2015_base_s,br_wp2020,br_wp2020_base_s) %>% pivot_longer(!date)) + 
  geom_line(aes(x=date,y=value,group=name,color=name)) + # facet_wrap(~name,scales="free") + 
  scale_x_date(expand=expansion(0.01,0),date_breaks="2 month") + theme_bw() + standard_theme

# baseline of daily burials july-november 2019
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
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*5) + 
  ggtitle(expression("weekly excess burials compared to confirmed COVID19 deaths (weekly sum)"))
# SAVE
# ggsave("simul_output/somalia/satellite_burials_reported_deaths_weekly.png",units="cm",height=15,width=25)
# ggsave("simul_output/somalia/satellite_burials_reported_deaths_weekly_barplot.png",units="cm",height=15,width=25)

# Rt estimate: https://epiforecasts.io/covid/posts/national/somalia/

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Oxford Stringency Index
# truncate until a given timepoint
OxCGRT_url="https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"
OxCGRT_input=fcn_get_OxCGRT(OxCGRT_url,"Somalia") %>% mutate(OxCGRT_scaled_smoothed=roll_mean(OxCGRT_scaled,30,align="center",fill=NA))
# separate into 4 phases
NPI_phases=list(first=c("2020-03-19","2020-06-30"),second=c("2020-07-01","2020-08-29"),
                third=c("2020-08-30","2020-10-08"),fourth=c("2020-10-09","2020-11-01"))
NPIvals=sapply(NPI_phases,function(x) mean(OxCGRT_input$OxCGRT_scaled[OxCGRT_input$date>as.Date(x)[1]&OxCGRT_input$date<as.Date(x)[2]]))
npi_df=left_join(data.frame(t(data.frame(on_off=c("on","off"),NPI_phases))) %>% add_rownames(var="name") %>% 
                 filter(!name=="on_off") %>% rename(on=X1,off=X2) %>% mutate(on=as.Date(on),off=as.Date(off)),
                 data.frame(NPIvals) %>% rownames_to_column(var="name"),by="name") %>% mutate(name=factor(name,levels=unique(name)))
# plot
ggplot(OxCGRT_input) + geom_line(aes(x=date,y=OxCGRT_scaled)) +
  geom_segment(data=npi_df,aes(x=on,xend=off,y=NPIvals,yend=NPIvals,group=factor(name),color=factor(name)),size=2) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) + labs(color="NPI phases") +
  scale_x_date(limits=c(min(subset(OxCGRT_input,NPI_on>0)$date)-15,max(out_bdr_daily_estimates$date)+35),breaks="week",
               expand=expansion(0.01,0)) + scale_y_continuous(breaks=(0:10)/10) + 
  geom_vline(data=npi_df,aes(xintercept=on,color=name),linetype="dashed") + ylab("level of contacts if reduction ~ StringencyIndex")
# ggsave(paste0("simul_output/somalia/OxCGRT_input.png"),width=30,height=18,units="cm")

### Somalia population, IFR ------------
somalia_agegroups_IFR=fcn_merge_ifr_above_age(left_join(fcn_load_age_str("Somalia",n_year="2015",90),
                                            fcn_load_ifr("data/IFR_by_age_imperial.csv"),by=c("agegroup","agegroup_min")),75)
somalia_agegroups_IFR$ifr_mean[1]=3e-6; somal_popul_tot=sum(somalia_agegroups_IFR$agegroupsize); mogadishu_popul=2.2e6
# other IFR estimates
# from Sandmann 2021 cmmid paper
IFR_estimates_Sandmann2021<-read_csv("data/IFR_estimates_Sandmann2021.csv")
 if (any(IFR_estimates_Sandmann2021$value_percent>1)) {n_cols<-2:ncol(IFR_estimates_Sandmann2021)
   IFR_estimates_Sandmann2021[,n_cols]<-IFR_estimates_Sandmann2021[,n_cols]/1e2}
ggplot(data.frame(age=factor(somalia_agegroups_IFR$agegroup,levels=unique(somalia_agegroups_IFR$agegroup)),
                  imper=somalia_agegroups_IFR$ifr_mean,sandmann=IFR_estimates_Sandmann2021$value_percent) %>%
  pivot_longer(!age),aes(x=age,y=value,group=name,color=name)) + geom_line() + geom_point() + theme_bw() + scale_y_log10()
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# COVIDM
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
cm_force_rebuild=F; cm_build_verbose=T; cm_version=1; source(file.path(cm_path,"R","covidm.R"))
countryval="Somalia"; params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval),
                                                 date_start="2019-11-01",date_end="2020-11-01")
# set population: Somalia --> Mogadishu
params$pop[[1]]$name=countryval
params$pop[[1]]$size=somalia_agegroups_IFR$agegroupsize*(mogadishu_popul/sum(somalia_agegroups_IFR$agegroupsize))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### SEEDING ---
npi_on_day=min(OxCGRT_input$date[OxCGRT_input$NPI_on>0]); introd_day=as.Date("2019-12-01")
seeding_t_window=sapply(c(introd_day,introd_day),function(x) as.numeric(x-as.Date(params$date0)))
seedsize_per_day=50
params$pop[[1]]$seed_times=rep(seeding_t_window[1]:seeding_t_window[2],each=seedsize_per_day) # x new infections/day for n days
# infections start in individuals aged 20-50, 1 introd in each age group
params$pop[[1]]$dist_seed_ages=cm_age_coefficients(20,60,5*(0:length(params$pop[[1]]$size)))
### add death process to model ------
params$processes <- list(cm_multinom_process("Ip",outcomes=data.table(death=somalia_agegroups_IFR$ifr_mean), 
                                             delays=data.table(death=cm_delay_gamma(22,22,60,1/4)$p), report="o"))
# suscept and clinical fraction age dependent
suscept_clinfract_posteriors<-read_csv("data/suscept_clinfract_posteriors_davies2010.csv") %>% 
  mutate(agegroup=factor(agegroup,levels=unique(agegroup)))
ggplot(suscept_clinfract_posteriors %>% mutate(name=gsub("clin_fract","clinical fraction",name)),
    aes(x=agegroup,y=value,group=1)) + geom_line() + geom_point() + ylab("proportion/normalised value") +
    facet_wrap(~name,scales="free") + theme_bw() + scale_y_continuous(breaks=(2:20)/20)
params$pop[[1]]$y=fun_lin_approx_agedep_par(agegroups=somalia_agegroups_IFR,min_val=0.25,max_val=0.7,rep_min=6,rep_max=2)
# change susceptibility to get R0
# susceptibility estimates from (warwick): 
# susceptibility_warvick_model <- left_join(read_csv("data/susceptibility_warvick_model.csv") %>% mutate(value=value/max(value)),
#   fcn_load_age_str("Somalia",n_year="2015",100),by="agegroup") %>% mutate(agegroup_weight=agegroupsize/sum(agegroupsize),
#   agegroup_min=ifelse(agegroup_min>=75,75,agegroup_min),agegroup_max=ifelse(agegroup_max>=75,100,agegroup_max)) %>% 
#   group_by(agegroup_min) %>% summarise(agegroup=unique(agegroup)[1],value=sum(value*agegroup_weight/sum(agegroup_weight)),
#     agegroupsize=sum(agegroupsize),agegroup_min=unique(agegroup_min),agegroup_max=unique(agegroup_max)) %>% arrange(agegroup_min)
target_R0=1.8; params$pop[[1]]$u=c(rep(0.38,4),rep(0.8,12)) # susceptibility_warvick_model$value # 
params$pop[[1]]$u=params$pop[[1]]$u*(target_R0/cm_calc_R0(params,1))
# R0: cm_calc_R0(params,1)
### ### ### ### ###
# input NPIs
k_compl=0.3
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
# df for seeding
seeding_df=data.frame(seed_date=unique(covidm_simul$date)[unique(params$pop[[1]]$seed_times)]) %>% summarise(min=min(seed_date),max=max(seed_date))
# make the plot
ggplot(subset(covidm_simul,!dynam_type %in% "preval")) + geom_area(aes(x=date,y=value,fill=compartment),color="black",size=0.3) +
  facet_wrap(dynam_type~compartm_type,scales="free") + theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) +
  scale_x_date(limits=as.Date(c(params$date0,params$time1)),date_breaks="2 weeks",expand=expansion(0,0)) + ylab("number") + 
  scale_y_continuous(expand=expansion(0.01,0)) +geom_vline(data=npi_df,aes(xintercept=on,color=name),size=1,linetype="dashed") +
  geom_rect(data=seeding_df,aes(xmin=min,xmax=max,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.4) + labs(fill="")
# SAVE
# ggsave(paste0("simul_output/somalia/covidm",cm_version,"_suscept_warvick","_output.png"),width=30,height=20,units="cm")

# PLOT incident deaths with data
# out_bdr_daily_estimates %>% select(date,new_graves_best_ipol,daily_baseline_subtr,rollmeanweek)
fitting_date_window=as.Date(c("2020-01-15","2020-10-01"))
fcn_covidm_singlesim_error(covidm_simul,introd_day,seedsize_per_day,out_bdr_daily_estimates,fitting_date_window)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### MCMC --------------
# define fitting parameters
# IFR: linear regression -> generate new values
somalia_agegroups_IFR = somalia_agegroups_IFR %>% mutate(log_ifr=log(ifr_mean),logit_ifr=log(ifr_mean/(1-ifr_mean)))
# ggplot(somalia_agegroups_IFR,aes(x=agegroup_mean)) + geom_line(aes(y=log_ifr)) + geom_point(aes(y=log_ifr)) + 
#   geom_line(aes(y=logit_ifr),color="red") + geom_point(aes(y=logit_ifr),color="red")+ scale_x_continuous(breaks=2.5+(0:16)*5)
linregr=lm(logit_ifr~agegroup_mean,data=somalia_agegroups_IFR %>% select(agegroup_mean,logit_ifr) )
ggplot(somalia_agegroups_IFR %>% 
 mutate( #uk_estim=IFR_estimates_Sandmann2021$value_percent,
  pred_ifr=inv.logit(linregr$coefficients["(Intercept)"] + linregr$coefficients["agegroup_mean"]*somalia_agegroups_IFR$agegroup_mean)) %>%
 select(agegroup_mean,ifr_mean,pred_ifr) %>% rename(estimate=ifr_mean,fit=pred_ifr) %>% pivot_longer(!c(agegroup_mean)),
 aes(x=agegroup_mean,y=value*1e2,group=name,color=name)) + geom_line(size=1.05) + geom_point(size=2) +
 theme_bw() + standard_theme + scale_x_continuous(breaks=2.5+(0:16)*5) + labs(color="") + xlab("age (year)") + ylab("IFR %") +
 theme(axis.text.x=element_text(vjust=0.5,size=12),axis.text.y=element_text(size=12)) +
 scale_y_log10(limits=c(1e-4,10^1.1),breaks=scales::trans_breaks("log10",function(x) 10^x))
# ggsave("simul_output/somalia/IFR_consensus_estimate_fit.png",width=15,height=10,units="cm")
# predicted IFR: exp(-10.8 + 0.1*c(2.5+(0:14)*5,80.255))
### ### ### ### ### ### ### ###
# define parameters func, which interprets a proposal for the posterior distribution as a parameter set usable by the underlying model.
fitting_params <- c("R0_fit","introd_date","seed_size","ifr_logit_intercept", "compliance")
pf <- function(parameters, x){x=as.list(x); n_groups=length(parameters$pop[[1]]$size);
    # R0
    target_R0=2 # cm_calc_R0(params,1) # params$pop[[1]]$u
    parameters$pop[[1]]$u=c(rep(0.0145,4),rep(0.0305,12))*(x$R0_fit/target_R0)
    # seed size and introd date
    parameters$pop[[1]]$seed_times=rep(x$introd_date:x$introd_date,each=x$seed_size)
  # IFR
  agegroupmeans=c(2.5,7.5,12.5,17.5,22.5,27.5,32.5,37.5,42.5,47.5,52.5,57.5,62.5,67.5,72.5,80.255)
  slope_val=as.numeric(linregr$coefficients["agegroup_mean"])
  parameters$processes<-list(cm_multinom_process("Ip",outcomes=data.table(death=inv.logit(x$ifr_logit_intercept+slope_val*agegroupmeans)),
                                          delays=data.table(death=data.table(death=cm_delay_gamma(22,22,60,1/4)$p)),report="o"))
  # compliance
  t_npi=list(first=c("2020-03-19","2020-06-30"),second=c("2020-07-01","2020-08-29"),
              third=c("2020-08-30","2020-10-08"),fourth=c("2020-10-09","2020-11-01")); npi_vals=c(0.455,0.736,0.593,0.624)
    for (k in 1:length(npi_vals)) { if (k==1) {iv=cm_iv_build(parameters)}
      cm_iv_contact(iv, t_npi[[k]][1], t_npi[[k]][2], 1 - (1-as.numeric(rep(npi_vals[k],4)) )*x$compliance ) 
      if (k==length(npi_vals)) {parameters$pop[[1]]$schedule=NULL; parameters=cm_iv_apply(parameters,iv)} }
    return (parameters) }
# priors
priors=list(R0_fit="N 3 1 T 1 5", introd_date="N 50 10 T 10 90",seed_size="U 1 5",ifr_logit_intercept="U -12 -6",compliance="U 0 1")
# define likelihood function
likelihood = function(parameters, dynamics, data, x){
  inc = data; inc[, t := as.numeric(date - ymd(parameters$date0))]
  eval = merge(dynamics[compartment == "death_o", .(model_case = sum(value)), by = t], inc, by = "t");
  ll = sum(dpois(eval$new_deaths, lambda = pmax(0.1, eval$model_case), log = T)); return (ll) }
###  fitting data - scaling by CDR ----------------
# ESTIMATE of Mogadishu death rate: 0.2-0.6 deaths/10^4 person-days. 0.4 deaths/(10K ppl*day) -> Mogadishu=2.2M -> 616 deaths/week
# CDR_val=0.2; # mogad_daily_death_rate=mogadishu_popul*CDR_val/1e4
CDR_vals=c(0.0422654,0.1,0.2,0.4); scale_factor=(mogadishu_popul*CDR_vals[1]/1e4)/baseline_daily_burials
fitting_incidence <- data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
                mutate(daily_baseline_subtr=round(daily_baseline_subtr*scale_factor)) %>% 
                filter(date>=fitting_date_window[1]&date<=fitting_date_window[2]) %>% rename(new_deaths=daily_baseline_subtr))
### fitting -------------------
fit=cm_fit(base_parameters=params, priors=priors, parameters_func=pf, likelihood_func=likelihood,
  data=fitting_incidence, mcmc_burn_in=500, mcmc_samples=2000, mcmc_init_opt=F, opt_maxeval=25 )

### ### ### ### ### ### ### ### ### ###
# load RDS file with several fits (different CDR values)
mcmc_filename="simul_output/somalia/fits_death_scaling.rds"; fits_death_scaling <- readRDS(mcmc_filename)
foldertag=gsub("mcmc_","",paste0(paste0(names(fits_death_scaling[[1]]$options[c("mcmc_burn_in","mcmc_samples")]),
       as.numeric(fits_death_scaling[[1]]$options[c("mcmc_burn_in","mcmc_samples")]),collapse = "_"),"_fittingwindow_",
  paste0(c(min(fits_death_scaling[[1]]$data$date),max(fits_death_scaling[[1]]$data$date)),collapse = "_")))
mcmc_foldername=paste0("simul_output/somalia/fit_",foldertag); dir.create(mcmc_foldername); file.copy(mcmc_filename,mcmc_foldername)
# CDR_vals=c(baseline_daily_burials*1e4/mogadishu_popul,0.1,0.2,0.4); slope_val=linregr$coefficients["agegroup_mean"]
# fitting_date_window
for (k in 1:length(CDR_vals)) {
  df_posteriors=fits_death_scaling[[k]]$posterior %>% select(chain,trial,lp,all_of(fitting_params)) %>% 
    mutate(`IFR sympt. infections (%)`=1e2*sapply(ifr_logit_intercept,
    function(x) sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc) ),
    `IFR all infections (%)`=1e2*sapply(ifr_logit_intercept,function(x) 
      sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc*params$pop[[1]]$y))) %>%
    mutate(CDR=round(CDR_vals[k],3))
  if (k==1){df_posteriors_comb=df_posteriors} else {df_posteriors_comb = rbind(df_posteriors_comb,df_posteriors)} }

# extract prior
# cm_evaluate_distribution(fits_death_scaling[[1]]$priors$R0_fit)
# cm_evaluate_distribution(fits_death_scaling[[1]]$priors$seed_size)
# # mean of prior
# sum(cm_evaluate_distribution(fits_death_scaling[[1]]$priors$seed_size)$x*
#       cm_evaluate_distribution(fits_death_scaling[[1]]$priors$seed_size)$p)/sum(
#         cm_evaluate_distribution(fits_death_scaling[[1]]$priors$seed_size)$p)

# calculate mean, median, CIs
posterior_CI95 = df_posteriors_comb %>% mutate(n=row_number()) %>% pivot_longer(!c(n,CDR,chain,trial,lp)) %>% 
  group_by(name,CDR) %>% summarise(mean=mean(value),median=median(value),
            ci95_low=quantile(value,probs=c(2.5,97.5)/1e2)[1],ci95_up=quantile(value,probs=c(2.5,97.5)/1e2)[2],
            ci50_low=quantile(value,probs=c(25,75)/1e2)[1],ci50_up=quantile(value,probs=c(25,75)/1e2)[2] ) %>%
  mutate(name=ifelse(name %in% "introd_date","introduction (days after 01/Nov/2019)",name),
         name=ifelse(name %in% "compliance","NPI compliance (0 to 1)",name))
# plot median and CIs
ggplot(posterior_CI95 %>% filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")),
       aes(x=factor(CDR),group=CDR,color=factor(CDR))) + scale_y_continuous(expand=expansion(0.1,0)) +
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge2(width=0.25),alpha=0.3,size=3) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge2(width=0.25),alpha=0.5,size=3) +
  geom_point(aes(y=median),pch="-",size=10,color="black") + facet_wrap(~name,scales="free") + theme_bw() + standard_theme + xlab("") + 
  ylab("mean (CI50, CI95)") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + labs(color="CDR",
  caption=paste0("Burn-in: ",fits_death_scaling[[1]]$options$mcmc_burn_in,", Iterations: ",fits_death_scaling[[1]]$options$mcmc_iter," (MCMC)")) +
  geom_text(aes(x=factor(CDR),y=median,label=round(median,2)),nudge_x=0.3,color="black",size=3.5)
# SAVE
ggsave(paste0(mcmc_foldername,"/posteriors_mean_CIs_CDRscan.png"),width=30,height=18,units="cm")

# plot IFR estimates
i_col=3:6 
ifr_estimates = bind_rows(lapply(1:length(CDR_vals), function(k) cbind(somalia_agegroups_IFR %>% 
  select(agegroup_mean),CDR=round(CDR_vals[k],3),
  sapply(as.numeric(array(subset(posterior_CI95,name %in% "ifr_logit_intercept" & CDR==round(CDR_vals[k],3))[,i_col])), function(x) 
 inv.logit(x+linregr$coefficients["agegroup_mean"]*somalia_agegroups_IFR$agegroup_mean)*params$pop[[1]]$y)))) %>% mutate(datatype="fit")
colnames(ifr_estimates)[i_col] <- colnames(posterior_CI95)[i_col]
ifr_estimates=bind_rows(ifr_estimates,somalia_agegroups_IFR %>% select(agegroup_mean,ifr_mean) %>% 
        mutate(ifr_mean=ifr_mean*params$pop[[1]]$y) %>% rename(median=ifr_mean) %>% 
          mutate(datatype="estimate from data" ) ) %>% mutate(CDR=ifelse(is.na(CDR),"estimate from literature",CDR))
# plot
ggplot(ifr_estimates,aes(x=agegroup_mean)) + 
  geom_line(aes(y=median*1e2,group=CDR,color=factor(CDR),linetype=factor(datatype)),size=1.05) + 
  geom_point(aes(y=median*1e2,group=CDR,color=factor(CDR),linetype=factor(datatype)),size=2) + labs(color="CDR",linetype="",fill="CDR") +
  geom_ribbon(aes(ymin=ci95_low*1e2,ymax=ci95_up*1e2,group=CDR,fill=factor(CDR)),alpha=0.2) + theme_bw() + standard_theme + 
  scale_x_continuous(breaks=2.5+(0:16)*5) + theme(axis.text.x=element_text(vjust=0.5,size=12),axis.text.y=element_text(size=12)) + 
  scale_y_log10(breaks=scales::trans_breaks("log10",function(x) 10^x,n=12),labels=scales::trans_format("log10",scales::math_format(10^.x))) +
  xlab("median age per age group (year)") + ylab("IFR %")
# SAVE
ggsave(paste0(mcmc_foldername,"/ifr_mcmc_estimates_CDRscan.png"),width=20,height=12,units="cm")

# plot traces from MCMC
ggplot(df_posteriors_comb %>% select(!`IFR sympt. infections (%)`,`IFR all infections (%)`,CDR) %>% 
         pivot_longer(!c(chain,trial,lp,CDR))) + geom_line(aes(x=trial,y=value,group=chain,color=factor(chain))) + 
  facet_grid(name~CDR,scales="free",labeller=labeller(CDR=label_both)) + theme_bw() + standard_theme + labs(color="chains")
# SAVE
ggsave(paste0(mcmc_foldername,"/MCMC_convergence.png"),width=30,height=20,units="cm")

### plot PRIORS + posteriors
for (k in 1:length(unique(df_posteriors_comb$CDR))) {
  cm_plot_posterior_mod(fits_death_scaling[[k]],plot_params = list(n_bin=50,line_width=0.7,cdr_val=CDR_vals[k]))
  ggsave(paste0(mcmc_foldername,"/priors_posteriors_CDR_",round(CDR_vals[k],3),".png"),width=30,height=20,units="cm") }

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# use posterior to generate sample dynamics from the model
for (k in 1: length(fits_death_scaling)){ dyn=cm_sample_fit(fits_death_scaling[[k]], 25) %>% mutate(CDR=CDR_vals[k])
if (k==1) {dyn_all=dyn} else {dyn_all=rbind(dyn_all,dyn)}; print(k) }

# summarize these runs
summ = dyn_all %>% filter(compartment %in% c("S","death_o")) %>% group_by(t,run,compartment,CDR) %>% summarise(value=sum(value)) %>%
  group_by(run,compartment,CDR) %>% mutate(value=ifelse(compartment %in% "S",(value[t==0]-value)/value[t==0],value), # /value[t==0]
         compartment=ifelse(compartment %in% "S","attack_rate",as.character(compartment))) %>% group_by(t,compartment,CDR) %>% 
  summarise(lower=hdi(value)[[1]],upper=hdi(value)[[2]],mean=mean(value) ) %>% 
  mutate(date=as.Date(seq(as.Date(params$date0),as.Date(params$date0)+max(t),1)[t+1]),CDR=round(CDR,5))

# fitting data (deaths incidence)
# CDR_vals=c(0.0422654,0.1,0.2,0.4) # scale_factor=(mogadishu_popul*CDR_vals[k]/1e4)/baseline_daily_burials
fitting_incidence_modelcompare = bind_rows(lapply((mogadishu_popul*CDR_vals/1e4)/baseline_daily_burials, function(x)
  data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
  mutate(daily_baseline_subtr=round(daily_baseline_subtr*x),CDR=x*1e4*baseline_daily_burials/mogadishu_popul,CDR=round(CDR,5),
  date_within_fitting_t=ifelse(date>=min(fits_death_scaling[[1]]$data$date)&date<=max(fits_death_scaling[[1]]$data$date),TRUE,FALSE)) %>%
           rename(new_deaths=daily_baseline_subtr)) ))
# calculate likelihood, deviance, DIC
logllk_values <- right_join(summ %>% filter(compartment=="death_o") %>% select(!c(lower,upper)),
                            subset(fitting_incidence_modelcompare, date_within_fitting_t) %>% select(!date_within_fitting_t),by=c("CDR","date")) %>% 
  mutate(logllk=dpois(new_deaths,lambda=mean,log=T)) %>% group_by(CDR) %>% 
  summarise(sum_logllk=sum(logllk),deviance=-2*sum_logllk,d_p=var(-2*logllk)/2,DIC=deviance+d_p)
# sum(dpois(eval$new_deaths, lambda = pmax(0.1, eval$model_case), log = T))

# show model fit
ggplot(summ %>% filter(compartment=="death_o")) + 
 geom_line(aes(x=date,y=mean), colour="blue") + geom_ribbon(aes(x=date,ymin=lower,ymax=upper), fill="blue",alpha=0.2) + 
 geom_line(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),linetype="dashed",size=0.3) +
 geom_point(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),fill=NA,shape=1) +
 facet_wrap(~CDR,scales="free",labeller=labeller(CDR=label_both)) + 
  geom_rect(data=fits_death_scaling[[1]]$data %>% summarise(min_date=min(date),max_date=max(date)),
            aes(xmin=min_date,xmax=max_date,ymin=-Inf,ymax=Inf),fill="pink",alpha=0.2,linetype="dashed",color="red") + 
  xlab("") + ylab("deaths (daily)") + standard_theme + theme_bw() + theme(axis.text.x=element_text(vjust=0.5,angle=90)) + 
 scale_x_date(date_breaks="week",limits=as.Date(c("2020-01-01","2020-10-15")),expand=expansion(0.0)) + 
 scale_y_continuous(breaks=(0:20)*10,expand=expansion(0,0))
# save
ggsave(paste0(mcmc_foldername,"/dynamics_fit_deaths_CDR_scan.png"),width=30,height=16,units="cm")

# attack rate plot
p <- ggplot(summ %>% filter(compartment=="attack_rate")) + 
  geom_line(aes(x=date,y=mean,group=CDR,color=factor(CDR)),size=1.2) + facet_wrap(~CDR,nrow=length(CDR_vals),labeller=labeller(CDR=label_both)) +
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,group=CDR,fill=factor(CDR)),alpha=0.2) + labs(color="CDR",fill="CDR") + 
  standard_theme + theme_bw() + xlab("") + ylab("cumulative attack rate") + theme(axis.text.x=element_text(vjust=0.5,angle=90)) +
  scale_x_date(date_breaks="week",limits=as.Date(c("2020-01-10","2020-11-01")),expand=expansion(0.0)) + 
  scale_y_continuous(breaks=(0:20)/20,expand=expansion(0,0)); p
# SAVE
att_rate_filename<-paste0(mcmc_foldername,"/dynamics_cumulattackrate_deaths_CDRscan.png")
if (length(p$facet$params$facets)>0) {att_rate_filename <- gsub(".png","_faceted.png",att_rate_filename)}
ggsave(att_rate_filename,width=22,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ACLED data on deaths from polit violence/terrorism
somalia_acled_data <- read_csv("data/somalia/terrorism_data/2018-04-23-2021-04-28-Somalia.csv")
somalia_acled_fatalities <- somalia_acled_data %>% group_by(event_date,admin1) %>% 
  summarise(date=unique(dmy(event_date)),fatalities=sum(fatalities),year=unique(year)) %>% ungroup() %>% group_by(admin1) %>% 
  complete(date=seq.Date(min(date),max(date),by="day")) %>% mutate(fatalities_missing_as_0=ifelse(is.na(fatalities),0,fatalities),
  rollingmean=roll_mean(fatalities_missing_as_0,7,align="center",fill=NA,na.rm=T),
  rollingsum=roll_sum(fatalities_missing_as_0,7,fill=NA,align="right")) %>% mutate(year=year(date),month=month(date),week=week(date))
# fitting_date_window=as.Date(c("2020-01-15","2020-10-01"))
p<-ggplot(subset(somalia_acled_fatalities, admin1=="Banadir" & date>=as.Date("2019-01-01") & date<=as.Date("2021-01-01"))) + # 
  geom_bar(aes(x=date,y=fatalities),stat="identity") + geom_line(aes(x=date,y=rollingmean)) + # facet_wrap(~admin1,scales="free") + 
  ggtitle("rolling 7-day mean of deaths in political violence/protest") +
  geom_rect(aes(xmin=as.Date("2020-01-15"),xmax=as.Date("2020-10-01"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,size=8)) +
  scale_x_date(date_breaks="month",expand=expansion(0,0)) + scale_y_continuous(expand=expansion(0.01,0)); p # breaks=(0:20),
# SAVE
acled_plot_name <- "simul_output/somalia/acled_banadir_rollingmean_fatalities.png"
if (length(p$facet$params$facets)>0) {acled_plot_name<-gsub("banadir","somalia",acled_plot_name)}
ggsave(acled_plot_name,width=30,height=20,units="cm")

# monthly
acled_monthly_fatalities=somalia_acled_fatalities %>% group_by(year,month,admin1) %>% 
  summarise(year=unique(year),date=min(date),fatalities=sum(fatalities_missing_as_0)) %>% 
  mutate(year_month=ifelse(nchar(as.character(month))==1,paste0(year,"/0",month),paste0(year,"/",month))) %>%
  mutate(year_month=factor(year_month,levels=unique(year_month)))
ggplot(subset(acled_monthly_fatalities,date<=as.Date("2020-11-01") & date>as.Date("2018-11-01")),aes(x=year_month,y=fatalities)) + geom_bar(stat="identity") + facet_wrap(~admin1,scales = "free") +
  geom_rect(aes(xmin="2020/01",xmax="2020/10",ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + ylab("deaths/month") + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,size=6))
# SAVE
ggsave("simul_output/somalia/acled_somalia_monthly_fatalities.png",width=30,height=18,units="cm")

# compare to burials
acled_burial_comparison = left_join(subset(somalia_acled_fatalities,admin1=="Banadir") %>% 
                    select(admin1,date,fatalities_missing_as_0,rollingmean) %>% rename(fatalities=fatalities_missing_as_0), 
  burial_data[!rowSums(is.na(burial_data))==(ncol(burial_data)-1),
  !colSums(is.na(burial_data))==nrow(burial_data)] %>% select(date,new_graves_best_ipol),by="date") 

# weekly monthly compare
weekly_acled_burial_comparison = acled_burial_comparison %>% mutate(week=week(date),year=year(date)) %>% group_by(year,week) %>% 
  summarise(date=min(date),fatalities=sum(fatalities),rollingmean=sum(rollingmean),new_graves_best_ipol=sum(new_graves_best_ipol)) %>%
  mutate(year_week=factor(paste0(year,"/",week),levels = unique(paste0(year,"/",week))))
# WEEKLY compare plot
ggplot(subset(weekly_acled_burial_comparison,date>=as.Date("2019-01-01") & date<=as.Date("2020-10-01")) %>% 
         select(!rollingmean) %>% pivot_longer(!c(date,year,date,week,year_week))) + #
  geom_bar(aes(x=year_week,y=value,group=name,fill=name),stat="identity",position=position_dodge(width=0.75),size=0.2) +
  scale_fill_discrete(labels=c("deaths due to political violence (ACLED)","burials")) +
  # geom_point(aes(x=year_week,y=value,group=name,color=name),pch="-",size=12) + 
  # geom_vline(xintercept=(1:length(unique(weekly_acled_burial_comparison$year_week)))-0.5,size=0.1,linetype="dashed") +
  # scale_color_discrete(labels=c("deaths due to political violence (ACLED)","burials")) +
  theme_bw() + standard_theme + theme(axis.text.x = element_text(vjust=0.5),legend.position="bottom") + labs(fill="",color="") +
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*10) + xlab("year/week") + ylab("number/week")
# SAVE
ggsave("simul_output/somalia/ACLED_data/acled_banadir_burials_comparison_WEEKLY.png",width=30,height=16,units="cm")

# daily compare plot
ggplot(acled_burial_comparison %>% pivot_longer(!c(admin1,date)) %>% filter(date<=as.Date("2020-10-01")),aes(x=date)) + # geom_line() +
  geom_bar(aes(y=fatalities,color="ACLED: daily deaths (armed conflicts/terrorism)"),stat="identity",alpha=0.1) +  
  geom_line(aes(y=rollingmean,color="ACLED: 7-day average")) + ggtitle("daily burials and deaths in armed conflicts/terror attacks") +
  geom_line(aes(y=new_graves_best_ipol,color="daily burials")) + scale_color_manual(values=c("red","pink","black")) + 
  scale_x_date(date_breaks="month",expand=expansion(0.001,0)) + scale_y_continuous(breaks=5*(0:20),expand = expansion(0.01,0)) +
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5),legend.position = "bottom")
# SAVE
ggsave("simul_output/somalia/acled_banadir_burials_comparison.png",width=30,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### 
### flight data
# flightlist_20191101_20191130 <- read_csv("data/somalia/flight_data/flightlist_20191101_20191130.csv")