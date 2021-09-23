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
         "scales","wpp2019","GGally","corrr","ungeviz"), library,character.only=TRUE)
# detach("package:fitdistrplus", unload = TRUE); detach("package:MASS", unload = TRUE) # "foreach","parallel","doParallel"
# functions and plotting theme
source("somalia_data_model_fcns.R")

### JHU global covid19 data ----------------
#' ## JHU global covid19 data
data("coronavirus")
# age structure
N_tot=fun_cntr_agestr("Somalia",i_year="2020",age_groups=data.frame(age_group=c(1:16),age_low=c(seq(0,75,5)),
                                                                    age_high=c(seq(4,74,5),100)))
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
# satellite data vs numbers by Barakaat Cemetery Developm Committee
satellite_bcdc_registry <- read_csv("data/satellite_bcdc_registry.csv")
ggplot(satellite_bcdc_registry) + 
  geom_point(aes(x=`Satellite imagery`,y=`Barakaat Development Committee`,shape=factor(Year),color=factor(Month)),size=4) +
  theme_bw() + standard_theme + ylab("") + scale_y_continuous(expand=expansion(0.02,0)) + # scale_fill_discrete() +
  labs(title="satellite data compared to cemetery registry")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### SATELLITE DATA on cemeteries ----------------
#' satellite image data
burial_data=read_csv("data/somalia_data/Mogadishu_data/mogadishu_burial_analysis-main/out_bdr_daily_estimates.csv")
# pop_wp2015,pop_wp2020,new_graves_best_ipol,
ggplot(burial_data %>% select(date,br_wp2015,br_wp2015_base_s,br_wp2020,br_wp2020_base_s) %>% pivot_longer(!date)) + 
  geom_line(aes(x=date,y=value,group=name,color=name)) + # facet_wrap(~name,scales="free") + 
  scale_x_date(expand=expansion(0.01,0),date_breaks="2 month") + theme_bw() + standard_theme

# baseline of daily burials july-november 2019
baseline_daily_burials <- mean(subset(burial_data,date>="2019-07-01" & date<="2019-11-01")$new_graves_best_ipol)
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
ggplot(subset(out_bdr_daily_estimates %>% pivot_longer(col=c(new_graves_best_ipol,rollmeanweek,rollsumweek)),grepl("roll",name)), # 
          aes(x=date,y=value)) + facet_wrap(~name,scales="free") + geom_line() + geom_point(size=0.3) + theme_bw() + standard_theme +
  scale_x_date(date_breaks="2 weeks",expand=expansion(0.0)) + theme(axis.text.x=element_text(vjust=0.5))
# compare # of burials, daily
ggplot(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr), # 
    aes(x=date,y=daily_baseline_subtr)) + geom_line() + geom_point(size=0.3) + theme_bw() + standard_theme +
  scale_x_date(date_breaks="2 weeks",expand=expansion(0.0)) + theme(axis.text.x=element_text(vjust=0.5))

# compare to reported deaths
weekly_deaths_reported=data.frame(subset(covid_somal,name %in% "death")[,c("date","value")],datasource="reported") %>% 
  mutate(datasource=as.character(datasource)) %>% rename(value_report_daily=value) %>%
  mutate(rollmeanweek=roll_mean(value_report_daily,7,align="center",fill=NA),
         rollsumweek=roll_sum(value_report_daily,7,fill=NA,align="right"))
# plot together
df_compare_report_satell <- bind_rows(weekly_deaths_reported,
  out_bdr_daily_estimates %>% select(date,daily_baseline_subtr,rollmeanweek,rollsumweek) %>% mutate(datasource="satellite") ) %>% 
  pivot_longer(cols=!c(date,datasource)) %>% filter(!is.na(value)) %>% #  & name=="rollsumweek"
  mutate(type=ifelse(!grepl("sum|mean",name),"daily",name),week=format(date,"%Y/%W")) %>%
  filter(type=="daily") %>% group_by(week,datasource) %>% 
  summarise(date=min(date),datasource=unique(datasource),name=unique(name),value=sum(value)) %>% mutate(week=gsub("/","/w",week))
# plot # library(ungeviz)
p <- ggplot(subset(df_compare_report_satell,date>"2020-01-15" & date<"2020-10-07") %>% mutate(week=gsub("2020/","",week)),
            aes(x=week,y=value,group=datasource)) + geom_hpline(aes(x=week,y=value,group=datasource,color=datasource),width=0.9) + # 
  geom_vline(xintercept=(1:length(unique(df_compare_report_satell$week)))-0.5,size=0.1,linetype="dashed") + 
  # geom_bar(aes(fill=datasource),stat="identity",position=position_dodge(width=0.75),color="black",size=0.2) + 
  scale_fill_discrete(labels=c("reported COVID19 deaths","excess burials")) + 
  scale_color_discrete(labels=c("reported COVID19 deaths","excess burials")) + labs(fill="",color="") + theme_bw() + standard_theme +
  theme(axis.text.x=element_text(vjust=0.5,size=15),axis.text.y=element_text(size=15),axis.title=element_text(size=17),
        legend.position=c(0.855,0.925),legend.text=element_text(size=15),legend.box.background=element_rect(colour="black",size=0.8),
        legend.background=element_blank(),legend.spacing.y=unit(0,"mm")) + 
  scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*5) + xlab("calendar week") + ylab("number per week"); p
# SAVE
if (any(grepl("Bar",class(p$layers[[1]]$geom)))) {plotfilename<-"satellite_burials_reported_deaths_weekly_barplot"} else {
  plotfilename<-"satellite_burials_reported_deaths_weekly" }
ggsave(paste0("simul_output/somalia/",plotfilename,".png"),units="cm",height=18,width=30)

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
            data.frame(NPIvals) %>% rownames_to_column(var="name"),by="name") %>% mutate(name=factor(name,levels=unique(name))) %>%
  rename(contact_level=NPIvals) %>% mutate(contact_reduction=1-contact_level)
# plot
ggplot(OxCGRT_input) + geom_line(aes(x=date,y=1-OxCGRT_scaled)) + geom_segment(data=npi_df,aes(x=on,xend=off,y=contact_reduction,
      yend=contact_reduction,group=factor(name),color=factor(name)),size=2) + theme_bw() + standard_theme + labs(color="NPI phases") + 
  theme(axis.text.x=element_text(vjust=0.5,size=14),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16)) +
  scale_x_date(limits=c(min(subset(OxCGRT_input,NPI_on>0)$date)-15,max(out_bdr_daily_estimates$date)+35),breaks="week",
               expand=expansion(0.01,0)) + scale_y_continuous(breaks=(0:10)/10) + xlab("") + ylab("stringency of NPIs (=StringencyIndex)") +
  geom_vline(data=npi_df,aes(xintercept=on,color=name),linetype="dashed")
# SAVE
# ggsave(paste0("simul_output/somalia/OxCGRT_input_contactreduction.png"),width=30,height=18,units="cm")

### Somalia population, IFR ------------
# IFR estimates from Imperial
somalia_agegroups_IFR<-fcn_merge_ifr_above_age(left_join(fcn_load_age_str("Somalia",n_year="2015",90),
  fcn_load_ifr("repo_data/IFR_by_age_imperial.csv"),by=c("agegroup","agegroup_min")),75) %>% 
  mutate(ifr_mean=ifelse(ifr_mean==0,min(ifr_mean[ifr_mean>0]),ifr_mean),log_ifr=log(ifr_mean),logit_ifr=log(ifr_mean/(1-ifr_mean)))
somal_popul_tot=sum(somalia_agegroups_IFR$agegroupsize); mogadishu_popul=2.2e6 # somalia_agegroups_IFR$ifr_mean[1]=3e-6; 
# other IFR estimates
# from Sandmann 2021 cmmid paper
IFR_estimates_Sandmann2021<-read_csv("repo_data/IFR_estimates_Sandmann2021.csv")
 if (any(IFR_estimates_Sandmann2021$value_percent>1)) {n_cols<-2:ncol(IFR_estimates_Sandmann2021)
   IFR_estimates_Sandmann2021[,n_cols]<-IFR_estimates_Sandmann2021[,n_cols]/1e2; 
   IFR_estimates_Sandmann2021 <- left_join(IFR_estimates_Sandmann2021 %>% rename(agegroup=Age,ifr_mean=value_percent), 
      somalia_agegroups_IFR %>% select(!c(ifr_mean,log_ifr,logit_ifr)),by="agegroup") %>% mutate(logit_ifr=log(ifr_mean/(1-ifr_mean)))
   }
ggplot(data.frame(age=factor(somalia_agegroups_IFR$agegroup,levels=unique(somalia_agegroups_IFR$agegroup)),
                  imper=somalia_agegroups_IFR$logit_ifr,sandmann=IFR_estimates_Sandmann2021$logit_ifr) %>% pivot_longer(!age),
  aes(x=age,y=value,group=name,color=name)) + geom_line() + geom_point() + theme_bw() + ylab("logit(IFR)")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# COVIDM
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
cm_force_rebuild=F; cm_build_verbose=T; cm_version=1; source(file.path(cm_path,"R","covidm.R"))
countryval="Somalia"; params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval),
                                                 date_start="2019-11-01",date_end="2020-10-01")
# set population: Somalia --> Mogadishu
params$pop[[1]]$name=countryval
params$pop[[1]]$size=somalia_agegroups_IFR$agegroupsize*(mogadishu_popul/sum(somalia_agegroups_IFR$agegroupsize))
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### SEEDING ---
npi_on_day=min(OxCGRT_input$date[OxCGRT_input$NPI_on>0]); introd_day=as.Date("2019-11-21")
seeding_t_window=sapply(rep(introd_day,2),function(x) as.numeric(x-as.Date(params$date0))); seed_size=50
params$pop[[1]]$seed_times=rep(seeding_t_window[1]:seeding_t_window[2],each=seed_size) # x new infections/day for n days
# infections start in individuals aged 20-50, 1 introd in each age group
params$pop[[1]]$dist_seed_ages=cm_age_coefficients(20,60,5*(0:length(params$pop[[1]]$size)))
# set approximate clin fract values (From Davies 2020 Nat Med)
params$pop[[1]]$y=fun_lin_approx_agedep_par(agegroups=somalia_agegroups_IFR,min_val=0.25,max_val=0.7,rep_min=6,rep_max=2)
# suscept and clinical fraction age dependent
suscept_clinfract_posteriors<-read_csv("data/suscept_clinfract_posteriors_davies2010.csv") %>% 
  mutate(agegroup=factor(agegroup,levels=unique(agegroup)))
susc_clinfract_plot=rbind(suscept_clinfract_posteriors %>% mutate(type="literature"),
    rbind(data.frame(name="clin_fract",agegroup=unique(suscept_clinfract_posteriors$agegroup),value=params$pop[[1]]$y,type="approximation"),
  data.frame(name="susceptibility",agegroup=unique(suscept_clinfract_posteriors$agegroup),value=c(rep(0.38,4),rep(0.8,12)),type="approximation"))) %>%
  mutate(name=gsub("clin_fract","clinical fraction",name))
ggplot(susc_clinfract_plot,aes(x=agegroup,y=value,color=type,group=type)) + geom_line() + geom_point() + facet_wrap(~name,scales="free") + 
  theme_bw() + standard_theme + scale_y_continuous(breaks=(2:20)/20) + theme(axis.text.x=element_text(vjust=0.5,size=14),
    axis.text.y=element_text(size=14),axis.title.x=element_text(size=16),axis.title.y=element_text(size=16),strip.text=element_text(size=16)) + 
  labs(color="") + xlab("age group (yr)") + ylab("proportion/normalised value")
# SAVE
# ggsave("simul_output/somalia/clinfract_susc_lit_approx.png",width=30,height=16,units="cm")
# plot IFR original vs adjusted
ggplot(cbind(data.frame(age=factor(somalia_agegroups_IFR$agegroup,levels=unique(somalia_agegroups_IFR$agegroup)),
                  literature=IFR_estimates_Sandmann2021$ifr_mean),
  data.frame(sapply(c(1,2,3), function(x) {inv.logit(logit(IFR_estimates_Sandmann2021$ifr_mean) + x)}))) %>% 
    rename(`logit(IFR)+1`=X1,`logit(IFR)+2`=X2,`logit(IFR)+3`=X3) %>% pivot_longer(!age),
  aes(x=age,y=value,group=name,color=name)) + geom_line() + geom_point() + theme_bw() + ylab("IFR") + labs(color="") + 
  scale_y_log10( breaks = scales::trans_breaks("log10", function(x) 10^x),
                 labels = scales::trans_format("log10", scales::math_format(10^.x)) ) + annotation_logticks(sides="lr") +
  theme(axis.text=element_text(size=13),axis.text.x=element_text(angle=90,vjust=0.5),axis.title=element_text(size=15),
        legend.text=element_text(size=15))
# SAVE
ggsave("simul_output/somalia/IFR_shifted_logit.png",width=25,height=16,units="cm")
### ### ### ### ### ### ### ### ### ###
# change susceptibility to get R0=x
target_R0=1.5; params$pop[[1]]$u=c(rep(0.38,4),rep(0.8,12)) # susceptibility_warvick_model$value # 
params$pop[[1]]$u=params$pop[[1]]$u*(target_R0/cm_calc_R0(params,1))
### add death process to model ------
ifr_logit <- c(-10.414283,-11.512915,-11.512915,-10.414283,-9.721106,-8.947846,-8.334632,-7.823646,-7.194687,
               -6.715924,-6.178135,-5.732038,-5.385862,-4.522041,-4.073073,-2.026972)
# sum(inv.logit(ifr_logit + 2)*somalia_agegroups_IFR$agegroupsize/sum(somalia_agegroups_IFR$agegroupsize))*1e2
# sum(IFR_estimates_Sandmann2021$ifr_mean*somalia_agegroups_IFR$agegroupsize/sum(somalia_agegroups_IFR$agegroupsize))*100
time_inf_to_death<-15 # IFR_estimates_Sandmann2021$ifr_mean
params$processes <- list(cm_multinom_process("Ip",outcomes=data.table(death=inv.logit(ifr_logit+2.5)/params$pop[[1]]$y),
                               delays=data.table(death=cm_delay_gamma(time_inf_to_death,time_inf_to_death,60,1/4)$p), report="o"))
# R0: cm_calc_R0(params,1)
# df for seeding
seeding_df=data.frame(seed_date=unique(covidm_simul$date)[unique(params$pop[[1]]$seed_times)]) %>%
  summarise(min=min(seed_date),max=max(seed_date))
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
    dynamics_type=list(cumul=c("D","R","S"),incid=c("cases","subclinical","death_o"),preval=c("Ia","Ip","Is")),populval=mogadishu_popul,params)
# PLOT
ggplot(subset(covidm_simul,!dynam_type %in% "preval")) + geom_area(aes(x=date,y=value,fill=compartment),color="black",size=0.3) +
  facet_wrap(dynam_type~compartm_type,scales="free") + theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) +
  scale_x_date(limits=as.Date(c("2019-12-01",params$time1)),date_breaks="2 weeks",expand=expansion(0,0)) + ylab("number") +
  scale_y_continuous(expand=expansion(0.01,0)) + geom_vline(data=npi_df,aes(xintercept=on,color=name),size=1,linetype="dashed") +
   geom_rect(data=seeding_df,aes(xmin=min,xmax=max,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.4) + labs(fill="")
# SAVE
# ggsave(paste0("simul_output/somalia/agegroups_summed_output_introddate_noninteg.png"),width=30,height=20,units="cm")

# sanity check of case dynamics
# SUMMER across age groups
plot_sum_t_series <- left_join(run$dynamics,IFR_estimates_Sandmann2021 %>% mutate(ifr_sympt=params$processes[[1]]$prob[1,]) %>% 
                                 rename(group=agegroup) %>% select(group,ifr_mean,ifr_sympt),by="group") %>% 
  filter(compartment %in% c("E","Is","cases","death_o")) %>% group_by(group,compartment) %>% #"Is" & t<200 & t>=max(seeding_t_window)-2
  mutate(deaths_calc_lag15d_cases=ifelse(compartment %in% "cases",lag(value,order_by=t,n=15)*ifr_sympt,NA),date=t+as.Date(params$date0)) %>%
  group_by(compartment,date,t) %>% summarise(value=sum(value),deaths_calc_lag15d_cases=sum(deaths_calc_lag15d_cases)) %>% ungroup() 
plot_sum_t_series <- bind_rows(plot_sum_t_series %>% select(compartment,date,t,value),plot_sum_t_series %>% 
                                 select(compartment,date,t,deaths_calc_lag15d_cases) %>% 
            rename(value=deaths_calc_lag15d_cases) %>% mutate(compartment="deaths_calc_lag15d_cases") %>% filter(!is.na(value) )) %>%
  mutate(compartment=ifelse(compartment=="death_o","deaths (simulated)",
                            ifelse(compartment %in% "deaths_calc_lag15d_cases","deaths (% of 15-day lagged presympt. cases)",compartment)))
# AGEGROUPS summed
sel_case_date<-plot_sum_t_series %>% filter(compartment=="cases") %>% 
                  filter(date==min(date[value>2^9]) | date==min(date[value>8]) | date==min(date[value>16]) ) # )$date
sel_n_cases_deaths<-plot_sum_t_series %>% filter((date %in% sel_case_date$date[3] & compartment %in% "cases") |
                               (date %in% (sel_case_date$date[3]+10) & compartment %in% "deaths (simulated)") )
cases_8<-plot_sum_t_series %>% filter(compartment=="cases") %>% filter(date==max(date[abs(value-8)<0.5]))
ifr_pop_aver=sum(params$processes[[1]]$prob[1,]*params$pop[[1]]$y*somalia_agegroups_IFR$agegroupsize/sum(somalia_agegroups_IFR$agegroupsize))
# plot cases - deaths
title_str<-expr(paste(R[0], "=",!!target_R0,", IFR (all cases)=",!!round(ifr_pop_aver*1e2,2),"%, date of introduction: ",!!as.character(seeding_df$min+1)))
ggplot() + geom_line(data=plot_sum_t_series %>% filter(!compartment %in% "Ip") %>% mutate(value=ifelse(compartment=="E"&date<=seeding_df$min,NA,value),
        compartment=ifelse(compartment=="Is","symptomatic cases (prevalence)",compartment),
        compartment=ifelse(compartment=="cases","symptomatic cases (incidence)",ifelse(compartment=="E","E (prevalence)",compartment))), 
        aes(x=date,y=value,color=compartment),size=1.1) + 
  annotate("rect",xmin=sel_case_date$date[1],xmax=sel_case_date$date[2],ymin=sel_case_date$value[1],ymax=sel_case_date$value[2],fill=NA,
           color="black",linetype="dashed",size=1/3) + annotate("rect",xmin=sel_n_cases_deaths$date[1],xmax=sel_n_cases_deaths$date[2],
           ymin=sel_n_cases_deaths$value[1],ymax=sel_n_cases_deaths$value[2],fill=NA,color="black",linetype="dashed",size=1/3) +
  geom_text(aes(x=sel_case_date$date[1]+5,y=20,label=paste0("≈",as.numeric(diff(sel_case_date$date)[1])," days")),size=4.5) +
  geom_point(data=sel_case_date %>% filter(date<as.Date("2020-03-01")),aes(x=date,y=value)) +
  geom_text(aes(x=sel_n_cases_deaths$date[1]+16,y=700,label="incident cases≈512"),size=4.5) + guides(color=guide_legend(nrow=2,byrow=TRUE)) +
  geom_text(aes(x=sel_n_cases_deaths$date[2]+18,y=sel_n_cases_deaths$value[2],label=paste0("incident deaths=",round(sel_n_cases_deaths$value[2]))),size=4.5) + 
  geom_point(aes(x=sel_n_cases_deaths$date,y=sel_n_cases_deaths$value)) + geom_point(aes(x=seeding_df$min+1,y=seed_size)) +
  geom_text(aes(x=seeding_df$min+14,y=seed_size*1.25),label=paste0(seed_size," introductions"),size=4.5) +
  geom_segment(aes(x=cases_8$date,xend=min(sel_n_cases_deaths$date),y=cases_8$value,yend=cases_8$value),arrow=arrow(length=unit(0.2,"cm"))) +
  geom_segment(aes(x=sel_n_cases_deaths$date[1],xend=sel_n_cases_deaths$date[2],y=2^7,yend=2^7),arrow=arrow(length=unit(0.2,"cm"))) +
  geom_text(aes(x=sel_n_cases_deaths$date[1]+8,y=2^7.4),label=expression(paste(t[d],"=10 days"))) + geom_point(aes(x=cases_8$date,y=cases_8$value)) + 
  geom_text(aes(x=sel_case_date$date[1]+25,y=7),label=paste0(as.numeric(sel_n_cases_deaths$date[1]-cases_8$date)," days, 6 doublings"),size=4.5) +
  scale_y_log10(limits=c(2,5e4),expand=expansion(0,0),breaks=2^(-1:15)) + theme_bw() + standard_theme+ylab("")+xlab("")+labs(color="") +
   theme(axis.text.x=element_text(vjust=0.5,size=14),axis.text.y=element_text(size=14),legend.position="bottom") +
  ggtitle(title_str) + scale_x_date(date_breaks="week",limits=c(seeding_df$min,as.Date("2020-07-06")),expand=expansion(0.01,0))
# save
ggsave(paste0("simul_output/somalia/cases_deaths_simul_ifr",round(ifr_pop_aver*1e2,2),".png"),width=30,height=20,units="cm")

# agegr separated
# ggplot(run$dynamics %>% filter(compartment %in% c("E","cases","death_o") & t<200 & t>=max(seeding_t_window)-2 & as.numeric(group)>5) %>% 
#          mutate(date=t+as.Date(params$date0))) + geom_line(aes(x=date,y=value,color=compartment)) + facet_wrap(~group,scales="free") + 
#   theme_bw() + standard_theme + scale_y_log10(limits=c(0.1,5e3)) + scale_x_date(date_breaks = "month")

# PLOT incident deaths with data
fitting_date_window=as.Date(c("2020-03-01","2020-08-04"))
fcn_covidm_singlesim_error(covidm_simul,introd_day,seedsize=50,out_bdr_daily_estimates,fitting_date_window)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ACLED data on deaths from polit violence/terrorism
somalia_acled_data <- read_csv("data/somalia/terrorism_data/2018-04-23-2021-04-28-Somalia.csv")
somalia_acled_fatalities <- somalia_acled_data %>% group_by(event_date,admin1) %>% 
  summarise(date=unique(dmy(event_date)),fatalities=sum(fatalities),year=unique(year)) %>% ungroup() %>% group_by(admin1) %>% 
  complete(date=seq.Date(min(date),max(date),by="day")) %>% mutate(fatalities_missing_as_0=ifelse(is.na(fatalities),0,fatalities),
  rollingmean=roll_mean(fatalities_missing_as_0,7,align="center",fill=NA,na.rm=T),
  rollingsum=roll_sum(fatalities_missing_as_0,7,fill=NA,align="right"),
  rollingsum_month_acled=roll_sum(fatalities_missing_as_0,30,fill=NA,align="right")) %>%
  mutate(year=year(date),month=month(date),week=week(date))
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
ggplot(subset(acled_monthly_fatalities,date<=as.Date("2020-11-01") & date>as.Date("2018-11-01")),aes(x=year_month,y=fatalities)) + 
  geom_bar(stat="identity") + facet_wrap(~admin1,scales = "free") +
  geom_rect(aes(xmin="2020/01",xmax="2020/10",ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + ylab("deaths/month") + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,size=6))
# SAVE
ggsave("simul_output/somalia/acled_somalia_monthly_fatalities.png",width=30,height=18,units="cm")

# compare to burials
acled_burial_comparison = left_join(subset(somalia_acled_fatalities,admin1=="Banadir") %>%
                    select(admin1,date,fatalities_missing_as_0,rollingmean,rollingsum,rollingsum_month_acled) %>% 
                    rename(fatalities=fatalities_missing_as_0,rollingmean_acled=rollingmean,rolling_week_sum_acled=rollingsum),
  burial_data[!rowSums(is.na(burial_data))==(ncol(burial_data)-1),!colSums(is.na(burial_data))==nrow(burial_data)] %>% 
  select(date,new_graves_best_ipol) %>% mutate(rolling_week_sum_burials=roll_sum(new_graves_best_ipol,7,fill=NA,align="right"),
  rolling_month_sum_burials=roll_sum(new_graves_best_ipol,30,fill=NA,align="right")),   by="date") 

# ggplot rolling sums
ggplot(acled_burial_comparison %>% ungroup() %>% select(date,rollingsum_month_acled,rolling_month_sum_burials) %>% pivot_longer(!date),
       aes(date)) + geom_line(aes(x=date,y=value,group=name,color=name)) + theme_bw() + standard_theme + theme(legend.position="top") + 
  labs(color="") + scale_x_date(limits=as.Date(c("2019-01-01","2020-09-25")),date_breaks="2 weeks",expand=expansion(0.01,0)) + 
  scale_y_continuous(expand=expansion(0.01,0))
ggsave("simul_output/somalia/acled_burials_rolling_monthly_sum.png",width=30,height=18,units="cm")

# weekly monthly compare (SI Figure 3)
weekly_acled_burial_comparison = acled_burial_comparison %>% mutate(week=week(date),year=year(date)) %>% group_by(year,week) %>% 
  summarise(date=min(date),fatalities=sum(fatalities),rollingmean_acled=sum(rollingmean_acled),new_graves_best_ipol=sum(new_graves_best_ipol)) %>%
  mutate(year_week=factor(paste0(year,"/",week),levels=unique(paste0(year,"/",week)))) %>% 
  mutate(year_week_merged=ifelse(week>=52,paste0(year,"/52_53"),as.character(year_week))) %>% 
  mutate(year_week_merged=factor(year_week_merged,levels=unique(year_week_merged))) %>% group_by(year_week_merged) %>% 
  summarise(year=unique(year),date=max(date),fatalities=sum(fatalities),new_graves_best_ipol=sum(new_graves_best_ipol))
# WEEKLY compare plot
p<-ggplot(subset(weekly_acled_burial_comparison,date>=as.Date("2019-01-01") & date<=as.Date("2020-10-01")) %>% pivot_longer(!c(date,year,year_week_merged))) +
#    geom_bar(aes(x=year_week_merged,y=value,group=name,fill=name),stat="identity",position=position_dodge(width=0.75),size=0.2) +
#  scale_fill_discrete(labels=c("deaths due to political violence (ACLED)","burials")) +
   geom_point(aes(x=year_week_merged,y=value,group=name,color=name),pch="-",size=12) + 
   geom_vline(xintercept=(1:length(unique(weekly_acled_burial_comparison$year_week_merged)))-0.5,size=0.1,linetype="dashed") +
   scale_color_discrete(labels=c("deaths due to political violence (ACLED)","burials")) +
   theme_bw() + standard_theme + theme(axis.text.x = element_text(vjust=0.5),legend.position="bottom") + labs(fill="",color="") +
   scale_y_continuous(expand=expansion(0.01,0),breaks=(0:20)*10) + xlab("year/week") + ylab("number/week"); p
# SAVE
# 
if (any(grepl("Point",class(p$layers[[1]]$geom)))) {plotfilename<-"acled_banadir_burials_comparison_WEEKLY_geompoint"} else {
  plotfilename<-"acled_banadir_burials_comparison_WEEKLY" }
ggsave(paste0("simul_output/somalia/ACLED_data/",plotfilename,".png"),width=30,height=16,units="cm")

# daily compare plot
ggplot(acled_burial_comparison %>% filter(date<=as.Date("2020-10-01")),aes(x=date)) + #  %>% pivot_longer(!c(admin1,date))
  geom_bar(aes(y=fatalities,color="ACLED: daily deaths (armed conflicts/terrorism)"),stat="identity",alpha=0.1) +  
  geom_line(aes(y=rollingmean_acled,color="ACLED: 7-day average")) + 
  ggtitle("daily burials and deaths in armed conflicts/terror attacks") +
  geom_line(aes(y=new_graves_best_ipol,color="daily burials")) + scale_color_manual(values=c("red","pink","black")) + 
  scale_x_date(date_breaks="month",expand=expansion(0.001,0)) + scale_y_continuous(breaks=5*(0:20),expand = expansion(0.01,0)) +
  theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5),legend.position = "bottom")
# SAVE
ggsave("simul_output/somalia/acled_banadir_burials_comparison.png",width=30,height=16,units="cm")