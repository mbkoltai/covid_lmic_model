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
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","countrycode","coronavirus",
         "RcppRoll","scales","dttr2","wpp2019","foreach","parallel","doParallel"), library,character.only=TRUE)
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
# functions and plotting theme
source("covid_LIC_fcns.R"); source("covid_IFR_fcns.R")

### JHU global covid19 data ----------------
#' ## JHU global covid19 data
data("coronavirus")
#' data for somalia
N_tot=fun_cntr_agestr("Somalia",i_year="2020",
                      age_groups=data.frame(age_group=c(1:16),age_low=c(seq(0,75,5)),age_high=c(seq(4,74,5),100)))
covid_somal=coronavirus %>% filter(country %in% "Somalia") %>% mutate(rollingmean=roll_mean(cases,7,align="center", fill=NA)) %>%
  mutate(per_million=rollingmean/(sum(N_tot)/1e6),name=str_replace(type,"confirmed","confirmed cases"))

#' ### Plot timecourse of cases/deaths
ggplot(subset(covid_somal,!type %in% "recovered"),aes(x=date)) + 
  geom_line(aes(y=rollingmean)) + # geom_point(aes(y=cases),size=0.5,fill=NA,shape=1) + #
  scale_x_date(name="", limits=c(as.Date("2020-03-01"), NA_Date_),date_breaks="weeks") + facet_wrap(~name,scales="free",nrow=2) +
  geom_rect(aes(xmin=as.Date("2020-04-20"),xmax=as.Date("2020-04-24"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_rect(aes(xmin=as.Date("2020-05-10"),xmax=as.Date("2020-05-15"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  geom_rect(aes(xmin=as.Date("2020-05-28"),xmax=as.Date("2020-06-01"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) +
  theme_bw() + standard_theme + ylab("") +
  labs(title="COVID19 in Somalia: confirmed cases and deaths (7-day rolling mean)",caption="source: JHU CCSE") 
#' No lag in deaths compared to cases
# save
# tcourse_filename="simul_output/somalia_output/reported_cases_deaths.png" # _dots
# ggsave(tcourse_filename,width=30,height=18,units="cm") # _dots

### Linelist data ----------------
#' ## LINELIST DATA
#' ### Testing data on national level
national_charts_testing=read_csv("somalia_linelist/national_charts_testing.csv")
selcolnames=c("Total Tested","7 day moving average","Number of test for 1 confirmed case","TOTAL TESTED (%)",
              "Confirmed cases","TESTED POSITIVE from total (%)","CUMULATIVE CASES","Death","CUMULATIVE DEATHS","Recovery",
              "CUMULATIVE RECOVERY","Positivity rate (%)","CFR (%)")
national_charts_testing = national_charts_testing %>% pivot_longer(cols=all_of(selcolnames))
# plot
ggplot(subset(national_charts_testing,!Month %in% 'TOTAL' & 
                !grepl("Recovery|CUMUL|TOTAL TESTED|from total",name)),aes(x=`Week Number`,y=value,group=1)) + 
  geom_line() + geom_point() + facet_wrap(~name,scales="free",nrow=2) + theme_bw() + standard_theme + xlab("") +
  geom_rect(aes(xmin=4,xmax=6,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + 
  labs(title="National-level COVID19 metrics, weeks 12-32",caption="source: line list")
#' Testing, detected cases and deaths went up together in weeks 15-16
# save
# ggsave("simul_output/somalia_output/national_charts_testing.png",width=32,height=18,units="cm")

#' ### By state
state_charts=read_csv("somalia_linelist/state_charts.csv")
state_charts = state_charts %>% pivot_longer(cols = !c("Month","Week Number","State"))
state_charts$`Week Number`=gsub("Week-","w",state_charts$`Week Number`)
state_charts[,"week"]=as.numeric(gsub("w","",state_charts$`Week Number`))
# variables are
# c("TOTAL TESTED", "Weekly Average", "Confirmed Case", "CUMULATIVE CASES", "Death", "CUMULATIVE DEATHES",
# "CFR (%)", "Recovery", "CUMULATIVE RECOVERY", "Positivity (%)")

# plot
# ggplot(subset(state_charts, !grepl("Recovery|CUMUL|from total|Weekly",name) & week<33),
#        aes(x=`Week Number`,y=value,group=1)) + geom_line() + geom_point(size=0.5,fill=NA,shape=1) +
#   facet_grid(name~State,scales = "free",switch="y") + theme_bw() + standard_theme + theme(axis.text.x = element_text(size=6))
# # save
# ggsave("simul_output/somalia_output/state_charts.png",width=32,height=18,units="cm")

# plot with facet_wrap
ggplot(subset(state_charts, !grepl("Recovery|CUMUL|from total|Weekly",name) & week<33),
       aes(x=`Week Number`,y=value,group=1)) + geom_line() + geom_point(size=0.5,fill=NA,shape=1) +
  facet_wrap(name~State,scales="free",ncol=length(unique(state_charts$State))) +
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=7)) + 
  labs(title="COVID19 metrics by state, weeks 12-32",caption = "source: line list")
#' most testing in Banadir and Somaliland
# save
# ggsave("simul_output/somalia_output/state_charts_facetwrap.png",width=30,height=20,units="cm")

#' ### Attack rate by district
attackrate_by_district=read_csv("somalia_linelist/attackrate_by_district.csv")
selcolnames=c("Total Pop","Positivity","Attack rate per 100K","Number of Cases","Total Tested",
              "Increase or decrease from mean","tests_per_100k","tests_per_case")
attackrate_by_district=attackrate_by_district %>% 
  mutate(tests_per_100k=(`Total Tested`/`Total Pop`)*1e5,tests_per_case=`Total Tested`/`Number of Cases`) %>% 
  pivot_longer(cols=all_of(selcolnames)) %>% 
  mutate(jointname=ifelse(name %in% c("Attack rate per 100K pop","tests_per_100k"),"attackrate_tests_per100k",name))
# plot
p <- ggplot(subset(attackrate_by_district,!grepl("Increase|Total Pop|Number of|Total Tested|tests_per_case",name)),
            aes(x=`DISRICT CODE`,y=value,group=1)) + geom_line() + geom_point(fill=NA,shape=1) + 
  facet_wrap(name~STATE,scales="free",ncol=7) + ylab("") + # facet_grid(name~STATE,scales="free",switch="y") + 
  theme_bw() + standard_theme + theme(axis.text.x = element_text(size=7),legend.position="top",legend.title=element_blank()) +
  labs(title="Cumulative attack rate, positivity and testing rate by state",caption = "source: line list")
p
#' Most data is from Banadir, some from Somaliland. Only in Banadir have most districts a non-zero/NA value.
#
# if (class(p$facet)[1] %in% "FacetWrap") {filenme="attackrate_by_district_facetwrap.png"} else {filenme="attackrate_by_district.png"}
# ggsave(paste0("simul_output/somalia_output/",filenme),width=32,height=16,units="cm")

#' ### Age distribution
# "summary_by_age_grp.csv" | "cases_by_age_gender.csv" | "death_by_age_gender.csv"
summary_by_age_grp=read_csv("somalia_linelist/summary_by_age_grp.csv")
# death_by_age_gender=read_csv("somalia_linelist/death_by_age_gender.csv")
somal_agestr=data.frame(`Age Group`=popF[popF$name %in% "Somalia","age"],
            value=popF[popF$name %in% "Somalia","2020"]+popM[popM$name %in% "Somalia","2020"])
summary_by_age_grp$total_popul=c(c(sapply(seq(1,19,2),function(x) {sum(somal_agestr$value[c(x,x+1)])})*1e3,
                                   somal_agestr$value[nrow(somal_agestr)]*1e3),
      sum(c(sapply(seq(1,19,2),function(x) {sum(somal_agestr$value[c(x,x+1)])})*1e3,
            somal_agestr$value[nrow(somal_agestr)]*1e3)));
summary_by_age_grp$total_popul_perc=c(N_tot,sum(N_tot))/sum(N_tot)

summary_by_age_grp=summary_by_age_grp %>% mutate(total_tested_perc=`Total Tested`/`Total Tested`[`Age Group` %in% "National"]) %>% 
  pivot_longer(cols=!c("Age Group"))
summary_by_age_grp$`Age Group`=factor(summary_by_age_grp$`Age Group`,levels=unique(summary_by_age_grp$`Age Group`))

# plot
ggplot(subset(summary_by_age_grp,!`Age Group` %in% "National" & !name %in% "total_popul"),aes(x=`Age Group`,y=value)) + 
  geom_bar(stat="identity",color="black",width =0.7) + facet_wrap(~name,scales="free") + theme_bw() + standard_theme + 
  labs(title="Age distribution of cumulative COVID19 metrics (weeks 12-32)",caption = "source: line list") + xlab("") + ylab("")
# SAVE
# ggsave("simul_output/somalia_output/summary_by_age_grp.png",width=32,height=24,units="cm")

# rmarkdown::render("somalia_data_model.R",output_dir = "simul_output/somalia_output/")

### SATELLITE DATA on cemeteries ----------------
#' satellite image data
# image is "somalia_linelist/out_overall_weekly_trends_area.png"
# reading data from the figure by eye
# grey curve (new surface area): 
# 07/2019 to 11/2019: ~580 (jump in mid-november)
# 12/2019-01/2020: 650
# 15/01 to 01/02: jump to ~1000
# 01/03 to 01/04: rising to ~1300, until 01/05
# 01/05: drops to 900
# 01/06 to 01/07: drops to 600
# to 01/09: stays at 600
# import data from Francesco
out_bdr_weekly_new_area=read_csv("somalia_linelist/Mogadishu_data/out_bdr_weekly_new_area.csv")
out_bdr_monthly_new_area=read_csv("somalia_linelist/Mogadishu_data/out_bdr_monthly_new_area.csv")
 
newsurface_weekly=subset(out_bdr_weekly_new_area,date>"2019-10-14" & date<max(date)); newsurface_weekly=newsurface_weekly[order(newsurface_weekly$date),]
ggplot(newsurface_weekly,aes(x=date,y=new_area_ipol)) + geom_line() + geom_point() + theme_bw() + standard_theme + 
  scale_x_date(date_breaks="4 weeks")
# manual approx from plot:
# data.frame(date=seq(as.Date("2019-09-01"),as.Date("2020-08-01"),14),
#      new_area_ipol=c(rep(580,5),rep(650,6),rep(1000,3),1150,rep(1300,3),rep(900,3),750,rep(600,2)),datasource="satellite")
newsurface_weekly[,"value_norm"] = (newsurface_weekly$new_area_ipol-min(newsurface_weekly$new_area_ipol))/(max(newsurface_weekly$new_area_ipol)-min(newsurface_weekly$new_area_ipol))
newsurface_weekly = newsurface_weekly %>% mutate(value_norm_rollmean=roll_mean(value_norm,5,align="center", fill=NA),datasource="satellite")
# cemetery newsurface weekly
ggplot(newsurface_weekly %>% pivot_longer(cols=contains("value")),aes(x=date,y=value,group=name,color=name)) + geom_line() + 
  theme_bw() + standard_theme + scale_x_date(date_breaks="2 weeks") + ggtitle("new surface area index of Mogadishu cemeteries (m2)")
# SAVE
# ggsave("simul_output/somalia_output/satellite_newsurface_weekly.png",units="cm",height=15,width=25)

# compare to reported deaths
weekly_deaths_sel=data.frame((subset(covid_somal,name %in% "death") %>% 
  mutate(sum_weekly=roll_sum(cases,7,fill=NA,align="right")))[,c("date","sum_weekly")],datasource="reported"); colnames(weekly_deaths_sel)[2]="value"
weekly_deaths_sel[,"new_area_ipol"]=(max(newsurface_weekly$new_area_ipol)-min(newsurface_weekly$new_area_ipol))*(
    weekly_deaths_sel$value/max(weekly_deaths_sel$value,na.rm=T) ) + min(newsurface_weekly$new_area_ipol)
weekly_deaths_sel[,"value_norm"]=(weekly_deaths_sel$value-min(weekly_deaths_sel$value,na.rm=T))/(max(weekly_deaths_sel$value,na.rm=T)-min(weekly_deaths_sel$value,na.rm=T))
weekly_deaths_sel = weekly_deaths_sel %>% mutate(value_norm_rollmean=roll_mean(value_norm,7,align="center", fill=NA))
# plot together
selcols=c("date","value_norm","value_norm_rollmean","datasource")
ggplot(bind_rows(weekly_deaths_sel[,selcols],newsurface_weekly[,selcols]) %>% pivot_longer(cols=!c(date,datasource)),
  aes(x=date,y=value,group=datasource,color=datasource)) + geom_line() + geom_point(aes(shape=datasource)) + facet_wrap(~name) +
  theme_bw() + standard_theme + scale_x_date(date_breaks="2 weeks") + 
  ggtitle(expression(paste("new surface area index (m2) - reported deaths (scaled)")))
# SAVE
# ggsave("simul_output/somalia_output/satellite_reported_deaths.png",units="cm",height=15,width=25)

# Rt estimate: https://epiforecasts.io/covid/posts/national/somalia/

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### simple SEIR model for Somalia satellite data ----------------------
# ODEs in fcn file

# Oxford Stringency Index
# Gov stringency index
OxCGRT_somalia=subset(read_csv("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"),CountryName %in% "Somalia")
OxCGRT_somalia[,"date"]=as.Date(paste(sapply(strsplit(as.character(OxCGRT_somalia$Date), "(?<=[0-9]{4})", perl=TRUE), "[[",1),
          sapply(strsplit(sapply(strsplit(as.character(OxCGRT_somalia$Date), "(?<=[0-9]{4})", perl=TRUE), "[[",2), 
          "(?<=[0-9]{2})", perl=TRUE),function(x){paste(x,collapse="-")}),sep="-"))
OxCGRT_somalia[,"NPI_on"]=0; OxCGRT_somalia$NPI_on[min(which(OxCGRT_somalia$StringencyIndex>0)):nrow(OxCGRT_somalia)]=1
# timespan of model
OxCGRT_input=merge(data.frame(date=seq(as.Date("2019-11-01"),max(OxCGRT_somalia$date),1)),
                   OxCGRT_somalia[,c("date","StringencyIndex","NPI_on")],by="date",all=TRUE)
OxCGRT_input$StringencyIndex[1:which.min(is.na(OxCGRT_input$StringencyIndex))-1]=0
OxCGRT_input=OxCGRT_input[!is.na(OxCGRT_input$StringencyIndex),]; OxCGRT_input$NPI_on[is.na(OxCGRT_input$NPI_on)]=0
# need to convert it into [0,1] to scale susceptibility (assume: pre-data period had not restrictions)
OxCGRT_input[,"OxCGRT_scaled"]=1-(OxCGRT_input$StringencyIndex)/100; # timespan_dates=OxCGRT_input$date

# var names
sir_varnames=c("S","E","I_R","I_D","R","D","newinfvar")
# duration of simul
timesteps<-seq(0,dim(OxCGRT_input)[1],by=0.25)
# somalia population
somal_popul_tot=sum(N_tot); mogadishu_popul=2.2e6
# transmission parameters
# d_e (pre-infectious days) | d_p (presympt inf.) | d_c (Duration of symptomatic infectiousness in days) | d_death (timescale of deaths)
d_e=mean(rgamma(1e4,shape=4,scale=4/4)); d_p=mean(rgamma(1e4,shape=4,scale=1.5/4)); d_c=mean(rgamma(n=1e4,shape=4,scale=3.5/4))
d_death=1/21; sigma=1/d_e; gamma=1/d_c
# TRANSMISSION parameter (same as susceptibility) (without popul normalisation, done in the fcn)
beta_val=0.5; R0=beta_val/gamma
# IFR
cntr_ex="Somalia"; pop_struct=popF[popF$name %in% cntr_ex,"2020"] + popM[popM$name %in% cntr_ex,"2020"]
IFR_estim=0.5/100 # IFR_estim=c((ifr_all_age_groups$mean %*% pop_struct/sum(pop_struct) ))
# scan values for date of introd
introd_date_scanvals=seq(as.Date("2019-11-21"),as.Date("2020-02-01"),16)
# scan values for seed size
seedsize_scanvals=seq(2,62,15); seeding_duration=10
# compliance with NPIs
compliance_val=0.5
var_categ_list=list(name_vars=c("t", "S", "E","I_R","I_D","R","D","newinfvar","beta"),
    sel_vars=c("t","S","I_R","R","D","new_infections","new_deaths","symptom_cases"), # "I_D",
    case_vars=c("S","I_R","R","symptom_cases","new_infections"),cumul_var=c("S","R","D"), 
    delta_var=c("new_infections","new_deaths"))
# param table of all permutations
scan_params=list(seedsize_scanvals=seedsize_scanvals,introd_date_scanvals=introd_date_scanvals,seeding_duration=seq(7,35,7),
                 betaval=seq(0.35,0.55,0.05),compliance=seq(0.15,0.8,0.15),IFR=seq(0.1,0.9,0.2))
param_table=expand.grid(scan_params); for (k in 1:length(scan_params)) {param_table[,k]=scan_params[[k]]}; colnames(param_table)=names(scan_params)
# SET UP CORES for parallel
cores=detectCores(); cl<-makeCluster(cores[1]-2); registerDoParallel(cl)
# RUN
ptm<-proc.time()
df_ode_solution_scan <- foreach(k=1:nrow(param_table),.combine=rbind,.packages=c("tidyr","deSolve","dplyr")) %dopar% {
# for (k in 1:nrow(param_table[1:100,])) {
# init conds
initvals_S_I=c(mogadishu_popul,rep(0,length(sir_varnames)-1),beta_val); names(initvals_S_I)=c(sir_varnames,"beta_dyn")
# RUN & process output
num_params=c(beta=param_table$betaval[k],gamma=gamma,d_death=d_death,IFR=param_table$IFR[k])
df_ode_solution=fcn_somal_sir_singlesimul(sir_varnames,var_categ_list,timesteps,day0=param_table$introd_date_scanvals[k],
  OxCGRT_input$date,num_params,y_val,seed_size_val=param_table$seedsize_scanvals[k],seeding_duration=param_table$seeding_duration[k],
  OxCGRT_input,compliance=param_table$compliance[k],N_tot*mogadishu_popul/sum(N_tot),c(),"")[[1]] %>% filter(name %in% "new_deaths")
# df_ode_solution_scan[[k]]=df_ode_solution 
# if (k %% 10 == 0) {print(round(k/100,3))}
}
time_sim=proc.time()-ptm; time_sim # df_ode_solution_scan=do.call(rbind,df_ode_solution_scan)
stopCluster(cl)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### plot one simulation -----
# timecourse for one param set
beta_val=0.4; num_params=c(beta=beta_val,gamma=gamma,d_death=d_death,IFR=IFR_estim); introd_date=introd_date_scanvals[5]
singlesimul <- fcn_somal_sir_singlesimul(sir_varnames,var_categ_list,timesteps,day0=introd_date,OxCGRT_input$date,
      num_params,y_val,seed_size_val=20,seeding_duration,OxCGRT_input,compliance=0.5,
      N_tot*mogadishu_popul/sum(N_tot),xlimvals=c(30,400),"plot"); singlesimul[[2]]
# SAVE
# filetag_tcourse=paste0(c(paste0(names(num_params)[c(1,4)],"_",round(num_params,3)[c(1,4)]*c(1,1e2)),"compliance",compliance_val,"day",introd_date),collapse="_")
# ggsave(paste0("simul_output/somalia_output/SIR_timecourse_",filetag,".png"),width=25,height=20,units="cm")
# PEAK cases and deaths
peak_dates=subset(singlesimul[[1]],name %in% c("new_infections","new_deaths")) %>% group_by(name) %>% summarise(maxval_date=date[value==max(value)])
ggplot(subset(singlesimul[[1]], name %in% c("new_infections","new_deaths")),aes(x=date,y=value,group=name,color=name)) + 
  geom_line() + geom_point() + scale_x_date(limits=as.Date(c("2020-01-01","2020-08-01")),date_breaks="1 week") + scale_y_log10(limits=c(2,2e4)) +
  geom_vline(data=peak_dates,aes(xintercept=maxval_date,color=name)) + theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5))
# ggsave(paste0("simul_output/somalia_output/peak_cases_deaths",filetag,".png"),width=25,height=15,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

### discretise number of new deaths (for parscan) -----
# ESTIMATE of death rate for mogadishu: 0.4 deaths/10K ppl/day -> Mogadishu is 2.2M -> 88 deaths/day -> ~1250/2 weeks
crude_deathrate_per1e4personday=0.4; mogad_weekly_death_rate=(mogadishu_popul/1e4)*crude_deathrate_per1e4personday*7
newdeaths_pred=subset(df_ode_solution_scan,name %in% c("new_deaths")) %>% 
  mutate(sum_weekly=roll_sum(value,7,fill=NA,align="right")) %>% group_by(introd_date,seed_size) %>%
  mutate(value_norm=sum_weekly/mogad_weekly_death_rate, datasource="model (excess deaths %)")
# concatenate data types + simul
pred_data_comp=left_join(left_join(subset(newdeaths_pred, date %in% newsurface_weekly$date),newsurface_weekly,by="date", 
                         suffix=c("","_satel_data")),weekly_deaths_sel,by="date",suffix=c("","_report_data")) 
#[,colnames(newsurface_weekly)]
# calculate ERROR
mse_pred_data_comp=pred_data_comp %>% group_by(introd_date,seed_size) %>% 
  summarise(mse_satel=mean((value_norm - value_norm_satel_data)^2,na.rm=T), 
            mse_reported=mean((value_norm - value_norm_report_data)^2,na.rm=T)) %>% pivot_longer(cols=!c(introd_date,seed_size))

### PLOT ALL trajectories ------------------------
ggplot(newdeaths_pred) + # subset(newdeaths_pred, date %in% newsurface_weekly$date)
  geom_line(aes(x=date,y=value_norm,color=datasource)) + # ,size=1.2 geom_bar(stat="identity",fill=NA) + # model
  geom_line(data=newsurface_weekly,aes(x=date,y=value_norm,color=datasource)) + # ,size=1.05 satellite
  # geom_line(data=weekly_deaths_sel[,c("date","value_norm","datasource")],aes(x=date,y=value_norm,color=datasource),size=1.05) +
  facet_grid(introd_date~seed_size,labeller=labeller(seed_size=label_both),scales="free") + 
  theme_bw() + standard_theme + theme(legend.position="top",axis.text.x=element_text(size=5),axis.text.y=element_text(size=5)) +
  geom_vline(aes(xintercept=introd_date),color="red") + geom_vline(xintercept=as.Date(min(OxCGRT_input$date[OxCGRT_input$NPI_on>0])),color="green") +
 scale_x_date(date_breaks="2 weeks",limits=c(min(OxCGRT_input$date)+15,max(newsurface_weekly$date)+30)) + xlab('time (days)') + ylab('number') +
 labs(title="SIR model ~ f(introduction date & seed size)",
       color="",subtitle=paste0("beta=",beta_val,", compliance (NPI)=",compliance_val,", IFR=",IFR_estim*1e2,"%"))
# SAVE
# filetag=paste0(c(paste0(names(num_params)[c(1,4)],"_",round(num_params,3)[c(1,4)]*c(1,1e2)),"compliance",compliance_val),collapse="_")
# ggsave(paste0("simul_output/somalia_output/param_scan_satellite_SIRmodel_",filetag,".png"),width=40,height=24,units="cm")

### ### ### ### ### ### ### ### ### ### ### ###
### heatmap of MSE ------------------------
minval_mse=min(subset(mse_pred_data_comp,grepl("satel",name))[,"value"])
highl_tiles=subset(mse_pred_data_comp,grepl("satel",name) & value<minval_mse*1.1) %>% select(introd_date,seed_size)
midpoint_val=median(subset(mse_pred_data_comp,grepl("satel",name))$value)
ggplot(subset(mse_pred_data_comp,grepl("satel",name)),aes(x=seed_size,y=introd_date,fill=value)) + 
 geom_tile(color="black") + geom_text(aes(label=round(value,3)),color="black",size=4) + 
 geom_rect(data=highl_tiles,size=1,fill=NA,colour="black", aes(xmin=seed_size-unique(diff(seedsize_scanvals))/2,xmax=seed_size+unique(diff(seedsize_scanvals))/2,
                    ymin=introd_date-unique(diff(introd_date_scanvals))/2,ymax=introd_date+unique(diff(introd_date_scanvals))/2)) +
 theme_bw() + standard_theme + scale_fill_gradient2(midpoint=midpoint_val,low="blue",mid="white",high="red") + 
 scale_y_date(breaks=unique(mse_pred_data_comp$introd_date),expand=c(0,0)) + 
 scale_x_continuous(breaks=unique(mse_pred_data_comp$seed_size),expand=c(0,0)) +
 labs(title="MSE for grid search",subtitle=paste0("beta=",beta_val,", compliance (NPI)=",compliance_val,", IFR=",IFR_estim*100,"%")) + 
 xlab("seed size (for 10 days)") + ylab("introd. date")
# SAVE
# filetag=paste0(c(paste0(names(num_params)[c(1,4)],"_",round(num_params,3)[c(1,4)]*c(1,1e2)),"compliance",compliance_val),collapse="_")
# ggsave(paste0("simul_output/somalia_output/mse_param_scan_satellite_SIR_heatmap_",filetag,".png"),width=32,height=20,units="cm")
### ### ### ### ### ### ### ### ### ### ### ###
### lineplot of MSE ------------------------
ggplot(subset(mse_pred_data_comp, grepl("satel",name)),aes(x=introd_date,y=value,group=seed_size,color=factor(seed_size))) + 
  geom_line(size=1.05) + geom_point() + theme_bw() + standard_theme + scale_x_date(breaks=unique(mse_pred_data_comp$introd_date)) + 
  scale_y_log10() + xlab("introduction date") + ylab("MSE") + labs(color="seed size",title="MSE for grid search",
      subtitle=paste0("beta=",beta_val,", compliance (NPI)=",compliance_val,", IFR=",IFR_estim*100,"%"))
# ggsave(paste0("simul_output/somalia_output/mse_param_scan_satellite_SIR_lineplot_",filetag,".png"),width=32,height=20,units="cm")


### MCMC fitting --------------
# remotes::install_github('sbfnk/fitR') # package from S Funk for MCMC
# library('fitR')
