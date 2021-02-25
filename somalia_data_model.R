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
  geom_rect(aes(xmin=as.Date("2020-05-28"),xmax=as.Date("2020-06-01"),ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + theme_bw() + standard_theme + 
  labs(title="COVID19 in Somalia: confirmed cases and deaths (7-day rolling mean)",caption="source: JHU CCSE") + ylab("")
#' No lag in deaths compared to cases
# save
# tcourse_filename="simul_output/somalia_output/reported_cases_deaths.png" # _dots
# ggsave(tcourse_filename,width=30,height=18,units="cm") # _dots

### SATELLITE DATA on cemeteries ----------------
#' satellite image data
out_bdr_weekly_burials=read_csv("somalia_linelist/Mogadishu_data/out_bdr_weekly_burials.csv")
newburials_weekly=subset(out_bdr_weekly_burials,date>"2019-10-14" & date<max(date) & !is.na(new_graves_best_ipol))
newburials_weekly=newburials_weekly[order(newburials_weekly$date),]
ggplot(newburials_weekly,aes(x=date,y=new_graves_best_ipol)) + geom_line() + geom_point() + theme_bw() + standard_theme + 
  scale_x_date(date_breaks="2 weeks")
basal_burial_rate=min(newburials_weekly$new_graves_best_ipol) # [newburials_weekly$date<as.Date("2020-01-26")]
newburials_weekly[,"value_norm"] = (newburials_weekly$new_graves_best_ipol - basal_burial_rate)/basal_burial_rate
newburials_weekly = newburials_weekly %>% mutate(value_norm_rollmean=roll_mean(value_norm,5,align="center", fill=NA),datasource="satellite")
# cemetery newsurface weekly
ggplot(newburials_weekly %>% pivot_longer(cols=matches("value|graves")) %>% mutate(type=ifelse(grepl("norm",name),"norm","abs")),
  aes(x=date,y=value,group=name,color=name)) + geom_line() + facet_wrap(~type,scales="free") +
  theme_bw() + standard_theme + scale_x_date(date_breaks="2 weeks") + ggtitle("new burials in 6 Mogadishu cemeteries")
# SAVE
# ggsave("simul_output/somalia_output/satellite_newburials_weekly.png",units="cm",height=12.5,width=25)

# compare to reported deaths
weekly_deaths_sel=data.frame((subset(covid_somal,name %in% "death") %>%
  mutate(sum_weekly=roll_sum(cases,7,fill=NA,align="right")))[,c("date","cases","sum_weekly")],datasource="reported") %>% 
  mutate(datasource=as.character(datasource)); colnames(weekly_deaths_sel)[2]="value_report_daily"
weekly_deaths_sel[,"value_norm_report_weekly"]=weekly_deaths_sel$sum_weekly/max(weekly_deaths_sel$sum_weekly,na.rm = T)
weekly_deaths_sel = weekly_deaths_sel %>% mutate(value_norm_rollmean_weekly=roll_mean(value_norm_report_weekly,7,align="center", fill=NA))
# plot together
df_satell_report=weekly_deaths_sel[,!grepl("daily",colnames(weekly_deaths_sel))]
df_satell_report=bind_rows(df_satell_report,newburials_weekly[,!grepl("epi_year|week",colnames(newburials_weekly))]) %>%
  pivot_longer(cols=!c(date,datasource)) %>% filter(!is.na(value)) %>% mutate(type=ifelse(grepl("norm",name),"norm","abs"),
                                                                        daily_rollmean=ifelse(grepl("roll",name),"rollmean","daily"))
ggplot(df_satell_report,aes(x=date,y=value,group=datasource,color=datasource)) + geom_line() + geom_point(aes(shape=datasource),size=0.9) + # 
 facet_grid(type~daily_rollmean,scales = "free") + theme_bw() + standard_theme + theme(axis.text.x=element_text(size=7)) +
 scale_x_date(date_breaks="2 weeks") + ggtitle(expression("new burials vs reported deaths"))
# SAVE
# ggsave("simul_output/somalia_output/satellite_burials_reported_deaths.png",units="cm",height=15,width=25)

# Rt estimate: https://epiforecasts.io/covid/posts/national/somalia/

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### simple SEIR model for Somalia satellite data ----------------------
# ODEs in fcn file

# Oxford Stringency Index
# Gov stringency index
OxCGRT_somalia=subset(read_csv("https://raw.githubusercontent.com/OxCGRT/covid-policy-tracker/master/data/OxCGRT_latest.csv"),
                      CountryName %in% "Somalia")
# PLOT
sel_cols=(unlist(lapply(OxCGRT_somalia,class)) %in% "numeric") & !grepl("Confirmed|_|Legacy|ForDisplay",colnames(OxCGRT_somalia))
ggplot(OxCGRT_somalia[,sel_cols] %>% pivot_longer(cols=!Date),
  aes(x=Date,y=value,group=name,color=name)) + geom_line() + theme_bw() + standard_theme +  facet_wrap(~name)
#   theme(legend.text=element_text(size=6),legend.position = "bottom")
#
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
d_e=mean(rgamma(1e4,shape=4,scale=4/4)); d_c=mean(rgamma(n=1e4,shape=4,scale=3.5/4)) # d_p=mean(rgamma(1e4,shape=4,scale=1.5/4)); 
d_death=21; # sigma=1/d_e; gamma=1/d_c
# TRANSMISSION parameter (same as susceptibility) (without popul normalisation, done in the fcn)
beta_val=0.5; R0=beta_val/(1/d_c)
# IFR
cntr_ex="Somalia"; pop_struct=popF[popF$name %in% cntr_ex,"2020"] + popM[popM$name %in% cntr_ex,"2020"]
IFR_estim=0.5/100 # IFR_estim=c((ifr_all_age_groups$mean %*% pop_struct/sum(pop_struct) ))
# scan values for date of introd: introd_date_scanvals=seq(as.Date("2019-11-21"),as.Date("2020-02-01"),16)
# scan values for seed size: seedsize_scanvals=seq(2,62,15); seeding_duration=10
# compliance with NPIs: compliance_val=0.5
var_categ_list=list(name_vars=c("t", "S", "E","I_R","I_D","R","D","newinfvar","beta"),
    sel_vars=c("t","S","I_R","R","D","new_infections","new_deaths","symptom_cases"), # "I_D",
    case_vars=c("S","I_R","R","symptom_cases","new_infections"),cumul_var=c("S","R","D"), 
    delta_var=c("new_infections","new_deaths"))
# param table of all permutations
scan_params=list(seedsize_scanvals=seq(1,20,4),introd_date_scanvals=seq(as.Date("2019-10-01"),as.Date("2019-12-21"),14),
    seeding_duration=c(seq(21,70,12),90),betaval=seq(0.25,0.4,0.05),compliance=c(0,1,5,10,20,50)/100,IFR=seq(1,5,1)/100)
param_table=expand.grid(scan_params); colnames(param_table)=names(scan_params)
# RUN PARSCAN with parallel (foreach) or sequentially (for) | write "full" for k_lim to scan whole parameter table, otherwise a number
cores=detectCores(); cl<-makeCluster(cores[1]-1); registerDoParallel(cl) # clusterExport(cl, c("fcn_somal_sir_singlesimul")); 
ptm<-proc.time()
df_ode_solution_scan=fcn_paramscan_parallel_seq("paral",k_lim="full",sir_varnames,var_categ_list,timesteps=seq(0,nrow(OxCGRT_input),by=0.25),
            param_table,transm_pars=c(d_c,d_e,d_death),clinical_fraction,OxCGRT_input,N_tot,subpopul=mogadishu_popul,newburials_weekly)
time_sim=proc.time()-ptm; time_sim
stopCluster(cl)

write_csv(df_ode_solution_scan,paste0("simul_output/somalia_output/df_ode_solution_scan",prod(sapply(scan_params, length)),".csv"))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### plot one simulation -----
# timecourse for one param set
singlesimul <- fcn_somal_sir_singlesimul(sir_varnames,var_categ_list,timesteps,
            day0=as.Date("2019-11-01"), num_params=c("beta"=0.325,"gamma"=1/d_c,"sigma"=1/d_e,"d_death"=1/d_death,"IFR"=2.5/100),
            clin_fract=clinical_fraction,seed_size_val=16,seeding_duration=50, # param_table$seeding_duration[k]
            OxCGRT_input,compliance=0.1, # param_table$compliance[k],
            N_tot*mogadishu_popul/sum(N_tot),xlimvals=c(1,350),"plot"); singlesimul[[2]]
# SAVE
# filetag_tcourse=paste0(c(paste0(names(num_params)[c(1,4)],"_",round(num_params,3)[c(1,4)]*c(1,1e2)),"compliance",compliance_val,"day",introd_date),collapse="_")
# ggsave(paste0("simul_output/somalia_output/SIR_timecourse_",filetag,".png"),width=25,height=20,units="cm")
# PEAK cases and deaths
peak_dates=subset(singlesimul[[1]],name %in% c("new_infections","new_deaths")) %>% group_by(name) %>% 
  summarise(maxval_date=date[value==max(value,na.rm=T)])
ggplot(subset(singlesimul[[1]], name %in% c("new_infections","new_deaths")),aes(x=date,y=value,group=name,color=name)) + 
  geom_line() + geom_point() + scale_x_date(limits=as.Date(c("2020-01-01","2020-09-21")),date_breaks="1 week") + scale_y_log10() +
  geom_vline(data=peak_dates,aes(xintercept=maxval_date,color=name)) + theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5))
# ggsave(paste0("simul_output/somalia_output/peak_cases_deaths",filetag,".png"),width=25,height=15,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

### discretise number of new deaths (for parscan) -----
# ESTIMATE of death rate for mogadishu: 0.4 deaths/10K ppl/day -> Mogadishu is 2.2M -> 88 deaths/day -> ~1250/2 weeks
crude_deathrate_per1e4personday=0.4; mogad_weekly_death_rate=(mogadishu_popul/1e4)*crude_deathrate_per1e4personday*7
df_ode_solution_scan=df_ode_solution_scan %>% mutate(value_norm=sum_weekly/mogad_weekly_death_rate, datasource="model (excess deaths %)",
  parset_ID=group_indices_(df_ode_solution_scan,.dots=rev(c("introd_date","seed_size","IFR","npi_compliance","beta","seeding_duration"))))
# basal_burial_rate=mean(newburials_weekly$new_graves_best_ipol[newburials_weekly$date<="2020-01-19"])
# concatenate data types + simul --> calculate ERROR
# NOT fitting the period before the rise of burials begin -> as.Date("2020-01-19")
error_pred_data_comp=left_join(subset(df_ode_solution_scan,!is.na(sum_weekly) & date>as.Date("2020-01-19")), # when model predicts 0 no fitting
  newburials_weekly[,!grepl("epi_year|week",colnames(newburials_weekly))],by="date",suffix=c("_model","_satel_data")) %>%
  group_by(introd_date,seed_size,seeding_duration,IFR,npi_compliance,beta,parset_ID) %>% # ,weekly_deaths_sel,by="date",suffix=c("","_report_data"))
  summarise(rmse=sqrt(mean((round((new_graves_best_ipol-basal_burial_rate)*(mogad_weekly_death_rate/basal_burial_rate))-round(sum_weekly))^2,na.rm=T)),
    loglh_poiss=-sum(dpois(x=round((new_graves_best_ipol-basal_burial_rate)*(mogad_weekly_death_rate/basal_burial_rate)),lambda=round(sum_weekly),log=TRUE) ),
    loglh_negbinom=-sum(dnbinom(x=round((new_graves_best_ipol-basal_burial_rate)*(mogad_weekly_death_rate/basal_burial_rate)),mu=round(sum_weekly),log=TRUE,size=10) )) %>% 
  pivot_longer(cols=!c(introd_date,seed_size,IFR,npi_compliance,beta,seeding_duration,parset_ID))
# mse_reported=mean((value_norm_model - value_norm_report_weekly)^2,na.rm=T)
# compare likelihood with MSE
err_types=data.frame(sq_err=error_pred_data_comp$value[grepl("mse",error_pred_data_comp$name)],
     poisson=error_pred_data_comp$value[grepl("poiss",error_pred_data_comp$name)],
     negbinom=error_pred_data_comp$value[grepl("negbinom",error_pred_data_comp$name)])
ggplot(err_types) + geom_point(aes(x=sq_err,y=poisson)) + geom_point(aes(x=sq_err,y=negbinom),color="red") + theme_bw() + standard_theme +
  scale_x_log10() + scale_y_log10() + ylab("neg. log-likelihood") + # limits=c(2e4,1e5) | limits=c(2e2,5e3)
  xlab("mean squared error") + labs(title="MSE vs negative log-likelihood",subtitle="red=neg. binom., black=poisson")
# ggsave("simul_output/somalia_output/errors/logll_mse_comparison.png",width=32,height=20,units="cm")

### PLOTALL trajectories for best fits ------------------------
err_type="negbinom"; n_err_round=0
best_fit=subset(error_pred_data_comp, grepl(err_type,name)) %>% group_by() %>% slice(which.min(value))
parcol_names=c("introd_date","seed_size","IFR","npi_compliance","beta","seeding_duration","parset_ID")
df_plot=left_join(subset(df_ode_solution_scan, IFR==best_fit$IFR & seeding_duration==best_fit$seeding_duration & beta==best_fit$beta ),
                  subset(error_pred_data_comp,grepl(err_type,name)),by=parcol_names,suffix=c("","_error"))
min_error_vals=df_plot %>% group_by(seed_size,introd_date) %>% summarise(min_error=min(value_error))
ggplot(df_plot) + geom_line(aes(x=date,y=value_norm,group=npi_compliance,color=factor(npi_compliance))) +
 geom_line(data=newburials_weekly,aes(x=date,y=value_norm),color="black",linetype="dashed") + # ,size=1.05 satellite
 facet_grid(introd_date~seed_size,labeller=labeller(seed_size=label_both),scales="free") + 
 theme_bw() + standard_theme + theme(legend.position="top",axis.text.x=element_text(size=7),axis.text.y=element_text(size=8)) +
 geom_vline(aes(xintercept=introd_date),color="red") + geom_vline(xintercept=as.Date(min(OxCGRT_input$date[OxCGRT_input$NPI_on>0])),color="green") +
 geom_text(data=min_error_vals,aes(label=round(min_error,n_err_round)),x=as.Date("2020-08-10"),y=1.4,size=3) + 
 scale_x_date(date_breaks="2 weeks",limits=c(min(OxCGRT_input$date)+15,max(newburials_weekly$date)+30)) + 
 xlab('time (days)') + ylab('normalised rate') + labs(title="SIR model ~ f(introduction date & seed size)",
      subtitle=paste0("beta=",unique(df_plot$beta),", IFR=",unique(df_plot$IFR)*1e2,"%, seeding duration=",best_fit$seeding_duration,
                      ", error=",err_type), color="NPI compliance")
# SAVE
# filetag=paste0("_beta_",unique(df_plot$beta),"_IFR_",unique(df_plot$IFR)*1e2,"_seeding_t_",best_fit$seeding_duration)
# ggsave(paste0("simul_output/somalia_output/param_scan_satellite_SIRmodel_bestfit_",err_type,filetag,".png"),width=40,height=24,units="cm")

### ### ### ### ### ### ### ### ### ### ### ###
### heatmap of MSE ------------------------
# PLOT
err_type="negbinom"; n_text=2; n_round=0
for (k in 1:length(scan_params$IFR)){
df_heatmap=subset(error_pred_data_comp,IFR==scan_params$IFR[k] & grepl(err_type,name)) #  & beta==scan_params$betaval[2]
colnames(df_heatmap)[grepl("npi",colnames(df_heatmap))]="npi_c"; colnames(df_heatmap)[grepl("durat",colnames(df_heatmap))]="seed_t"
# best 1% param sets
# minval_mse=min(subset(error_pred_data_comp,grepl(err_type,name))[,"value"]) 
min_error_limit=10^(floor(log10(sort(subset(error_pred_data_comp,grepl(err_type,name))$value)[round(0.01*nrow(error_pred_data_comp))])*10)/10 )
hmap_lims=summary(subset(error_pred_data_comp,grepl(err_type,name) & !is.infinite(value))$value)[c("Min.","Median","Max.")]*c(0.9,1,1.1)
highl_tiles=subset(df_heatmap, value<min_error_limit) # %>% select(introd_date,seed_size) # midpoint_val=2*min(df_heatmap$value)
# PLOT
ggplot(df_heatmap, aes(x=seed_size,y=introd_date,fill=value)) + geom_tile(color="black") + 
 geom_text(aes(label=round(value,n_round)),color="black",size=n_text) + geom_rect(data=highl_tiles,size=1,fill=NA,colour="green",
  aes(xmin=seed_size-unique(diff(scan_params$seedsize_scanvals))/2,xmax=seed_size+unique(diff(scan_params$seedsize_scanvals))/2,
  ymin=introd_date-unique(diff(scan_params$introd_date_scanvals))/2,ymax=introd_date+unique(diff(scan_params$introd_date_scanvals))/2)) +
  theme_bw() + standard_theme + theme(axis.text.y=element_text(size=6)) + 
  scale_fill_gradientn(colors=c("blue","white","red"),values=scales::rescale(array(hmap_lims))) + #,limits=c(0,hmap_lims[3]) 
  facet_grid(npi_c~seed_t~beta,labeller=labeller(npi_c=label_both,seed_t=label_both,beta=label_both)) + 
  scale_y_date(breaks=unique(df_heatmap$introd_date),expand=c(0,0)) + scale_x_continuous(breaks=unique(df_heatmap$seed_size),expand=c(0,0)) +
  labs(title="Error for grid search",subtitle=paste0("IFR=",unique(df_heatmap$IFR)*100,"%, error: ",
  gsub(gsub(err_type,"poiss","NLL, poisson"),"negbinom","NLL, neg.binom."))) + xlab("seed size") + ylab("introd. date")
# SAVE
filetag=paste0("IFR_",1e2*unique(df_heatmap$IFR),"_error_",err_type) 
ggsave(paste0("simul_output/somalia_output/errors/error_parscan_SIR_heatmap_",filetag,".png"),width=32,height=40,units="cm")
print(k)
}

### ### ### ### ### ### ### ### ### ### ### ###
### profile likelihood from grid search
profile_likelihood=error_pred_data_comp %>% mutate(introd_date=as.numeric(introd_date)) %>% # %>% select(!contains("name")) 
 rename(error_type=name,error=value) %>% pivot_longer(cols=c(introd_date,seed_size,seeding_duration,IFR,npi_compliance,beta)) %>% # 
 group_by(error_type,name,value) %>% summarise(best_fit=min(error)) %>% mutate(value=ifelse(name %in% c("IFR","npi_compliance"),value*1e2,value))
# to have dates as strings not numbers
# profile_likelihood=left_join(profile_likelihood,data.frame(date_orig=as.character(unique(df_ode_solution_scan$introd_date)),
#              value=as.numeric(unique(df_ode_solution_scan$introd_date))),by="value")
# profile_likelihood$date_orig=as.character(profile_likelihood$date_orig)
# profile_likelihood$date_orig[is.na(profile_likelihood$date_orig)]=profile_likelihood$value[is.na(profile_likelihood$date_orig)]
# univals=profile_likelihood$date_orig[1:(nrow(profile_likelihood)/3)] 
# profile_likelihood$date_orig=factor(profile_likelihood$date_orig)

ggplot(profile_likelihood,aes(x=value,y=best_fit,group=1)) + geom_line() + geom_point() + 
  facet_grid(error_type~name,scales="free") + 
  # facet_wrap(error_type~name,scales="free",ncol=length(scan_params)) + 
  theme_bw() + standard_theme + xlab("parameter") +
  scale_y_log10() + labs(title="Profile likelihood (MSE, negLogL-Poisson, negLogL-negat.binom.)",
                         subtitle="best fit for given parameter value (from grid search)")
# SAVE
ggsave(paste0("simul_output/somalia_output/errors/profile_likelihood_october_facetgrid.png"),width=32,height=18,units="cm") # _facetgrid

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### optimisation --------------
# ODE in fcns file

# DATA | # basal_burial_rate=min(newburials_weekly$new_graves_best_ipol)
fitting_data=data.frame(date=newburials_weekly$date, t=as.numeric(newburials_weekly$date)-as.numeric(newburials_weekly$date)[1],
            new_burials_weekly=(newburials_weekly$new_graves_best_ipol-basal_burial_rate)*mogad_weekly_death_rate/basal_burial_rate)
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
