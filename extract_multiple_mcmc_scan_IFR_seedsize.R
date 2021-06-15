### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### multiple fits with fixed seed size and compliance levels
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
lapply(c("tidyverse","deSolve","qs","gtools","rstudioapi","wpp2019","countrycode","coronavirus","wesanderson","dttr2","RcppRoll",
         "scales","wpp2019","GGally","corrr","ungeviz"), library,character.only=TRUE)
# functions and plotting theme
source("somalia_data_model_fcns.R")
# load folder (change to "repo_data" or own folder)
parscan_mcmc_dirname=
"simul_output/somalia/3param_fits_seedsize_IFR_fixed/scan_seedsize_ifr_introddate_N_182_20_fitperiod_20200223_20200824/fits_seeding_30_70y/"
parfit_scan_files<-list.files(parscan_mcmc_dirname,pattern = ".rds"); # slope_val=round(as.numeric(linregr$coefficients[2]),4)
# need to have IFR estimates from Sandmann: IFR_estimates_Sandmann2021$logit_ifr
# how many CDR values were used?
onefit=readRDS(paste0(parscan_mcmc_dirname,parfit_scan_files[1]))[[1]]; x_dodge_val=0.6; fitting_params <- names(onefit$priors)
scan_params<-c("seedsize","ifr_logit_increm","IFR all infections (%)")
# load burial data
burial_data=read_csv("repo_data/out_bdr_daily_estimates.csv")
baseline_daily_burials=mean(subset(burial_data,date>="2019-07-01" & date<="2019-11-01")$new_graves_best_ipol)
out_bdr_daily_estimates=burial_data[!rowSums(is.na(burial_data))==(ncol(burial_data)-1),
                                    !colSums(is.na(burial_data))==nrow(burial_data)] %>% filter(date>"2019-11-01") %>%
  mutate(daily_baseline_subtr=ifelse(new_graves_best_ipol-baseline_daily_burials>0,new_graves_best_ipol-baseline_daily_burials,0),
         rollmeanweek=roll_mean(daily_baseline_subtr,7,align="center", fill=NA), # rolling mean BASELINE subtracted
         rollsumweek=roll_sum(daily_baseline_subtr,7,align="left",fill=NA),
         rollmeanweek_no_subtr=roll_mean(new_graves_best_ipol,7,align="center", fill=NA),
         rollsumweek_no_subtr=roll_sum(new_graves_best_ipol,7,align="left",fill=NA))
# CDR (crude death rate) value: fitting only with estimate from data
mogadishu_popul=2.2e6
CDR_vals=c(baseline_daily_burials*1e4/mogadishu_popul,0.1,0.2,0.4)[1:length(readRDS(paste0(parscan_mcmc_dirname,parfit_scan_files[1])))]
# load IFR
data(pop)
somalia_agegroups_IFR=fcn_merge_ifr_above_age(left_join(fcn_load_age_str("Somalia",n_year="2015",90),
   fcn_load_ifr("repo_data/IFR_by_age_imperial.csv"),by=c("agegroup","agegroup_min")),75) %>% 
  mutate(ifr_mean=ifelse(ifr_mean==0,min(ifr_mean[ifr_mean>0]),ifr_mean),log_ifr=log(ifr_mean),logit_ifr=log(ifr_mean/(1-ifr_mean)))
# from Sandmann 2021 cmmid paper
IFR_estimates_Sandmann2021<-read_csv("repo_data/IFR_estimates_Sandmann2021.csv")
if (any(IFR_estimates_Sandmann2021$value_percent>1)) {n_cols<-2:ncol(IFR_estimates_Sandmann2021)
IFR_estimates_Sandmann2021[,n_cols]<-IFR_estimates_Sandmann2021[,n_cols]/1e2; 
IFR_estimates_Sandmann2021 <- left_join(IFR_estimates_Sandmann2021 %>% rename(agegroup=Age,ifr_mean=value_percent), 
   somalia_agegroups_IFR %>% select(!c(ifr_mean,log_ifr,logit_ifr)),by="agegroup") %>% mutate(logit_ifr=log(ifr_mean/(1-ifr_mean))) }
### ### ### ### ### ### ### ### ### ### 
# read in NPI data Oxford Stringency Index
# truncate until a given timepoint
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
### ### ###
# load COVIDM parameters
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
cm_force_rebuild=F; cm_build_verbose=T; cm_version=1; source(file.path(cm_path,"R","covidm.R"))
countryval="Somalia"; params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval),
                                                 date_start="2019-10-01",date_end="2020-10-01")
# set population: Somalia --> Mogadishu (these should be the same as used for the fitting!)
params$pop[[1]]$name=onefit$base_parameters$pop[[1]]$name
params$pop[[1]]$size=onefit$base_parameters$pop[[1]]$size
params$pop[[1]]$dist_seed_ages=onefit$base_parameters$pop[[1]]$dist_seed_ages
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# read in results of fitting
scan_params<-c("seedsize","ifr_logit_increm","IFR all infections (%)")
CDR_vals=c(baseline_daily_burials*1e4/mogadishu_popul,0.1,0.2,0.4)[1:length(readRDS(paste0(parscan_mcmc_dirname,parfit_scan_files[1])))]
for (k in 1:length(parfit_scan_files)) {
  # for (n_CDR in 1:length(CDR_vals)) {
  x=readRDS(paste0(parscan_mcmc_dirname,parfit_scan_files[k]))[[1]]$posterior
  if (is.null(x)) {x=readRDS(paste0(parscan_mcmc_dirname,parfit_scan_files[k]))$posterior}
  df_posteriors_comb=x %>% select(chain,trial,lp,all_of(fitting_params)) %>% # mutate(CDR=round(CDR_vals[n_CDR],3)) %>%
      mutate(ifr_logit_increm=as.numeric(gsub("ifr_increm|_","",as.numeric(str_match(parfit_scan_files[k],"ifr_increm(.*?)_")[2]))),
          seedsize=as.numeric(str_match(parfit_scan_files[k],"seedsize(.*?).rds")[2]),`IFR all infections (%)`=sapply(ifr_logit_increm,
        function(x) 1e2*sum(inv.logit(IFR_estimates_Sandmann2021$logit_ifr + x)*somalia_agegroups_IFR$agegroup_perc) ) )
    # if (n_CDR==1){df_posteriors_comb=df_posteriors} else {df_posteriors_comb=rbind(df_posteriors_comb,df_posteriors)}
  # }
  if (k==1) {df_posteriors_parscan=df_posteriors_comb} else { df_posteriors_parscan=bind_rows(df_posteriors_parscan,df_posteriors_comb)} 
}
# save
write_csv(df_posteriors_parscan,paste0(parscan_mcmc_dirname,"df_posteriors_parscan.csv"))
### ###
# calc summary stats
posteriors_summary_stats = df_posteriors_parscan %>% mutate(n=row_number()) %>% 
  pivot_longer(!c(n,all_of(scan_params),chain,trial,lp)) %>% group_by(name,ifr_logit_increm,seedsize) %>% 
  summarise(mean=mean(value),median=median(value),ci95_low=quantile(value,probs=c(2.5,97.5)/1e2)[1],
          ci95_up=quantile(value,probs=c(2.5,97.5)/1e2)[2],ci50_low=quantile(value,probs=c(25,75)/1e2)[1],
          ci50_up=quantile(value,probs=c(25,75)/1e2)[2],post_lkl_mean=mean(lp),post_lkl_ci95_low=quantile(lp,probs=c(2.5,97.5)/1e2)[1],
    post_lkl_ci95_up=quantile(lp,probs=c(2.5,97.5)/1e2)[2]) %>% mutate(name=ifelse(name %in% "introd_date",
    paste0("introduction (days after ",format(as.Date(onefit$base_parameters$date0),"%d/%m/%y"),")"),name),
    name=ifelse(name %in% "NPI_scale","NPI NPI_scale (0 to 1)",name)) %>% mutate_if(is.numeric,round,3)
# save
write_csv(posteriors_summary_stats,paste0(parscan_mcmc_dirname,"posteriors_summary_stats.csv"))
# posteriors_summary_stats = posteriors_summary_stats %>% mutate()
### ### ### ### ### ### ### ### ### ### ### ###
# RUN simuls from param distribs of fits
for (k_seedsize in sort(unique(df_posteriors_parscan$seedsize))) {
  for (k_ifr_logit in unique(df_posteriors_parscan$ifr_logit_increm)) {
    sel_file=parfit_scan_files[grepl(paste0("seedsize",k_seedsize,".rds"),parfit_scan_files) &
                                 grepl(paste0("ifr_increm",k_ifr_logit,"_"),parfit_scan_files)]
    # read in file
    fits_death_scaling=readRDS(paste0(parscan_mcmc_dirname,sel_file)) # n_seedsize=k_seedsize; ifr_logit_increm=k_ifr_logit
    if (!is.null(fits_death_scaling[[1]]$posterior)){ fits_death_scaling=fits_death_scaling[[1]] }
        print(paste0("ifr logit increm=",k_ifr_logit,", seedsize=",k_seedsize)) # "CDR scale=",k,
        if (length(unique(fits_death_scaling$base_parameters$pop[[1]]$y))==1){ message(paste0(k_seedsize,",",k_ifr_logit))
        fits_death_scaling$base_parameters$pop[[1]]$y=fun_lin_approx_agedep_par(agegroups=somalia_agegroups_IFR,
                                                              min_val=0.25,max_val=0.7,rep_min=6,rep_max=2)}
        ####
        # these two vars we need to make global bc cm_sample_fit looks for them
        n_seedsize=k_seedsize; ifr_logit_increm=k_ifr_logit
        dyn=cm_sample_fit(fits_death_scaling,25) %>% mutate(ifr_logit_increm=k_ifr_logit,seedsize=k_seedsize) %>%
        filter(compartment %in% c("S","cases","death_o"))
    if (k_seedsize==sort(unique(df_posteriors_parscan$seedsize))[1] & k_ifr_logit==unique(df_posteriors_parscan$ifr_logit_increm)[1]){
      dyn_all=dyn} else {dyn_all=rbind(dyn_all,dyn)}  } }
# summarize these runs
summ = dyn_all %>% group_by(t,run,compartment,ifr_logit_increm,seedsize) %>% summarise(value=sum(value)) %>% # ,CDR
  group_by(run,compartment,ifr_logit_increm,seedsize) %>% mutate(value=ifelse(compartment %in% "S",
   (mogadishu_popul-value)/mogadishu_popul,value),compartment=ifelse(compartment %in% "S","attack_rate",as.character(compartment))) %>%
  group_by(t,compartment,ifr_logit_increm,seedsize) %>% summarise(lower=hdi(value)[[1]],upper=hdi(value)[[2]],mean=mean(value)) %>%
  mutate(date=as.Date(seq(as.Date(fits_death_scaling$base_parameters$date0),
                          as.Date(fits_death_scaling$base_parameters$date0)+max(t),1)[t+1])) %>%
  mutate(ifr_all_inf=round(sapply(ifr_logit_increm,
            function(x) 1e2*sum(inv.logit(IFR_estimates_Sandmann2021$logit_ifr+x)*somalia_agegroups_IFR$agegroup_perc)),2) )
# SAVE
write_csv(summ,paste0(parscan_mcmc_dirname,"summ.csv")); rm(dyn_all)

# # cumulative sympt cases
# ggplot(summ %>% filter(compartment=="cases") %>% group_by(ifr_logit_increm,seedsize) %>% arrange(date) %>%
#  mutate(cumul_cases=cumsum(mean)/2.2e6),aes(x=date,y=cumul_cases)) + geom_line() + facet_grid(ifr_logit_increm~seedsize,scale="free") +
#    standard_theme

# compare to data
fitting_incidence_modelcompare=data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
  mutate(new_deaths=round(daily_baseline_subtr),
         date_within_fitting_t=ifelse(date>=min(onefit$data$date) & date<=max(onefit$data$date),TRUE,FALSE))) 
# calculate likelihood, deviance, DIC
DIC_logllk_values <- left_join(right_join(summ %>% filter(compartment=="death_o") %>% select(!c(lower,upper)),
      subset(fitting_incidence_modelcompare,date_within_fitting_t) %>% select(!date_within_fitting_t),by=c("date")) %>% 
      mutate(logllk=ifelse(!is.infinite(dpois(new_deaths,lambda=mean,log=T)),dpois(new_deaths,lambda=mean,log=T),NA)) %>% 
      group_by(ifr_logit_increm,seedsize) %>% summarise(sum_logllk=sum(logllk,na.rm=T), # out_bdr_daily_estimates
      deviance=-2*sum_logllk,d_p=var(-2*logllk,na.rm=T)/2,DIC=deviance+d_p) %>% mutate(dummy=1), # JOINED WITH
      fitting_incidence_modelcompare %>% summarise(maxval=max(new_deaths)) %>% mutate(dummy=1),by="dummy") %>% select(!dummy) %>% 
   mutate(DIC_str=round(DIC) ) %>% mutate(ifr_all_inf=round(sapply(ifr_logit_increm,
          function(x) 1e2*sum(inv.logit(IFR_estimates_Sandmann2021$logit_ifr+x)*somalia_agegroups_IFR$agegroup_perc)),2) )
# paste0(round(DIC/10^floor(log10(DIC)),2),"e",floor(log10(DIC)))

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### correlations between parameters
for (k_seedsize in unique(df_posteriors_parscan$seedsize)) {
  for (k_ifr_logit in unique(df_posteriors_parscan$ifr_logit_increm)) {
    pairwise_corrs = # bind_rows( lapply(unique(df_posteriors_parscan$CDR), function(x)  CDR==x,
      cbind(df_posteriors_parscan %>% filter(seedsize==k_seedsize,ifr_logit_increm==k_ifr_logit) %>% select(all_of(fitting_params)) %>% 
        corrr::correlate() %>% corrr::shave() %>% corrr::stretch(),seedsize=k_seedsize,ifr_logit_increm=k_ifr_logit) %>% # )) 
 mutate(x=factor(x,levels=unique(x)),y=factor(y,levels=unique(y)),highlight_color=ifelse(abs(r)<0.8|is.na(r),0,1)) %>% filter(!is.na(r))
    if (k_seedsize==unique(df_posteriors_parscan$seedsize)[1] & k_ifr_logit==unique(df_posteriors_parscan$ifr_logit_increm)[1]) {
     pairwise_corrs_parscan = pairwise_corrs} else { pairwise_corrs_parscan = bind_rows(pairwise_corrs_parscan,pairwise_corrs)} } }

# plot pairwise correlations
ggplot(pairwise_corrs_parscan,aes(x,y,fill=r))+geom_tile(color="black") + # aes(color=factor(highlight_color)),size=2,width=0.96,height=0.96
  facet_grid(seedsize~ifr_logit_increm,labeller=labeller(seedsize=label_both,ifr_logit_increm=label_both)) + # CDR~
  geom_text(aes(label=round(r,2)),size=5) + scale_x_discrete(expand=c(0,0.02)) + scale_y_discrete(expand=c(0,0.02)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",limit=c(-1,1)) + # scale_color_manual(values=c("grey","black"),guide=FALSE) +
  standard_theme + theme_bw() + xlab("") + ylab("") +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.95),
        axis.ticks = element_blank())
# SAVE
if (!dir.exists(paste0(parscan_mcmc_dirname,"mcmc_diagnostics"))) {dir.create(paste0(parscan_mcmc_dirname,"mcmc_diagnostics"))}
ggsave(paste0(parscan_mcmc_dirname,"mcmc_diagnostics/param_correlations.png"),width=36,height=16,units="cm")

# scatterplots of fitted parameters: (R0_fit,ifr_logit_intercept) (R0_fit,introd_date) (introd_date,ifr_logit_intercept)
var_pairs <- list(fitting_params[c(1,2)],fitting_params[c(1,3)],fitting_params[c(2,3)])
lp_stats <- df_posteriors_parscan %>% summarise(min_lp=min(lp),mean_lp=mean(lp),max_lp=max(lp))
for (k in 1:length(var_pairs)){
  varpair_fl=gsub("_fit","",var_pairs[[k]])
  p <- ggplot(df_posteriors_parscan,aes(x=get(var_pairs[[k]][1]),y=get(var_pairs[[k]][2]))) + 
    geom_point(aes(color=lp)) + geom_smooth(size=0.5,color="black") + # ,method="loess"
    scale_color_gradient2(low="red",mid="white",high="blue",midpoint=lp_stats$mean_lp) +
  facet_grid(ifr_logit_increm~seedsize,labeller=labeller(seedsize=label_both,ifr_logit_increm=label_both),
    scales="free") + theme_bw() + standard_theme + # ,nrow=length(unique(df_posteriors_parscan$ifr_logit_increm))
    labs(color="log-posterior") + xlab(varpair_fl[1]) + ylab(varpair_fl[2]); p
  # SAVE
  ggsave(paste0(parscan_mcmc_dirname,"mcmc_diagnostics/scatterplot_",varpair_fl[1],"_",varpair_fl[2],".png"),
         width=30,height=18,units="cm") }

# plot traces of mcmc wrt post log-llh and param values
for (k in 3:6){
  p <- ggplot(df_posteriors_parscan %>% 
                mutate(rolling_mean=roll_mean(!!as.name(colnames(df_posteriors_parscan)[k]),n=50,align="center",fill=NA)),
    aes(x=trial)) + geom_line(aes(y=get(colnames(df_posteriors_parscan)[k]),color=factor(chain)),size=0.5) +
    # geom_line(aes(y=rolling_mean),size=0.3)+ 
    facet_grid(ifr_logit_increm~seedsize,labeller=labeller(seedsize=label_both,ifr_logit_increm=label_both),
             scales="free") + # ,nrow=length(unique(df_posteriors_parscan$ifr_logit_increm))
    theme_bw() + standard_theme + labs(color="chains") + ylab(colnames(df_posteriors_parscan)[k]); p
  # SAVE
  ggsave(paste0(parscan_mcmc_dirname,"mcmc_diagnostics/trace_",colnames(df_posteriors_parscan)[k],".png"),width=30,height=18,units="cm")
  }

# plot posterior likelihoods
for (k in 4:6){
  p <- ggplot(df_posteriors_parscan,aes(x=get(colnames(df_posteriors_parscan)[k]),y=lp)) + 
   geom_point(aes(color=lp),size=0.4) + scale_color_gradient(low="red",high="blue") + # mid="white", ,midpoint=lp_stats$mean_lp
  geom_smooth(aes(fill=lp),color="black",method="loess",size=0.25) +
  facet_grid(ifr_logit_increm~seedsize,labeller=labeller(seedsize=label_both,ifr_logit_increm=label_both),
             scales="free") + # ,nrow=length(unique(df_posteriors_parscan$ifr_logit_increm))
  theme_bw() + standard_theme + labs(color="log-posterior",fill="log-posterior") + xlab(colnames(df_posteriors_parscan)[k]) +
  ylab("log-posterior"); p
  # SAVE
 ggsave(paste0(parscan_mcmc_dirname,"mcmc_diagnostics/logpost_",colnames(df_posteriors_parscan)[k],".png"),width=30,height=18,units="cm")
}
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT DYNAMIC fits faceted by (CDR -seedsize-NPI_scale)
fitting_dates=fits_death_scaling$data %>% summarise(min_date=min(date),max_date=max(date))
ggplot(summ %>% filter(compartment=="death_o")) + 
  geom_line(aes(x=date,y=mean), colour="blue") + geom_ribbon(aes(x=date,ymin=lower,ymax=upper), fill="blue",alpha=0.2) + 
  geom_line(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),linetype="dashed",size=0.3) +
  geom_point(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),fill=NA,shape=1,size=0.2) +
  facet_grid(seedsize~ifr_logit_increm,labeller=labeller(seedsize=label_both,ifr_logit_increm=label_both),scales="free") +
  geom_rect(data=subset(DIC_logllk_values,DIC<1.02*min(DIC)),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill=NA,color="red",size=1.2) +
  xlab("") + ylab("deaths (daily)") + standard_theme + theme_bw() + theme(axis.text.x=element_text(vjust=0.5,angle=90)) + 
  geom_text(data=DIC_logllk_values,aes(x=fitting_dates$min_date,y=0.9*maxval,label=paste0("DIC=",round(DIC))),color="black",size=3.5)+
  scale_x_date(date_breaks="4 week",limits=as.Date(c("2020-01-01","2020-10-01")),expand=expansion(0.0)) + 
  scale_y_continuous(breaks=(0:10)*2,expand=expansion(0,0))
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"all_dynamic_fits.png"),width=38,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### 
# fits with DIC values
fitting_params_best_estim = left_join(posteriors_summary_stats %>% # mutate(CDR=ifelse(CDR==0.042,round(CDR_vals[1],5),CDR)) %>% 
  # check start date of simul, rewrite if needed!!
  filter(name %in% c("npi_scale","R0_fit","introduction (days after 01/09/19)")) %>% 
  select(name,seedsize,ifr_logit_increm,median),DIC_logllk_values %>% select(!c(deviance,d_p)),by=c("seedsize","ifr_logit_increm")) %>%
  pivot_wider(names_from=name,values_from=median) %>% # rename(IFR=`IFR all infections (%)`) %>%
  mutate(introd_date=as.Date(onefit$base_parameters$date0)+`introduction (days after 01/09/19)`) %>%
  mutate(ifr_all_inf=round(sapply(ifr_logit_increm,
                  function(x) 1e2*sum(inv.logit(IFR_estimates_Sandmann2021$logit_ifr+x)*somalia_agegroups_IFR$agegroup_perc) ),2))
# contact_red=round(max(npi_df$contact_reduction)*NPI_scale,2)

# dynamic fits faceted by (CDR-seedsize), NPI_scale by color
y_text=14.5; x_dodge_text=5; y_low=7.5
ggplot(summ %>% filter(compartment=="death_o")  ) +
  geom_line(aes(x=date,y=mean,color=factor(round(ifr_all_inf,2)))) + # ,group=interaction(NPI_scale,CDR),alpha=factor(CDR)
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=factor(round(ifr_all_inf,2))),alpha=0.2) + # alpha=factor(CDR)
  geom_line(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),linetype="dashed",size=0.3) +
  geom_point(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),fill=NA,shape=1,size=0.2) +
  facet_wrap(~seedsize,labeller=labeller(seedsize=label_both),scales="free",nrow=2) + # CDR~seedsize
  scale_alpha_manual(values=c(0.4,0.45,0.5,0.55)) + xlab("") + ylab("deaths (daily)") + standard_theme + theme_bw() +
  theme(axis.text.x=element_text(vjust=0.5,angle=90),panel.grid.minor=element_blank(),legend.position="top") +
  # DIC labels
  geom_text(data=DIC_logllk_values,aes(x=as.Date("2020-01-15")+x_dodge_text-3,y=as.numeric(factor(ifr_logit_increm))*1.4+y_text,
                                       label="DIC=",color=factor(ifr_all_inf)),size=3.25,show.legend=FALSE) + # DIC values
  geom_text(data=DIC_logllk_values,aes(x=as.Date("2020-01-15")+x_dodge_text+20,
      y=as.numeric(factor(ifr_all_inf))*1.4+y_text,label=paste0(DIC_str)),size=3.25,show.legend=FALSE) + # introd date labels
  geom_text(data=fitting_params_best_estim,aes(x=as.Date("2020-01-15")+x_dodge_text,y=as.numeric(factor(ifr_all_inf))*1.4+y_text-y_low,
                                               label="introd.=",color=factor(ifr_all_inf)),size=3.25,show.legend=FALSE) + 
  geom_text(data=fitting_params_best_estim,aes(x=as.Date("2020-01-15")+x_dodge_text+32, # introd dates
        y=as.numeric(factor(ifr_all_inf))*1.4+y_text-y_low,label=gsub("-","/",format(introd_date,"%d/%m/%y"))),size=3.25,show.legend=FALSE) +
  scale_x_date(date_breaks="2 week",limits=as.Date(c("2020-01-01","2020-10-01")),expand=expansion(0.0)) + 
  geom_vline(data=npi_df,aes(xintercept=on),linetype="dashed",show.legend=F,size=0.5,color="darkgrey") + # NPIs
  scale_y_continuous(breaks=(0:12)*2,limits=c(0,22.5),expand=expansion(0,0)) + labs(color="IFR",fill="IFR")
  # NPI labels # geom_text(data=npi_df[1:3,],aes(x=on+ifelse(name=="first",38,10),y=21.5,label=paste0(ifelse(name=="first","max(Stringency)=",""),
  #     round(contact_reduction,2))),size=3) + 
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"all_dynamic_fits_NPIcolorcode.png"),width=36,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT DIC as function of CDR~seedsize~ifr NPI_scale
label_y_dodge=12; label_y_factor=1.03; fontsize=4.5
ggplot(fitting_params_best_estim %>% filter(seedsize>20) %>% # ifr_all_inf<1.6 & 
         mutate(R0_label=paste0("R0=",round(R0_fit,1)),npi_max_effect=npi_scale*max(npi_df$contact_reduction)) ) + 
  geom_hpline(aes(x=factor(ifr_all_inf),y=DIC,group=seedsize,color=factor(seedsize)),width=0.24, # group=CDR,alpha=factor(CDR)
              position=position_dodge(width=1)) + scale_alpha_manual(values=c(0.4,0.55,0.7,0.85)) + 
  geom_vline(xintercept=(1:4)+0.5,size=0.2,linetype="dashed") + # facet_wrap(~seedsize,labeller=labeller(CDR=label_both),nrow = 3) +
  scale_x_discrete(expand=expansion(0.1,0)) + # scale_y_log10() + # scale_y_continuous(breaks=seq(5,100,by=1)*100) + 
  theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(vjust=0.5,hjust=0.95,size=17),axis.text.y=element_text(size=17),panel.grid.major.x=element_blank(),
        legend.position="top",strip.text=element_text(size=15),axis.title=element_text(size=20),legend.text=element_text(size=15),
        legend.title=element_text(size=17)) +
  # introd date
  geom_text(aes(x=factor(ifr_all_inf),y=DIC+label_y_dodge,group=seedsize,label=gsub("-","/",format(introd_date,"%d/%m/%y"))),
            position=position_dodge(width=1),size=fontsize) +
  # R0
  geom_text(aes(x=factor(ifr_all_inf),y=DIC+label_y_dodge+14,group=seedsize,label=R0_label),position=position_dodge(width=1),size=fontsize) + 
  geom_text(aes(x=factor(ifr_all_inf),y=DIC-label_y_dodge,group=seedsize,label=paste0(round(npi_max_effect*100),"%")),
            position=position_dodge(width=1),size=fontsize) + xlab("IFR (%)") + ylab("DIC") + labs(color="seedsize",alpha="CDR")
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"DIC_xaxis_NPIscaling_colorcode_seedsize_alpha_CDR.png"),width=40,height=24,units="cm")

### ### ### ### ### ### ###
# PLOT facet_grid NAME-NPI, x-axis CDR, color coded by SEEDSIZE
round_prec=1
df_summary_plot=left_join(posteriors_summary_stats %>% group_by(name) %>% mutate(name=ifelse(name=="R0_fit","R0",name),
  median_range=max(median)-min(median)),DIC_logllk_values %>% select(!c(maxval,deviance,d_p)),by=c("seedsize","ifr_logit_increm")) %>% 
  filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")) %>% mutate(DIC_str=round(DIC),
    introd_date_str=gsub("2019","20",gsub("-","/",format(as.Date(onefit$base_parameters$date0)+median,"%d/%m/%y")))) %>%
  mutate(ifr_all_inf=round(sapply(ifr_logit_increm,
                function(x) 1e2*sum(inv.logit(IFR_estimates_Sandmann2021$logit_ifr+x)*somalia_agegroups_IFR$agegroup_perc)),2) )
# paste0(round(DIC/10^floor(log10(DIC)),round_prec),"e",floor(log10(DIC)))
######
# plot
x_dodge_val=1; hpline_val=22
p <- ggplot(df_summary_plot %>% filter(seedsize>20),aes(x=factor(ifr_all_inf),group=seedsize)) + # ifr_all_inf<1.6 & 
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up,color=factor(seedsize)),position=position_dodge(width=x_dodge_val),alpha=0.3,
                 size=hpline_val,show.legend=FALSE) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up,color=factor(seedsize)),position=position_dodge(width=x_dodge_val),alpha=0.6,
                 size=hpline_val) + # ,show.legend=FALSE
  geom_hpline(aes(y=median),color="black",position=position_dodge(width=x_dodge_val),width=0.23,size=0.6) + # factor(seedsize)
  facet_wrap(~name,scales="free",labeller=labeller(ifr_all_inf=label_both),nrow=3) + # ~NPI_scale
  geom_vline(xintercept=0.5+(1:4),size=0.2,linetype="dashed") + theme_bw() + standard_theme + xlab("IFR (%)") + 
  ylab("median (CI50, CI95)") + labs(color="seedsize",caption=paste0("Burn-in: ",onefit$options$mcmc_burn_in,", Samples: ",
              onefit$options$mcmc_samples," (MCMC)")) + scale_x_discrete(expand=expansion(0,0.5)) +
  theme(axis.text.x=element_text(vjust=0.5,hjust=0.95,size=15),axis.text.y=element_text(size=15),panel.grid.major.x=element_blank(),
    legend.position="top",strip.text=element_text(size=15),legend.key.height=unit(0.8,'cm'),axis.title=element_text(size=20),
    legend.text=element_text(size=15),legend.title=element_text(size=17)) +
  guides(colour=guide_legend(override.aes=list(size=3))) + geom_text(aes(x=factor(ifr_all_inf),y=ifelse(!grepl("introd",name),NA,ci50_low-10),
        label=ifelse(!grepl("introd",name),"",introd_date_str )),size=4.5,position=position_dodge(width=x_dodge_val)); p
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"posteriors_mean_CIs_facet",ifelse(grepl("Wrap",class(p$facet)[1]),"wrap","grid"),
              "_name_NPIscale_colorcode_seedsize.png"),width=36,height=30,units="cm")
### ### ### ### ###
# plot attack rates
p <- ggplot(summ %>% filter(compartment=="attack_rate")) + 
  geom_line(aes(x=date,y=mean,color=factor(ifr_all_inf),alpha=factor(seedsize)),size=1.2) + 
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=factor(ifr_all_inf),alpha=factor(seedsize))) + 
  scale_alpha_manual(values=c(0.4,0.5,0.6,0.7,0.8)) + 
  standard_theme + theme_bw() + xlab("") + ylab("cumulative attack rate") + theme(axis.text.x=element_text(vjust=0.5,angle=90)) +
  scale_x_date(date_breaks="week",limits=as.Date(c("2020-01-10","2020-10-01")),expand=expansion(0.0)) + 
  labs(color="IFR (%)",fill="IFR (%)",alpha="seedsize") + scale_y_continuous(breaks=(0:20)/20,expand=expansion(0,0)) +
  geom_vline(data=npi_df,aes(xintercept=on),linetype="dashed",show.legend=F,size=0.5,color="darkgrey"); p
# SAVE
att_rate_filename<-paste0(parscan_mcmc_dirname,"dynamics_cumulattackrate_deaths_CDRscan.png")
if (length(p$facet$params$facets)>0) {att_rate_filename <- gsub(".png","_faceted.png",att_rate_filename)}
ggsave(att_rate_filename,width=26,height=16,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# case dynamics
ggplot(left_join(summ %>% filter(compartment %in% c("cases","death_o") & date<as.Date("2020-10-01")), fitting_params_best_estim %>% 
  select(seedsize,ifr_logit_increm,R0_fit),by=c("seedsize","ifr_logit_increm")) %>%filter(ifr_all_inf<1.6) %>%
  rename(IFR=ifr_all_inf) %>% mutate(IFR=paste0(IFR,"%"))) +
  # %>% filter(seedsize<500)
  geom_line(aes(x=date,y=mean,color=factor(seedsize), linetype=compartment),size=1.2) + # ,alpha=factor(seedsize)
  # geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=factor(seedsize),color=factor(seedsize),linetype=compartment),alpha=0.2,size=0.4) + 
  geom_point(data=fitting_incidence_modelcompare,aes(x=date,y=new_deaths),fill=NA,shape=1,size=1) +
  geom_line(data=fitting_incidence_modelcompare,aes(x=date,y=new_deaths),linetype="dashed",color="black",size=0.2) +
 facet_wrap(~IFR,scales="free",labeller=labeller(IFR=label_both)) + scale_alpha_manual(values=c(0.4,0.5,0.6,0.7,0.8)) +
  standard_theme + theme_bw() + xlab("") + ylab("cases/deaths") + 
  theme(axis.text=element_text(size=13),axis.text.x=element_text(vjust=0.5,angle=90),panel.grid.minor=element_blank(),
        axis.title=element_text(size=15),strip.text = element_text(size=15),legend.text=element_text(size=15)) +#legend.position="top",
  scale_x_date(date_breaks="3 week",expand=expansion(0.0),limits=as.Date(c("2019-09-22","2020-08-28"))) + # 
  scale_y_log10(expand=expansion(0,0),limits=c(0.6,8e3)) + labs(color="seed size",fill="seed size") + # color="IFR (%)",fill="IFR (%)"
  geom_vline(data=npi_df,aes(xintercept=on),linetype="dashed",show.legend=F,size=0.5,color="darkgrey")  # NPIs
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"dynamics_cases_deaths_withdata.png"),width=30,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# single seedsize dynamics plot
y_text=16; x_dodge_text=2; y_low=7; x_dodge_vals=c(5,24); x_start=as.Date("2020-01-10"); y_stepsize=1.1; lab_fontsize=4.5
seedsize_sel=200
ggplot(summ %>% filter(compartment=="death_o" & seedsize==seedsize_sel & ifr_all_inf<1.6) ) +
  geom_line(aes(x=date,y=mean,color=factor(round(ifr_all_inf,2)))) + 
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,fill=factor(round(ifr_all_inf,2))),alpha=0.2) + # data
  geom_line(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),linetype="dashed",size=0.3) +
  geom_point(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),fill=NA,shape=1,size=1) +
  facet_wrap(~seedsize,labeller=labeller(seedsize=label_both),scales="free",nrow=2) + # CDR~seedsize
  scale_alpha_manual(values=c(0.4,0.45,0.5,0.55)) + xlab("") + ylab("deaths (daily)") + # DIC labels
  geom_text(data=DIC_logllk_values %>% filter(seedsize==seedsize_sel & ifr_all_inf<1.6), 
        aes(x=x_start+x_dodge_text-3,y=as.numeric(factor(ifr_logit_increm))*y_stepsize+y_text,
        label="DIC=",color=factor(ifr_all_inf)),size=lab_fontsize,show.legend=FALSE) + # DIC values
  geom_text(data=DIC_logllk_values %>% filter(seedsize==seedsize_sel & ifr_all_inf<1.6),aes(x=x_start+x_dodge_text+x_dodge_vals[1],
        y=as.numeric(factor(ifr_all_inf))*y_stepsize+y_text,label=DIC_str),size=lab_fontsize,show.legend=FALSE) +
  geom_text(data=fitting_params_best_estim %>% filter(seedsize==seedsize_sel& ifr_all_inf<1.6), # introd date labels
          aes(x=x_start+x_dodge_text+4,y=as.numeric(factor(ifr_all_inf))*y_stepsize+y_text-y_low,
                label="introd. date=",color=factor(ifr_all_inf)),size=lab_fontsize,show.legend=FALSE) + 
  geom_text(data=fitting_params_best_estim %>% filter(seedsize==seedsize_sel&ifr_all_inf<1.6),aes(x=x_start+x_dodge_text+x_dodge_vals[2],
    y=as.numeric(factor(ifr_all_inf))*y_stepsize+y_text-y_low,label=gsub("-","/",introd_date)),size=lab_fontsize,show.legend=FALSE) +
  scale_x_date(date_breaks="1 week",limits=as.Date(c("2020-01-01","2020-10-01")),expand=expansion(0.0)) + 
  geom_vline(data=npi_df,aes(xintercept=on),linetype="dashed",show.legend=F,size=0.5,color="darkgrey") + labs(color="IFR",fill="IFR")+
  standard_theme + theme_bw() + labs(color="IFR",fill="IFR") + theme(axis.text=element_text(size=13),
        axis.text.x=element_text(vjust=0.5,angle=90),panel.grid.minor=element_blank(), axis.title=element_text(size=15),
        strip.text=element_text(size=15),legend.position="top",legend.text=element_text(size=15),legend.title=element_text(size=15))
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"all_dynamic_fits_NPIcolorcode_seedsize",seedsize_sel,".png"),width=36,height=22,units="cm")
