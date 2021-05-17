### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### multiple fits with fixed seed size and compliance levels
fitting_params <- c("R0_fit","introd_date","ifr_logit_intercept"); CDR_vals=c(baseline_daily_burials*1e4/mogadishu_popul,0.1,0.2,0.4)
parscan_mcmc_dirname="simul_output/somalia/3param_fits_seedsize_compliance_fixed/scan_seedsize_compliance_introd_date_N_120_20_simulstart1oct/"
parfit_scan_files<-list.files(parscan_mcmc_dirname); slope_val=0.0999
for (k in 1:length(parfit_scan_files)) {
  for (n_CDR in 1:length(CDR_vals)) {
  df_posteriors=readRDS(paste0(parscan_mcmc_dirname,parfit_scan_files[k]))[[n_CDR]]$posterior %>% 
    select(chain,trial,lp,all_of(fitting_params)) %>% mutate(`IFR sympt. infections (%)`=1e2*sapply(ifr_logit_intercept,
          function(x) sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc) ),
          `IFR all infections (%)`=1e2*sapply(ifr_logit_intercept,function(x) # this is bc we fitted IFR for symptom cases!!
          sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc*params$pop[[1]]$y))) %>%
    mutate(CDR=round(CDR_vals[n_CDR],3)) %>%
    mutate(NPI_scale=as.numeric(gsub("compliance|_","",str_extract(parfit_scan_files[k],"compliance.*_"))),
          seedsize=as.numeric(str_match(parfit_scan_files[k],"seedsize(.*?).rds")[2]))
  if (n_CDR==1){df_posteriors_comb=df_posteriors} else {df_posteriors_comb=rbind(df_posteriors_comb,df_posteriors)} }
  if (k==1) {df_posteriors_parscan=df_posteriors_comb} else {
    df_posteriors_parscan=bind_rows(df_posteriors_parscan,df_posteriors_comb)} }
# save
write_csv(df_posteriors_parscan,paste0(parscan_mcmc_dirname,"df_posteriors_parscan.csv"))
# calc summary stats# calc summary stats
posteriors_summary_stats = df_posteriors_parscan %>% mutate(n=row_number()) %>% 
  pivot_longer(!c(n,CDR,NPI_scale,seedsize,chain,trial,lp)) %>% group_by(name,CDR,seedsize,NPI_scale) %>% 
  summarise(mean=mean(value),median=median(value),
            ci95_low=quantile(value,probs=c(2.5,97.5)/1e2)[1],ci95_up=quantile(value,probs=c(2.5,97.5)/1e2)[2],
            ci50_low=quantile(value,probs=c(25,75)/1e2)[1],ci50_up=quantile(value,probs=c(25,75)/1e2)[2],
            post_lkl_mean=mean(lp),post_lkl_ci95_low=quantile(lp,probs=c(2.5,97.5)/1e2)[1],
            post_lkl_ci95_up=quantile(lp,probs=c(2.5,97.5)/1e2)[2]) %>%
  mutate(name=ifelse(name %in% "introd_date",paste0("introduction (days after ",
                     format(as.Date(onefit$base_parameters$date0),"%d/%m/%y"),")"),name),
         name=ifelse(name %in% "NPI_scale","NPI NPI_scale (0 to 1)",name))
# save
write_csv(posteriors_summary_stats,paste0(parscan_mcmc_dirname,"posteriors_summary_stats.csv"))

# PLOT medians and CIs for parameters
onefit=readRDS(paste0(parscan_mcmc_dirname,parfit_scan_files[k]))[[1]]; x_dodge_val=0.6; onefit$priors

### ### ### ### ### ### ### ### ### ### ### ###
# RUN simuls from param distribs of fits
for (k_seedsize in unique(df_posteriors_parscan$seedsize)) {
  for (k_npi_scale in unique(df_posteriors_parscan$NPI_scale)) {
    sel_file=parfit_scan_files[grepl(paste0("seedsize",k_seedsize,".rds"),parfit_scan_files)&
                                 grepl(paste0("compliance",k_npi_scale,"_"),parfit_scan_files)]
    # read in file
    fits_death_scaling=readRDS(paste0(parscan_mcmc_dirname,sel_file))
    n_seedsize=k_seedsize; n_compliance=k_npi_scale
    for (k in 1:length(fits_death_scaling)){ print(c(k,k_npi_scale,k_seedsize))
      dyn=cm_sample_fit(fits_death_scaling[[k]],25) %>% mutate(CDR=CDR_vals[k],NPI_scale=k_npi_scale,seedsize=k_seedsize) %>%
        filter(compartment %in% c("S","death_o"))
      if (k==1 & k_seedsize==unique(df_posteriors_parscan$seedsize)[1] & k_npi_scale==unique(df_posteriors_parscan$NPI_scale)[1]) {
        dyn_all=dyn} else {dyn_all=rbind(dyn_all,dyn)} }
  }
}

# summarize these runs
summ = dyn_all %>% group_by(t,run,compartment,CDR,NPI_scale,seedsize) %>% summarise(value=sum(value)) %>% 
  group_by(run,compartment,CDR,NPI_scale,seedsize) %>% mutate(value=ifelse(compartment %in% "S",(value[t==0]-value)/value[t==0],value),
  compartment=ifelse(compartment %in% "S","attack_rate",as.character(compartment))) %>%
  group_by(t,compartment,CDR,NPI_scale,seedsize) %>% summarise(lower=hdi(value)[[1]],upper=hdi(value)[[2]],mean=mean(value) ) %>% 
  mutate(date=as.Date(seq(as.Date(fits_death_scaling[[1]]$base_parameters$date0),
                          as.Date(fits_death_scaling[[1]]$base_parameters$date0)+max(t),1)[t+1]),CDR=round(CDR,5),
         lower=lower*(round(CDR_vals[1],5)/CDR),upper=upper*(round(CDR_vals[1],5)/CDR),mean=mean*(round(CDR_vals[1],5)/CDR))
# SAVE
write_csv(summ,paste0(parscan_mcmc_dirname,"summ.csv")); rm(dyn_all)

# compare to data # params
fitting_incidence_modelcompare = bind_rows(lapply((mogadishu_popul*CDR_vals/1e4)/baseline_daily_burials, function(x)
  data.table(out_bdr_daily_estimates %>% select(date,daily_baseline_subtr) %>% 
               mutate(new_deaths=round(daily_baseline_subtr),CDR=round(x*1e4*baseline_daily_burials/mogadishu_popul,5),
                      date_within_fitting_t=ifelse(date>=min(onefit$data$date) & date<=max(onefit$data$date),TRUE,FALSE))) )) 
#  %>% rename(new_deaths=daily_baseline_subtr)
# calculate likelihood, deviance, DIC
DIC_logllk_values <- left_join(right_join(summ %>% filter(compartment=="death_o") %>% select(!c(lower,upper)),
        subset(fitting_incidence_modelcompare, date_within_fitting_t) %>% select(!date_within_fitting_t),by=c("CDR","date")) %>% 
        mutate(logllk=ifelse(!is.infinite(dpois(new_deaths,lambda=mean,log=T)),dpois(new_deaths,lambda=mean,log=T),NA)) %>% 
        group_by(CDR,NPI_scale,seedsize) %>% summarise(sum_logllk=sum(logllk,na.rm=T),deviance=-2*sum_logllk,d_p=var(-2*logllk,na.rm=T)/2,
        DIC=deviance+d_p),fitting_incidence_modelcompare %>% group_by(CDR) %>% summarise(maxval=max(new_deaths)),by="CDR") %>%
        mutate(DIC_str=paste0(round(DIC/10^floor(log10(DIC)),1),"e",floor(log10(DIC))))

### ### ### ### ### ### ### ### ### ### ### ###
### correlations between parameters
for (k_seedsize in unique(df_posteriors_parscan$seedsize)) {
  for (k_npi_scale in unique(df_posteriors_parscan$NPI_scale)) {
    pairwise_corrs = bind_rows(lapply(unique(df_posteriors_parscan$CDR), function(x) cbind(df_posteriors_parscan %>% 
        filter(CDR==x,seedsize==k_seedsize,NPI_scale==k_npi_scale) %>% select(all_of(fitting_params)) %>% corrr::correlate() %>% 
        corrr::shave() %>% corrr::stretch(),CDR=x,seedsize=k_seedsize,NPI_scale=k_npi_scale))) %>% 
 mutate(x=factor(x,levels=unique(x)),y=factor(y,levels=unique(y)),highlight_color=ifelse(abs(r)<0.8|is.na(r),0,1)) %>% filter(!is.na(r))
    if (k_seedsize==unique(df_posteriors_parscan$seedsize)[1] & k_npi_scale==unique(df_posteriors_parscan$NPI_scale)[1]) {
      pairwise_corrs_parscan = pairwise_corrs} else { pairwise_corrs_parscan = bind_rows(pairwise_corrs_parscan,pairwise_corrs)}
  }
}

# plot pairwise correlations
ggplot(pairwise_corrs_parscan,aes(x,y,fill=r))+geom_tile(color="black") + # aes(color=factor(highlight_color)),size=2,width=0.96,height=0.96
  facet_grid(CDR~seedsize+NPI_scale,labeller=labeller(CDR=label_both,seedsize=label_both,NPI_scale=label_both)) + 
  geom_text(aes(label=round(r,2)),size=5) + scale_x_discrete(expand=c(0,0.02)) + scale_y_discrete(expand=c(0,0.02)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",limit=c(-1,1)) + # scale_color_manual(values=c("grey","black"),guide=FALSE) +
  standard_theme + theme_bw() + xlab("") + ylab("") +
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.95),
        axis.ticks = element_blank())
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"/param_correlations.png"),width=36,height=16,units="cm")

# scatterplots of fitted parameters: (R0_fit,ifr_logit_intercept) (R0_fit,introd_date) (introd_date,ifr_logit_intercept)
var_pairs <- list(c("R0_fit","ifr_logit_intercept"),c("R0_fit","introd_date"),c("introd_date","ifr_logit_intercept"))
for (k in 1:length(var_pairs)){
p <- ggplot(df_posteriors_parscan,aes(x=get(var_pairs[[k]][1]),y=get(var_pairs[[k]][2]))) + 
  geom_point(aes(color=factor(CDR))) + geom_smooth(aes(group=CDR),method="lm",size=0.5,color="black") + # 
  facet_wrap(NPI_scale~seedsize,labeller=labeller(CDR=label_both,seedsize=label_both,NPI_scale=label_both),scales="free") +
  theme_bw() + standard_theme + labs(color="CDR") + xlab(var_pairs[[k]][1]) + ylab(var_pairs[[k]][2]); p 
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"mcmc_scatterplot_",var_pairs[[k]][1],"_",var_pairs[[k]][2],".png"),width=30,height=18,units="cm")
}

# plots chains of mcmc wrt post log-llh and param values
# df_posteriors_parscan$
for (k in 3:6){
p <- ggplot(df_posteriors_parscan,aes(x=trial,y=get(colnames(df_posteriors_parscan)[k]))) + geom_line(aes(color=factor(CDR)),size=0.3) +
  facet_wrap(NPI_scale~seedsize,labeller=labeller(CDR=label_both,seedsize=label_both,NPI_scale=label_both),scales="free") +
  theme_bw() + standard_theme + labs(color="CDR") + ylab(colnames(df_posteriors_parscan)[k]); p
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"mcmc_",colnames(df_posteriors_parscan)[k],"_by_chain.png"),width=30,height=18,units="cm")
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT DYNAMIC fits faceted by (CDR-seedsize-NPI_scale)
fitting_dates=fits_death_scaling[[1]]$data %>% summarise(min_date=min(date),max_date=max(date))
ggplot(summ %>% filter(compartment=="death_o")) + 
  geom_line(aes(x=date,y=mean), colour="blue") + geom_ribbon(aes(x=date,ymin=lower,ymax=upper), fill="blue",alpha=0.2) + 
  geom_line(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),linetype="dashed",size=0.3) +
  geom_point(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),fill=NA,shape=1,size=0.2) +
  facet_grid(CDR~seedsize+NPI_scale,labeller=labeller(CDR=label_both,seedsize=label_both,NPI_scale=label_both),scales="free") +
  geom_rect(data=subset(DIC_logllk_values,DIC<1.1*min(DIC)),xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill=NA,color="red",size=1.2) +
  xlab("") + ylab("deaths (daily)") + standard_theme + theme_bw() + theme(axis.text.x=element_text(vjust=0.5,angle=90)) + 
  geom_text(data=DIC_logllk_values,aes(x=fitting_dates$min_date+70,y=1.4*maxval, # left_join(,summ,by=c("CDR","NPI_scale","seedsize"))
                                       label=paste0("DIC=",round(DIC))),color="black",size=3.5)+
  scale_x_date(date_breaks="4 week",limits=as.Date(c("2020-01-01","2020-10-15")),expand=expansion(0.0)) + 
  scale_y_continuous(breaks=(0:10)*2,expand=expansion(0,0))
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"all_dynamic_fits.png"),width=38,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### 
# fits with DIC values
fitting_params_best_estim = left_join(posteriors_summary_stats %>% mutate(CDR=ifelse(CDR==0.042,round(CDR_vals[1],5),CDR)) %>% 
        filter(name %in% c("IFR all infections (%)","R0_fit","introduction (days after 01/10/19)")) %>% 
        select(name,CDR,seedsize,NPI_scale,mean),DIC_logllk_values %>% select(!c(deviance,d_p)),by=c("CDR","seedsize","NPI_scale")) %>%
  pivot_wider(names_from=name,values_from=mean) %>% rename(IFR=`IFR all infections (%)`) %>%
  mutate(introd_date=as.Date(onefit$base_parameters$date0)+`introduction (days after 01/10/19)`)

# faceted by (CDR-seedsize), NPI_scale by color
ggplot(summ %>% filter(compartment=="death_o")) + 
  geom_line(aes(x=date,y=mean,group=interaction(NPI_scale,CDR),color=factor(NPI_scale))) + # ,alpha=factor(CDR)
  geom_ribbon(aes(x=date,ymin=lower,ymax=upper,group=interaction(NPI_scale,CDR),fill=factor(NPI_scale),alpha=factor(CDR))) + #alpha=0.2
  geom_line(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),linetype="dashed",size=0.3) +
  geom_point(data=fitting_incidence_modelcompare %>% mutate(t=as.numeric(date)),aes(x=date,y=new_deaths),fill=NA,shape=1,size=0.2) +
  facet_wrap(~seedsize,labeller=labeller(CDR=label_both,seedsize=label_both),scales="free",nrow=2) + # CDR~seedsize
  # geom_rect(data=fitting_dates,aes(xmin=min_date,xmax=max_date,ymin=-Inf,ymax=Inf),fill=NA,linetype="dashed",color="red",size=0.4) +
  scale_alpha_manual(values=c(0.4,0.45,0.5,0.55)) + xlab("") + ylab("deaths (daily)") + standard_theme + theme_bw() +
  theme(axis.text.x=element_text(vjust=0.5,angle=90),panel.grid.minor=element_blank()) +
  geom_text(data=DIC_logllk_values,aes(x=fitting_dates$min_date-4.5,y=as.numeric(factor(NPI_scale))*2.05+15,label="DIC=",
                                       color=factor(NPI_scale)),size=3.5) + # DIC values
  geom_text(data=DIC_logllk_values,aes(x=fitting_dates$min_date-6+as.numeric(factor(CDR))*19,y=as.numeric(factor(NPI_scale))*2.05+15,
                                       label=paste0(DIC_str,ifelse(as.numeric(factor(CDR))!=4,",",""))),size=3.5) + # ,color=factor(NPI_scale)
  geom_text(data=fitting_params_best_estim,aes(x=fitting_dates$max_date-68.5,y=as.numeric(factor(NPI_scale))*2.1+15,
                                               label="introd.=",color=factor(NPI_scale)),size=3.25) + # introd dates
  geom_text(data=fitting_params_best_estim,aes(x=fitting_dates$max_date-67+as.numeric(factor(CDR))*17.5,
                                               y=as.numeric(factor(NPI_scale))*2.1+15,
                                               label=paste0(gsub("-","/",gsub("2020-|2019-","",introd_date)),ifelse(as.numeric(factor(CDR))!=4,",",""))),
            size=3.25)+ scale_x_date(date_breaks="2 week",limits=as.Date(c("2020-01-01","2020-10-17")),expand=expansion(0.0)) + 
  scale_y_continuous(breaks=(0:12)*2,limits=c(0,22.5),expand=expansion(0,0)) + 
  labs(color="NPI scaling",fill="NPI scaling",alpha="CDR scaling")
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"all_dynamic_fits_NPIcolorcode.png"),width=36,height=18,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PLOT DIC as function of CDR~seedsize~NPI_scale
ggplot(fitting_params_best_estim %>% mutate(R0_label=paste0("R0=",round(R0_fit,2))) ) + 
  geom_hpline(aes(x=factor(NPI_scale),y=DIC,group=CDR,color=factor(seedsize),alpha=factor(CDR)),width=0.22,
              position=position_dodge(width=1)) + scale_alpha_manual(values=c(0.4,0.55,0.7,0.85)) + 
  geom_vline(xintercept=(2:3)-0.5,size=0.2,linetype="dashed") + # facet_wrap(~seedsize,labeller=labeller(CDR=label_both),nrow = 3) +
  scale_x_discrete(expand=expansion(0,0.1)) + theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5,hjust=0.95)) +
  geom_text(aes(x=factor(NPI_scale),y=ifelse(as.numeric(factor(seedsize)) %% 2==0,DIC-20,DIC+20),# ifelse(NPI_scale==0.5,-1,1)
      group=CDR,label=paste0(gsub("2019/|2020/","",gsub("-","/",introd_date)),", ",R0_label)),position=position_dodge(width=1),size=3,
      check_overlap=T) + xlab("NPI scaling") + ylab("DIC") + labs(color="seedsize",alpha="CDR") 
# caption="Introd. date not shown when equal to 2019-11-11"
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"DIC_xaxis_NPIscaling_colorcode_seedsize_alpha_CDR.png"),width=36,height=18,units="cm")

# PLOT facet_grid NAME-SEEDSIZE, x-axis CDR, color coded by NPIscale
x_dodge_val=1
p<-ggplot(left_join(posteriors_summary_stats %>% group_by(name) %>% mutate(median_range=max(median)-min(median)),
    DIC_logllk_values %>% mutate(CDR=round(CDR,3)) %>% select(!c(maxval,deviance,d_p)),by=c("CDR","seedsize","NPI_scale")) %>% 
    filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")), aes(x=factor(CDR),group=interaction(NPI_scale,CDR))) +
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up,color=factor(NPI_scale)),position=position_dodge2(width=x_dodge_val),alpha=0.3,size=8) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up,color=factor(NPI_scale)),position=position_dodge2(width=x_dodge_val),alpha=0.6,size=8) +
  geom_hpline(aes(y=median),color="black",position=position_dodge2(width=x_dodge_val),width=1/3,size=0.6) +
  facet_grid(name~seedsize,scales="free",labeller=labeller(CDR=label_both,seedsize=label_both)) + 
  geom_vline(xintercept=0.5+1:3,size=0.2,linetype="dashed") + theme_bw() + standard_theme + xlab("CDR scaling") + 
  ylab("mean (CI50, CI95)") + labs(color="NPI scaling",caption=paste0("Burn-in: ",onefit$options$mcmc_burn_in,", Samples: ",
      onefit$options$mcmc_samples," (MCMC). Labels in plot are DIC values.")) + scale_x_discrete(expand=expansion(0,0.5)) +
  theme(axis.text.x=element_text(vjust=0.5,hjust=0.95),panel.grid.major.x=element_blank(),legend.position="top") +
  geom_text(aes(x=factor(CDR),group=interaction(NPI_scale,CDR),y=median+0.08*median_range,label=DIC_str),
            size=2.5,position=position_dodge2(width=x_dodge_val)); p
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"posteriors_mean_CIs_facet",ifelse(grepl("Wrap",class(p$facet)[1]),"wrap","grid"),
  "_name_seedsize_colorcodeNPI",ifelse(length(p$layers)>4,"_geomtext",""),".png"),width=36,height=18,units="cm")

### ###
# PLOT facet_grid NAME-NPI, x-axis CDR, color coded by SEEDSIZE
x_dodge_val=1; hpline_val=8
p<-ggplot(left_join(posteriors_summary_stats %>% group_by(name) %>% mutate(median_range=max(median)-min(median)),
        DIC_logllk_values %>% mutate(CDR=round(CDR,3)) %>% select(!c(maxval,deviance,d_p)),by=c("CDR","seedsize","NPI_scale")) %>% 
            filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")) %>% 
            mutate(DIC_str=paste0(round(DIC/10^floor(log10(DIC)),1),"e",floor(log10(DIC)))),
          aes(x=factor(CDR),group=interaction(seedsize,CDR))) +
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up,color=factor(seedsize)),position=position_dodge2(width=x_dodge_val),alpha=0.3,size=hpline_val) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up,color=factor(seedsize)),position=position_dodge2(width=x_dodge_val),alpha=0.6,size=hpline_val) +
  geom_hpline(aes(y=median),color="black",position=position_dodge2(width=x_dodge_val),width=0.24,size=0.6) +
  facet_grid(name~NPI_scale,scales="free",labeller=labeller(CDR=label_both,NPI_scale=label_both)) + 
  geom_vline(xintercept=0.5+1:3,size=0.2,linetype="dashed") + theme_bw() + standard_theme + xlab("CDR scaling") + 
  ylab("mean (CI50, CI95)") + labs(color="seedsize",caption=paste0("Burn-in: ",onefit$options$mcmc_burn_in,", Samples: ",
    onefit$options$mcmc_samples," (MCMC). Labels in plot are DIC values.")) + scale_x_discrete(expand=expansion(0,0.5)) +
  theme(axis.text.x=element_text(vjust=0.5,hjust=0.95),panel.grid.major.x=element_blank()) +
  geom_text(aes(x=factor(CDR),group=interaction(seedsize,CDR),y=median+0.07*median_range,label=DIC_str),
            size=3,position=position_dodge2(width=x_dodge_val)); p # scale_y_continuous(expand=expansion(0.25,0))
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"posteriors_mean_CIs_facet",ifelse(grepl("Wrap",class(p$facet)[1]),"wrap","grid"),
  "_name_NPIscale_colorcode_seedsize",ifelse(length(p$layers)>4,"_geomtext",""),".png"),width=36,height=18,units="cm")

### ###
# PLOT facet_grid NAME-SEEDSIZE, x-axis NPI, color coded by CDR
x_dodge_val=1
p<-ggplot(left_join(posteriors_summary_stats %>% group_by(name) %>% mutate(median_range=max(median)-min(median)),
        DIC_logllk_values %>% mutate(CDR=round(CDR,3)) %>% select(!c(maxval,deviance,d_p)),by=c("CDR","seedsize","NPI_scale")) %>% 
         filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")) %>% 
         mutate(DIC_str=paste0(round(DIC/10^floor(log10(DIC)),1),"e",floor(log10(DIC)))),
       aes(x=factor(NPI_scale),group=interaction(CDR,NPI_scale))) +
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up,color=factor(CDR)),position=position_dodge2(width=x_dodge_val),alpha=0.3,size=8) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up,color=factor(CDR)),position=position_dodge2(width=x_dodge_val),alpha=0.6,size=8) +
  geom_hpline(aes(y=median),color="black",position=position_dodge2(width=x_dodge_val),width=0.24,size=0.6) +
  facet_grid(name~seedsize,scales="free",labeller=labeller(CDR=label_both,seedsize=label_both)) + 
  geom_vline(xintercept=c(1.5,2.5),size=0.2,linetype="dashed") + 
  theme_bw() + standard_theme + xlab("NPI scaling") + ylab("mean (CI50, CI95)") + # scale_y_continuous(expand=expansion(0.25,0)) + 
  labs(color="CDR scaling",caption=paste0("Burn-in: ",onefit$options$mcmc_burn_in,", Samples: ",onefit$options$mcmc_samples,
                  " (MCMC). Labels in plot are posterior likelihoods")) + scale_x_discrete(expand=expansion(0,0.5)) +
  theme(axis.text.x=element_text(vjust=0.5,hjust=0.95),panel.grid.major.x=element_blank()) +
  geom_text(aes(x=factor(NPI_scale),group=interaction(CDR,NPI_scale),y=median+0.05*median_range,label=DIC_str),
            size=2.5,position=position_dodge2(width=x_dodge_val)); p
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"posteriors_mean_CIs_facet",ifelse(grepl("Wrap",class(p$facet)[1]),"wrap","grid"),
      "_name_seedsize_colorcodeCDR",ifelse(length(p$layers)>4,"_geomtext",""),".png"),width=36,height=18,units="cm")

### ###
# PLOT facet_grid NAME-NPIscaling, x-axis seed size, color coded by CDR
x_dodge_val=1
p<-ggplot(left_join(posteriors_summary_stats %>% group_by(name) %>% mutate(median_range=max(median)-min(median)),
    DIC_logllk_values %>% mutate(CDR=round(CDR,3)) %>% select(!c(maxval,deviance,d_p)),by=c("CDR","seedsize","NPI_scale")) %>% 
    filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")) %>% 
      mutate(DIC_str=paste0(round(DIC/10^floor(log10(DIC)),1),"e",floor(log10(DIC)))),
    aes(x=factor(seedsize),group=interaction(CDR,seedsize))) +
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up,color=factor(CDR)),position=position_dodge2(width=x_dodge_val),alpha=0.3,size=8) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up,color=factor(CDR)),position=position_dodge2(width=x_dodge_val),alpha=0.6,size=8) +
  geom_hpline(aes(y=median),color="black",position=position_dodge2(width=x_dodge_val),width=0.23,size=0.6) +
  facet_grid(name~NPI_scale,scales="free",labeller=labeller(CDR=label_both,NPI_scale=label_both)) + 
  geom_vline(xintercept=(1:3)+0.5,size=0.2,linetype="dashed") + 
  theme_bw() + standard_theme + xlab("seed size") + ylab("mean (CI50, CI95)") + # scale_y_continuous(expand=expansion(0.25,0)) + 
  labs(color="CDR scaling",caption=paste0("Burn-in: ",onefit$options$mcmc_burn_in,", Samples: ",onefit$options$mcmc_samples,
            " (MCMC). Labels are DIC values.")) + scale_x_discrete(expand=expansion(0,0.5)) + 
  theme(axis.text.x=element_text(vjust=0.5,hjust=0.95),panel.grid.major.x=element_blank()) +
  geom_text(aes(x=factor(seedsize),group=interaction(CDR,seedsize),y=median+0.05*median_range,label=DIC_str),
            size=2.5,position=position_dodge2(width=x_dodge_val)); p
# SAVE
ggsave(paste0(parscan_mcmc_dirname,"posteriors_mean_CIs_facet",ifelse(grepl("Wrap",class(p$facet)[1]),"wrap","grid"),
              "_name_NPIscale_colorcodeCDR",ifelse(length(p$layers)>4,"_geomtext",""),".png"),width=36,height=22,units="cm")

# ### ### ### ### ### ### ### ###
# # plot log-posterior by NPIscaling-seedsize-CDR
# ggplot(posteriors_summary_stats %>% filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")) %>%
#   group_by(CDR,NPI_scale,seedsize) %>% summarise(post_lkl_mean=unique(post_lkl_mean),post_lkl_ci95_low=unique(post_lkl_ci95_low),
#    post_lkl_ci95_up=unique(post_lkl_ci95_up)),aes(x=factor(NPI_scale),group=interaction(NPI_scale,seedsize),color=factor(seedsize))) +
# geom_hpline(aes(y=post_lkl_mean,color=factor(seedsize)),width=0.98) +facet_wrap(~CDR,scales="free",labeller=labeller(CDR=label_both))+
# geom_vline(xintercept=(1:4)-0.5,size=0.2,linetype="dashed") + theme_bw() + standard_theme +xlab("NPI scaling")+ylab("log-posterior") +
#   scale_y_continuous(expand=expansion(0.1,0),trans=pseudolog10_trans) + # scale_y_continuous(expand=expansion(0.05,0)) +
# labs(color="seed size",caption=paste0("Burn-in: ",onefit$options$mcmc_burn_in,", Samples: ",onefit$options$mcmc_samples," (MCMC)")) +
#   theme(axis.text.x=element_text(vjust=0.5,hjust=0.95)) # scale_x_discrete(expand=expansion(0.06,0)) + 
# # SAVE
# ggsave(paste0("simul_output/somalia/scan_seedsize_compliance/posterior_likelihoods.png"),width=30,height=22,units="cm")


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# extract results of multiple mcmc fits (different priors / sample size)
# fit_dirs=list.dirs("simul_output/somalia",recursive=F)[grepl("fit_",list.dirs("simul_output/somalia",recursive=F))]
# l_pairwise_corrs=list()
# #######
# # start loop
# for (k_file in 1:length(fit_dirs)){
#   mcmc_filename=list.files(fit_dirs[k_file],pattern=".rds",full.names=T); fits_death_scaling <- readRDS(mcmc_filename)
#   foldertag=gsub("mcmc_","",
#             paste0(paste0(paste0(c(min(fits_death_scaling[[1]]$data$date),max(fits_death_scaling[[1]]$data$date)),collapse="_"),"_"),
#                   paste0("seedsize_",gsub(" ","_",fits_death_scaling[[1]]$priors$seed_size),"_"),
#                       paste0(names(fits_death_scaling[[1]]$options[c("mcmc_burn_in","mcmc_samples")]),
#                                as.numeric(fits_death_scaling[[1]]$options[c("mcmc_burn_in","mcmc_samples")]),collapse = "_") ) )
#   mcmc_foldername=paste0("simul_output/somalia/fit_",foldertag)
#   if (!dir.exists(mcmc_foldername)) {dir.create(mcmc_foldername); file.copy(mcmc_filename,mcmc_foldername)}
#   # CDR_vals=c(baseline_daily_burials*1e4/mogadishu_popul,0.1,0.2,0.4); slope_val=linregr$coefficients["agegroup_mean"]
#   # fitting_date_window
#   for (k in 1:length(CDR_vals)) { df_posteriors=fits_death_scaling[[k]]$posterior %>% select(chain,trial,lp,all_of(fitting_params)) %>%
#     mutate(`IFR sympt. infections (%)`=1e2*sapply(ifr_logit_intercept,
#                   function(x) sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc) ),
#           `IFR all infections (%)`=1e2*sapply(ifr_logit_intercept,function(x) 
#           sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc*params$pop[[1]]$y))) %>%
#     mutate(CDR=round(CDR_vals[k],3))
#   if (k==1){df_posteriors_comb=df_posteriors} else {df_posteriors_comb=rbind(df_posteriors_comb,df_posteriors)} }
#   
#   ### correlations between parameters ----------------
#   # scatterplots + corrs for all values of CDR
#   corr_scatterplot=0; if (corr_scatterplot){
#     png(filename=paste0(mcmc_foldername,"/param_corrs_CDR_scatterplot.png"),units="cm",width=30,height=24,res=500)
#     p <- ggpairs(df_posteriors_comb %>% select(all_of(fitting_params),CDR) %>% mutate(CDR=factor(CDR)),aes(color=CDR,alpha=0.4)) +
#       theme_bw() + standard_theme + labs(caption=foldertag); print(p); dev.off()}
#   # SAVE # ggsave(paste0(mcmc_foldername,"/param_corrs_CDR_scatterplot.png"),width=30,height=24,units="cm")
#   # correlations only
#   pairwise_corrs = bind_rows(lapply(unique(df_posteriors_comb$CDR), function(x) cbind(df_posteriors_comb %>% filter(CDR==x) %>%
#     select(all_of(fitting_params)) %>% corrr::correlate() %>% corrr::shave() %>% corrr::stretch() ,CDR=x))) %>%
# mutate(x=factor(x,levels=unique(x)),y=factor(y,levels=unique(y)),highlight_color=ifelse(abs(r)<0.7|is.na(r),0,1)) %>% filter(!is.na(r))
#   l_pairwise_corrs[[k_file]]= pairwise_corrs %>% mutate(folder_tag=gsub("2020-","",gsub("_burn.*","",foldertag)))
#   # plot
#   corr_scalar_plot=0
#   if (corr_scalar_plot){
#     ggplot(pairwise_corrs,aes(x,y,fill=r)) + geom_tile(aes(color=factor(highlight_color)),size=2,width=0.96,height=0.96) +
#       facet_wrap(~CDR,labeller=labeller(CDR=label_both)) + geom_text(aes(label=round(r,2)),size=5) +
#       scale_x_discrete(expand=c(0,0.02)) + scale_y_discrete(expand=c(0,0.02)) +
#       scale_fill_gradient2(low="blue",mid="white",high="red",limit=c(-1,1)) + scale_color_manual(values=c(NA,"black"),guide=FALSE) +
#       standard_theme + theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +xlab("") + ylab("")
#     # ## SAVE
#     ggsave(paste0(mcmc_foldername,"/param_corrs_CDR.png"),width=30,height=24,units="cm")
#   }
# } # end of loop
# 
# # plot for all fits
# ggplot(bind_rows(l_pairwise_corrs),aes(x,y,fill=r)) + geom_tile(aes(color=factor(highlight_color)),size=1,width=0.96,height=0.96) +
#   facet_grid(CDR~folder_tag,labeller=labeller(CDR=label_both)) + geom_text(aes(label=round(r,2)),size=4) +
#   scale_x_discrete(expand=c(0,0.02)) + scale_y_discrete(expand=c(0,0.02)) + scale_color_manual(values=c("grey","black"),guide=FALSE) +
#   scale_fill_gradient2(low="blue",mid="white",high="red",limit=c(-1,1)) + standard_theme + theme_bw() + 
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.95)) + 
#   xlab("") + ylab("")
# ggsave("simul_output/somalia/all_fits_param_corrs_CDR.png",width=40,height=20,units="cm")