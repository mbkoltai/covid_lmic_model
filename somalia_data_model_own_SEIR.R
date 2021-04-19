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
