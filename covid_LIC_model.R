# simulation study to compare potential mechanisms explaining the lower than expected case and death numbers in LIC

# to render as html: rmarkdown::render("RSV_kenya_SA_calculation.R",output_dir='output/cea_plots/')
# library(tidyverse); library(reshape2); library(matrixStats); library(rstudioapi); # library(fitdistrplus); library(qs); library(data.table)
rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","countrycode","coronavirus","RcppRoll","scales","dttr2","wpp2019"),
       library,character.only=TRUE)
cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
# functions and plotting theme
source("covid_LIC_fcns.R")
### ### ### ### ### ### ### ### ### ### ### ###
### SEIR age structured (covidm) --------------------------

# SEIR age-structured model (own) --------------------------
seir_agestr_ode <- function(t,X,parms){ # params=list(K_m,C_m,b_m_full,u_val,a_m,ind_all_suscept,ind_all_inf_vars)
  K_m=parms[[1]];C_m=parms[[2]];b_m_full=parms[[3]]; u_val=parms[[4]];a_m=parms[[5]]; 
  susc_ind=parms[[6]]; inf_vars_inds=parms[[7]]
  dXdt=b_m_full %*% diag(X[susc_ind]) %*% C_m %*% diag(u_val) %*% a_m %*% X[inf_vars_inds] + K_m%*%X; list(dXdt) }
# SYSTEM DIMENSIONS
n_age=16; vartype_list=c('S','E','Ip','Ic','Is'); infect_vartype=c('E','Ip','Ic','Is'); dim_sys=n_age*length(vartype_list)
# full list of vars
full_varname_list=fun_seir_agestr_varnames(vartype_list,n_age)
# name+age --> linear index: fun_sub2ind_seir_agestr(j_age=2,varname='Is',vartype_list,length(vartype_list),n_age=2)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### POPULATION (age struct resol as in covidm) --------------------------
age_groups <- data.frame(age_group=c(1:16), age_low=c(seq(0,75,5) ), age_high=c(seq(4,74,5),100))
# population data from wpp2019
countryval="Sudan"; N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups) # uniform: N_tot=rep(1e6,n_age)
# CONTACT MATRIX
if (!exists("covid_params")){cm_force_rebuild=F; cm_build_verbose=T; cm_version=2; source(file.path(cm_path,"R","covidm.R"))}
covid_params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval))
C_m_orig=Reduce('+',covid_params$pop[[1]]$matrices)
# Make matrix reciprocal
C_m=fun_contmatr_recipr(C_m_orig,N_tot)
# i=7;j=4; c(C_m[i,j]*N_tot[i],C_m[j,i]*N_tot[j])
# constant transmission parameters
d_e=1/3; d_p=1/2; d_c=1/3; d_s=1/5; infect_first_ord_pars=c(d_e,d_p,d_c,d_s)
########################################################
### single simulation with COVIDM ---------------------------
params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval))
target_R0=2.4; scale_r0=target_R0/cm_calc_R0(params,1)
# scale susceptibility
for (k in 1:length(params$pop)){ params$pop[[k]]$u=params$pop[[k]]$u*scale_r0}
# run the model
ptm<-proc.time(); run=cm_simulate(params,1); proc.time()-ptm; lmic_simul=run$dynamics
# simul output wrangling
lmic_simul=fcn_covidm_df(lmic_simul,sel_vars=c("S","E","Ip","Ic","Is","R"),params)
### PLOT timecourse ---
sel_metric=unique(lmic_simul$variable)[3]
ggplot(subset(lmic_simul, compartment %in% "cases" & variable %in% sel_metric & t>=ifelse(grepl("percent_cases",sel_metric),10,0))) + 
 geom_line(aes(x=t,y=ifelse(grepl("percent",sel_metric),1e2,1)*value,colour=group,group=group),size=1.1) + 
 geom_hline(aes(yintercept=pop_perc*1e2),linetype="dashed") + facet_wrap(~group,scales="free",ncol=3) +
 theme_bw() + standard_theme + ggtitle(paste0(countryval," epidemic timecourse (",sel_metric,")")) +
 ylab(paste0(sel_var," (",ifelse(grepl("percent",sel_metric),"% ",""),gsub("percent_","",sel_metric),")")) 
ggsave(paste0("simul_output/timecourse_",sel_metric,".png"),height=28,width=20,units = "cm")
# + xlim(c(0,120))
# save

### PLOT age distribution of cases ---
sel_metric=unique(lmic_simul$variable); sel_compartm="cumul_sympt_cases" # R cumul_sympt_cases
sel_time=100+(0:6)*20; # if (grepl("perc",sel_metric)) {scale_f=1e2} else {scale_f=1}
ggplot(subset(lmic_simul, compartment %in% sel_compartm & t %in% sel_time)) + #  (lmic_simul$population %in% sel_reg) &
  geom_line(aes(x=group,y=value,colour=t,group=t),size=1.1) + facet_wrap(~variable,scales="free") + 
  scale_color_gradientn(colours=wes_palette("Zissou1")) + theme_bw() + standard_theme + 
  xlab("age group") + ylab("cumulative cases") + ggtitle("Cumulative clinical attack rate by age group")
ggsave(paste0("simul_output/",countryval,"_cases_cumul_agedistrib.png"),width=30,height=12,units="cm")
########################################################
### AGE DEPENDENT parameters ---------------------------
min_val_susc=1e-2; maxval_susc=0.15; midpoint_susc=(min_val_susc+maxval_susc)/2; delta_susc=maxval_susc-midpoint_susc
midpoint_clinfract=0.4; delta_clinfr=0.35; rep_min_susc=3; rep_max_susc=5; rep_min_clinfr=5; rep_max_clinfr=3
### SINGLE SIMULATION ---------------------------
# set params
depval=1/2
u_val=fun_lin_approx_agedep_par(min_val=midpoint_susc-depval*delta_susc,max_val=midpoint_susc+depval*delta_susc,
                                rep_min_susc,rep_max_susc)
y_val=fun_lin_approx_agedep_par(midpoint_clinfract-depval*delta_clinfr,midpoint_clinfract+depval*delta_clinfr,
                                rep_min_clinfr,rep_max_clinfr)
f_val=rep(0.5,n_age)
# KINETIC MATRIX (LINEAR TERMS)
K_m=fun_seir_agestr_kinmatr(infect_vartype,vartype_list,n_age,infect_first_ord_pars,y_val)
# vectors to scale FORCE of INFECTION terms
list_bm_am=fun_force_inf_vects(vartype_list,forceinf_vars=c("S","E"),n_age,f_val,N_tot); b_m_full=list_bm_am[[1]];a_m=list_bm_am[[2]]
# inf vector is: b_m_full %*% diag(c(1,2)) %*% C_m %*% diag(u_val) %*% a_m %*% matrix(1,8,1)
l=fun_inds_vartypes(n_age,vartype_list,infect_vartype); ind_all_inf_vars=l[[1]]; ind_all_suscept=l[[2]]
# INITIAL CONDITIONS (seed epidemic by "E">0 in a given (or multiple) age groups)
initvals_seir_model=fcn_set_init_conds(inf_initval=10,init_inf_age_groups=7,init_inf_var="E",
                                       n_age,N_tot,vartype_list,ind_all_suscept)
### run ODEs ----------------------------
# time resolution and duration
n_days_year=365; n_years=1; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=0.25)
# define inputs for ODE fcn
params_ode=list(K_m, C_m, b_m_full, u_val, a_m, ind_all_suscept, ind_all_inf_vars)
# RUN simulation
ptm<-proc.time(); ode_solution<-lsoda(initvals_seir_model,timesteps,func=seir_agestr_ode,parms=params_ode); proc.time() - ptm
# PROCESS OUTPUT
# df_ode_solution = ode_solution %>% as.data.frame() %>% setNames(c("t",full_varname_list))
l_proc_sol=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars); df_ode_solution=l_proc_sol[[1]]
# final values
final_vals=fun_final_vals(df_ode_solution,N_tot,"fract"); sum((1-final_vals)*N_tot)/sum(N_tot)
# PLOT DYNAMICS (df_ode_solution_tidy=l_proc_sol[[2]])
fun_seir_agegroups_dyn(age_groups,l_proc_sol[[2]],xlim_val=220,abs_or_fract=2,timesteps,standard_theme)
ggsave(paste0("simul_output/",countryval,"_simul/SEIR_age_str_",n_age,"agegroups_",savetag,".png"),width=30,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
# average age: fun_meanage_cntr(popM,popF,"Pakistan",year_str="2020")
### PARAMETER SCAN ---------------------------
cntr_list=c("Nigeria","Sudan","Egypt","South Africa","Brazil","Thailand","Russian Federation","France","Italy")
scan_param_fullname=c("susceptibility","clinical fraction","asymptom. infectiousness")
### Attack rates ---------------------------
# Scan in one parameter & one country ---------------------------
# run scan in 1 param with own model
# age_dep_paramtable=fun_agedep_partable(agedep_param_lims,k_max=6)
# ptm<-proc.time(); attack_rates_singlecntr=fun_paramscan_singlepar(age_dep_paramtable,scan_param_fullname[3],infect_first_ord_pars,
#           infect_vartype,vartype_list,C_m,inf_initval=10,init_inf_age_groups=7,init_inf_var="E",N_tot,n_years,age_groups); 
# proc.time() - ptm
### with COVIDM
ptm<-proc.time(); attack_rates_singlecntr_CM=fun_paramscan_single_covidm(countryval="Egypt",age_dep_paramtable,scan_parameter="susc",
                                    sel_vars=c("S","E","Ip","Ic","Is","R"),n_years=1,age_groups); proc.time() - ptm
# attack_rates_singlecntr[,"popul"]=rep(N_tot,length(unique(attack_rates_singlecntr$col_scale)))
# sum_attackrate=attack_rates_singlecntr %>% group_by(col_scale) %>% summarise(sum_clin_attrate=sum(clinical_attackrate*popul)/sum(popul))
# Scan in multiple params and cntrs ---------------------------
agedep_param_lims=c(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,delta_clinfr,
                    rep_min_clinfr,rep_max_clinfr)
# with own script
# attack_rates_multipar_multicntr=fun_onedim_multipar_multicntr_scan(cntr_list,year_str="2020",age_groups,cm_path,scan_param_fullname,
#           n_loop=6,agedep_param_lims,
#           infect_first_ord_pars,infect_vartype,vartype_list,inf_initval=10,init_inf_age_groups=7,init_inf_var="E",n_years=1)
# with COVIDM
attack_rates_multipar_multicntr_CM=fun_onedim_multipar_multicntr_scan_covidm(cntr_list,year_str="2020",age_groups,
                            scan_param_fullname,n_loop=6, agedep_param_lims,sel_vars=c("S","E","Ip","Ic","Is","R"),n_years=1)
# PLOT attack rates by age groups
sel_var=c("fract_agegroup","fraction_cases","per_1000_popul")[2]; truthval_var=grepl("clin",attack_rates_multipar_multicntr$name)
ggplot(attack_rates_multipar_multicntr_CM[truthval_var,], aes_string(x="agegroup_names",y=sel_var,group="col_scale",color="col_scale")) + 
 geom_line() + geom_point() + facet_grid(country~scanpar,switch="y") + theme_bw() + standard_theme + 
 xlab("age group") + labs(color="age-dependence/level \nof transmission \nparameter") +
 ggtitle(paste0("clinical attack rates (",gsub("_"," ",gsub("value","fraction",sel_var)),")"))
# scale_y_continuous(breaks=seq(0,1,0.05))
# SAVE
ggsave(paste0("simul_output/attackrate_",sel_var,"_agecolored.png"),width=30,height=25,units="cm")
#######
# summed attackrates (across age groups)
attackrate_all_ages=attack_rates_multipar_multicntr %>% group_by(par_scale,scanpar,name,country) %>%
  summarise(sum_attackrate_per1000=sum(per_1000_popul)) %>% # ,sum_attackrate_agegr_frac=round(sum(n_case)/sum(agegr_pop),3)
  pivot_longer(cols=contains("attackrate"),names_to="vartype") # sum_attackrate_casenum=round(sum(n_case))
mean_age=unlist(lapply(cntr_list, function(x) {fun_meanage_cntr(popM,popF,x,year_str="2020")}))
attackrate_all_ages[,"aver_age"]=mean_age[match(attackrate_all_ages$country,cntr_list)]; 
attackrate_all_ages[,"cntrcode"]=countrycode(attackrate_all_ages$country,origin='country.name',destination='iso3c')
# plot summed attackrate as a fcn of param values
maxval=ceiling(max(attackrate_all_ages$value[grepl("clin",attackrate_all_ages$name)])/50)*50
sel_var=c("sum_attackrate_per1000","sum_attackrate_casenum","sum_attackrate_agegrfrac")[1] # sprintf("factor(%s)","par_scale")
cntr_table=attackrate_all_ages %>% group_by(cntrcode,scanpar,name) %>% filter(par_scale==max(par_scale)) %>% summarise(maxval_bypar=value)
ggplot(attackrate_all_ages[grepl("clin",attackrate_all_ages$name),],aes(x=par_scale,y=value,group=cntrcode,color=aver_age)) + 
  geom_line(size=1.3) + facet_grid(name~scanpar,switch="y",scales="free_y") + 
  scale_x_continuous(breaks=unique(attackrate_all_ages$par_scale)) + scale_color_gradientn(colours=pal) +
  geom_text(data=cntr_table[grepl("clin",cntr_table$name),],aes(x=6.2,y=maxval_bypar,label=cntrcode),size=3,color="black") +
  theme_bw() + standard_theme + xlab('age dependence (asymp. inf.: level)') + ylab("") + ggtitle('attack rates per 1000 popul') + 
  labs(color="<age>") # + guides(color=guide_legend())
# SAVE
cntrlistcd=paste0(countrycode(cntr_list,origin='country.name',destination='iso3c'),collapse="_")
ggsave(paste0("simul_output/clin_attackrate_sum_",cntrlistcd,"_byage.png"),width=25,height=10,units="cm")
# plot as fcn of mean age
cntr_table=attackrate_all_ages %>% group_by(country,scanpar,aver_age) %>% summarise(minval=min(value)) %>% 
              mutate(cntr_iso=countrycode(country,origin='country.name',destination='iso3c'))
ggplot(attackrate_all_ages[grepl("clin",attackrate_all_ages$name),],aes(x=aver_age,y=value,color=par_scale)) + # ,group=par_scale
  geom_line(aes(group=par_scale),size=1.3) + geom_point() + facet_grid(name~scanpar,switch="y",scales="free_y") + 
  scale_color_gradientn(colours=pal) + scale_x_continuous(breaks=seq(20,45,2.5)) +
  geom_text(data=cntr_table, aes(x=aver_age,y=minval-15,label=cntr_iso),size=3,color="black",angle=90) +
  theme_bw() + standard_theme + xlab('average age') + ylab("") + ggtitle('attack rates per 1000 popul') + labs(color="age dep")
ggsave(paste0("simul_output/clin_attackrate_sum_",cntrlistcd,"_ageplot.png"),width=25,height=10,units="cm")

### R0 and HIT ----------------------------
# for one cntr
countryval="Sudan"; N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups);covid_params=cm_parameters_SEI3R(gsub("Sudan","Ethiopia",countryval))
C_m_orig=Reduce('+',covid_params$pop[[1]]$matrices); C_m=fun_contmatr_recipr(C_m_orig,N_tot)
R0=fun_NGM_R0(C_m,N_tot,y_val,u_val,f_val,c(d_p,d_c,d_s))[[2]] # NGM=fun_NGM_R0(C_m,N_tot,y_val,u_val,f_val,c(d_p,d_c,d_s))[[1]]
# # covidm default: dIp=2.4; dIs=3.2; dIa=7 ; y_def=rep(0.5,length(N_tot)); u_def=rep(0.08,length(N_tot))
# PARSCAN for R0
min_val_susc=1e-2; maxval_susc=0.15; midpoint_susc=(min_val_susc+maxval_susc)/2; delta_susc=maxval_susc-midpoint_susc
midpoint_clinfract=0.4; delta_clinfr=0.35; rep_min_susc=3; rep_max_susc=5; rep_min_clinfr=5; rep_max_clinfr=3
# set params
agedep_param_lims=c(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,delta_clinfr,rep_min_clinfr,rep_max_clinfr)
# cntr_list=c("Sudan","Pakistan","South Africa","Brazil","Russian Federation","Italy")
# scan for n cntrs
k_max=6
for (k_cntr in 1:length(cntr_list)){countryval=cntr_list[k_cntr]
N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups); covid_params=cm_parameters_SEI3R(gsub("Sudan","Ethiopia",countryval))
C_m_orig=Reduce('+',covid_params$pop[[1]]$matrices); C_m=fun_contmatr_recipr(C_m_orig,N_tot)
R0_scan=fun_multidim_scan(countryval,C_m,N_tot,k_max,agedep_param_lims,c(d_p,d_c,d_s))
if (k_cntr==1){R0_scan_cntrs=R0_scan} else {R0_scan_cntrs=rbind(R0_scan_cntrs,R0_scan)} } 
# add cntr name and age
R0_scan_cntrs$country=factor(R0_scan_cntrs$country,levels=cntr_list);cntr_codes=countrycode(cntr_list,origin='country.name',destination='iso3c')
mean_age=unlist(lapply(cntr_list, function(x) {fun_meanage_cntr(popM,popF,x,year_str="2020")}))
R0_scan_cntrs[,"aver_age"]=mean_age[match(R0_scan_cntrs$country,cntr_list)]
# HIT line plot
plotval="R0" # HIT R0
# facet by cntr (colors=clin_fract)
ggplot(R0_scan_cntrs[R0_scan_cntrs$asympt_inf/k_max<1&R0_scan_cntrs$asympt_inf>0,],
 aes_string(x="susc",y=plotval,color="clinfract",group="clinfract")) + geom_line(size=1.2) + scale_color_gradientn(colours=pal) + 
 facet_grid(country~asympt_inf_str) + theme_bw() + standard_theme + theme(axis.text.y=element_text(size=8)) +
 xlab("age dependence of susceptibility") + ylab("HIT (=1-1/R0)") + guides(color=guide_legend(title="age dependence of\nclinical fraction")) +
 scale_x_continuous(breaks=unique(R0_scan_cntrs$susc)) + ggtitle(plotval) # 
# scale_color_brewer(palette="YlOrRd")
# SAVE
ggsave(paste0("simul_output/HITs/",plotval,"_",paste(cntr_codes,collapse="_"),".png"),width=30,height=25,units="cm")
# facet by clinfract and asympt infect (cntrs by different colors)
ggplot(R0_scan_cntrs[R0_scan_cntrs$asympt_inf>0,],aes(x=susc,y=HIT,color=country,group=interaction(clinfract_str,country))) +
  geom_line(size=1.2) + scale_color_manual(values=c("red","orange","pink","green","blue")) + 
  facet_grid(clinfract_str~asympt_inf_str,scales="free") + theme_bw() + standard_theme + # theme(axis.text.y=element_text(size=8)) +
  xlab("age dependence of susceptibility") + ylab("HIT (=1-1/R0)") +
  # guides(color=guide_legend(title="age dependence of\nclinical fraction")) +
  scale_x_continuous(breaks=unique(R0_scan_cntrs$susc)) + scale_y_continuous(breaks=(0:10)/10) + ggtitle("Herd immunity thresholds")
# save
ggsave(paste0("simul_output/HITs/HIT_",paste(cntr_codes,collapse="_"),"_facet_clinfract_asymptinfect.png"),
       width=30,height=24,units="cm")
############################
# plot by age
truthvals_R0=R0_scan_cntrs$asympt_inf>0&R0_scan_cntrs$asympt_inf<max(R0_scan_cntrs$asympt_inf)
cntr_table=R0_scan_cntrs[truthvals_R0,] %>% group_by(country,aver_age,asympt_inf,susc,asympt_inf_str,susc_str) %>%
  filter(clinfract==max(clinfract)) %>% summarise(hit_minval=HIT,r0_minval=R0) %>% 
  mutate(cntr_iso=countrycode(country,origin='country.name',destination='iso3c'))
plotval="R0" # HIT R0
ggplot(R0_scan_cntrs[truthvals_R0,],aes_string(x="aver_age",y=plotval,color="clinfract")) + geom_line(aes(group=clinfract),size=1.3) + 
  geom_point() + facet_grid(susc_str~asympt_inf_str,switch="y",scales="free_y") + scale_color_gradientn(colours=pal) +
  geom_text(data=cntr_table, aes(x=aver_age,y=r0_minval,label=cntr_iso),size=3,color="black",angle=90) +
  theme_bw() + standard_theme + xlab('average age') + ylab("") + ggtitle(plotval) # + labs(color="age dep")
# save
ggsave(paste0("simul_output/HITs/",plotval,"_",paste(cntr_codes,collapse="_"),"_byage.png"),
       width=30,height=24,units="cm")
### plot of parameter ranges
agedep_partable=fun_agedep_partable(agedep_param_lims,k_max); 
ggplot(agedep_partable,aes(x=agegroup,y=value,group=dep_fact,color=dep_fact)) + geom_line(size=1.25) + 
 scale_x_continuous(breaks=1:16) + facet_wrap(~name,scales="free") + scale_color_gradientn(colours=pal) +theme_bw()+ standard_theme
ggsave("simul_output/agedep_params_fval_uval_yval.png",width=30,height=12,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### Sudan linelist data ----------------------------
linelist_clean=readRDS("SDN_linelist_script/linelist_clean_2020-12-09.rds")
