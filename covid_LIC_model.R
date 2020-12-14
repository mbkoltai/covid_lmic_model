# simulation study to compare potential mechanisms explaining the lower than expected case and death numbers in LIC

# to render as html: rmarkdown::render("RSV_kenya_SA_calculation.R",output_dir='output/cea_plots/')
# library(tidyverse); library(reshape2); library(matrixStats); library(rstudioapi); # library(fitdistrplus)
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
library(tidyverse); library(wpp2019); library(countrycode); library(deSolve); # library(qs); library(data.table)
covidm_abs_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
# functions and plotting theme
source("covid_LIC_fcns.R")
### ### ### ### ### ### ### ### ### ### ### ###
### SEIR age structured (covidm) --------------------------

# SEIR age structured ODEs
seir_agestr_ode <- function(t,X,parms){ # params=list(K_m,C_m,b_m_full,u_val,a_m,ind_all_suscept,ind_all_inf_vars)
  K_m=parms[[1]];C_m=parms[[2]];b_m_full=parms[[3]]; u_val=parms[[4]];a_m=parms[[5]]; 
  susc_ind=parms[[6]]; inf_vars_inds=parms[[7]]
  dXdt=b_m_full %*% diag(X[susc_ind]) %*% C_m %*% diag(u_val) %*% a_m %*% X[inf_vars_inds] + K_m%*%X; list(dXdt) }
# SYSTEM DIMENSIONS
n_age=16; vartype_list=c('S','E','Ip','Ic','Is'); infect_vartype=c('E','Ip','Ic','Is'); dim_sys=n_age*length(vartype_list)
# full list of vars
full_varname_list=fun_seir_agestr_varnames(vartype_list,n_age)
# name+age --> linear index: fun_sub2ind_seir_agestr(j_age=2,varname='Is',vartype_list,length(vartype_list),n_age=2)
# POPULATION (age struct resol as in covidm)
age_groups <- data.frame(age_group=c(1:16), age_low=c(seq(0,75,5) ), age_high=c(seq(4,74,5),100))
# population data from wpp2019
countryval="Sudan"; N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups) # uniform: N_tot=rep(1e6,n_age)
# CONTACT MATRIX
if (!exists("covid_params")){setwd("covidm/"); cm_force_rebuild=F; cm_build_verbose=T; cm_version=2; 
source(file.path(covidm_abs_path,"R","covidm.R")); setwd("../")}
covid_params=cm_parameters_SEI3R(gsub("Sudan","Ethiopia",countryval))
C_m=Reduce('+',covid_params$pop[[1]]$matrices)
# constant transmission parameters
d_e=1/3; d_p=1/2; d_c=1/3; d_s=1/5; infect_first_ord_pars=c(d_e,d_p,d_c,d_s)
########################################################
### AGE DEPENDENT parameters ---------------------------
min_val_susc=1e-2; maxval_susc=0.1; midpoint_susc=(min_val_susc+maxval_susc)/2;delta_susc=maxval_susc-midpoint_susc; mean_susc_exp=0.05
mean_clinfract=0.307; midpoint_clinfract=0.4; delta_clinfr=0.305; min_val_clinfr=0.05
### SINGLE SIMULATION ---------------------------
# set params
depval=1; u_val=fun_lin_approx_agedep_par(min_val_susc,max_val=maxval_susc,rep_min=3,rep_max=5) # /sum(N_tot)
y_val=fun_lin_approx_agedep_par(min_val_clinfr,max_val=midpoint_clinfract+depval*delta_clinfr,rep_min=5,rep_max=3); f_val=rep(0.5,n_age)
# KINETIC MATRIX (LINEAR TERMS)
K_m=fun_seir_agestr_kinmatr(infect_vartype,vartype_list,n_age,infect_first_ord_pars,y_val)
# vectors to scale FORCE of INFECTION terms
list_bm_am=fun_force_inf_vects(vartype_list,forceinf_vars=c("S","E"),n_age,f_val,N_tot); b_m_full=list_bm_am[[1]];a_m=list_bm_am[[2]]
# inf vector is: b_m_full %*% diag(c(1,2)) %*% C_m %*% diag(u_val) %*% a_m %*% matrix(1,8,1)
l=fun_inds_vartypes(n_age,vartype_list,infect_vartype); ind_all_inf_vars=l[[1]]; ind_all_suscept=l[[2]]
# INITIAL CONDITIONS (seed epidemic by "E">0 in a given (or multiple) age groups)
initvals_seir_model=fcn_set_init_conds(inf_initval=10,init_inf_age_groups=7,init_inf_var="E",n_age,N_tot,vartype_list,ind_all_suscept)
### run ODEs ----------------------------
# time resolution and duration
n_days_year=365; n_years=1; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=0.25)
# define inputs for ODE fcn
params=list(K_m, C_m, b_m_full, u_val, a_m, ind_all_suscept, ind_all_inf_vars)
# RUN simulation
ptm<-proc.time(); ode_solution<-lsoda(initvals_seir_model,timesteps,func=seir_agestr_ode,parms=params); proc.time() - ptm
# PROCESS OUTPUT
# df_ode_solution = ode_solution %>% as.data.frame() %>% setNames(c("t",full_varname_list))
l_proc_sol=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars); df_ode_solution=l_proc_sol[[1]]
# final values
final_vals=fun_final_vals(df_ode_solution,N_tot,"fract"); sum((1-final_vals)*N_tot)/sum(N_tot)
# PLOT DYNAMICS (df_ode_solution_tidy=l_proc_sol[[2]])
fun_seir_agegroups_dyn(age_groups,l_proc_sol[[2]],xlim_val=220,abs_or_fract=2,timesteps,standard_theme)
ggsave(paste0("covid_simul_output/",countryval,"_simul/SEIR_age_str_",n_age,"agegroups_",savetag,".png"),width=30,height=20,units="cm")
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### PARAMETER SCAN ---------------------------
cntr_list=c("Sudan","South Africa","Brazil","Italy");scan_param_fullname=c("susceptibility","clinical fraction","asymptomatic infectiousness")
# average age: fun_meanage_cntr(popM,popF,"Pakistan",year_str="2020")
attack_rates_multipar_multicntr=fun_multipar_multicntr_scan(cntr_list,year_str="2020",age_groups,covidm_abs_path,scan_param_fullname,
          n_loop=6,midpoint_susc,delta_susc,midpoint_clinfract,delta_clinfr,mean_susc_exp,mean_clinfract,n_age,
          infect_first_ord_pars,infect_vartype,vartype_list,inf_initval=10,init_inf_age_groups=7,init_inf_var="E",n_years)
# PLOT attack rates as fcn of par values FOR ONE CNTR
sel_var=c("fract_agegroup","fraction_cases","per_1000_popul")[2]; truthval_cntr=grepl("clin",attack_rates_multipar_multicntr$name)
ggplot(attack_rates_multipar_multicntr[truthval_cntr,], aes_string(x="agegroup_names",y=sel_var,group="col_scale",color="col_scale")) + 
 geom_line() + geom_point() + facet_grid(country~scanpar,switch="y") + # ,scales='free_y'  facet_wrap(~name+scanpar,scales="free_y") +
 theme_bw() + standard_theme + xlab("age group") + ggtitle(paste0("clinical attack rates (",gsub("value","fraction",sel_var),")")) +
 labs(color="age-dependence/level \nof transmission \nparameter") # + scale_y_continuous(breaks=seq(0,1,0.05))
# SAVE
ggsave(paste0("covid_simul_output/attackrate_",sel_var,"_agedep.png"),width=30,height=15,units="cm")
#######
# summed attackrates (across age groups)
attackrate_all_ages=attack_rates_multipar_multicntr %>% group_by(par_scale,scanpar,name,country) %>%
  summarise(sum_attackrate_per1000=sum(per_1000_popul)) %>% # ,sum_attackrate_agegr_frac=round(sum(n_case)/sum(agegr_pop),3)
  pivot_longer(cols=contains("attackrate"),names_to="vartype") # sum_attackrate_casenum=round(sum(n_case))
# total attackrate
maxval=ceiling(max(attackrate_all_ages$value[grepl("clin",attackrate_all_ages$name)])/50)*50
sel_var=c("sum_attackrate_per1000","sum_attackrate_casenum","sum_attackrate_agegrfrac")[1] # sprintf("factor(%s)","par_scale")
# [grepl("clin",attackrate_all_ages$name),]
ggplot(attackrate_all_ages,aes(x=par_scale,y=value,color=country)) + geom_line(size=1.3) + 
  facet_grid(name~scanpar,switch="y",scales="free_y") + theme_bw() + standard_theme + xlab('age dependence (asymp. inf.: level)') +
  ylab("") + ggtitle('attack rates per 1000 popul') + labs(color="countries") + 
  scale_x_continuous(breaks=unique(attackrate_all_ages$par_scale)) # + scale_y_continuous(limits=c(50,maxval),breaks=seq(0,maxval,50))
ggsave(paste0("covid_simul_output/attackrates_sum_all_cntrs.png"),width=25,height=10,units="cm")

### linelist data ----------------------------
linelist_clean=readRDS("SDN_linelist_script/linelist_clean_2020-12-09.rds")

### R0 and HIT ----------------------------
countryval="Sudan"; N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups)
covid_params=cm_parameters_SEI3R(gsub("Sudan","Ethiopia",countryval)); C_m=Reduce('+',covid_params$pop[[1]]$matrices)
R0=fun_NGM_R0(N_tot,y_val,u_val,f_val,C_m,c(d_p,d_c,d_s))[[2]] # NGM=fun_NGM_R0(N_tot,y_val,u_val,f_val,C_m)[[1]]
# # covidm default: dIp=2.4; dIs=3.2; dIa=7 ; y_def=rep(0.5,length(N_tot)); u_def=rep(0.08,length(N_tot))
# NGM_def[i_row,j_col]=u_def[i_row]*C_m[i_row,j_col]*(y_def[j_col]*(2.4+3.2) + (1-y_def[j_col])*unique(f_val)*7) 
# PARSCAN for R0
min_val_susc=1e-2; maxval_susc=0.15; midpoint_susc=(min_val_susc+maxval_susc)/2; delta_susc=maxval_susc-midpoint_susc # mean_susc_exp=0.05
midpoint_clinfract=0.4; delta_clinfr=0.35 # mean_clinfract=0.307; min_val_clinfr=0.05
rep_min_susc=3; rep_max_susc=5; rep_min_clinfr=5; rep_max_clinfr=3
### SINGLE SIMULATION ---------------------------
# set params
agedep_param_lims=c(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,delta_clinfr,rep_min_clinfr,rep_max_clinfr)
R0_scan=fun_multidim_scan(N_tot,k_max,agedep_param_lims,c(d_p,d_c,d_s))
# mycolors <- colorRampPalette(brewer.pal(8, "YlOrRd"))(length(unique(R0_scan$susc)))
# HIT line plot
ggplot(R0_scan,aes(x=susc,y=HIT,color=clinfract,group=clinfract)) + geom_line() + facet_wrap(~asympt_inf_str) + 
  scale_color_distiller(palette="Spectral") + # scale_color_manual(values = mycolors) + # scale_color_brewer(palette="YlOrRd") + # 
  theme_bw() + standard_theme + xlab("susceptibility") + ylab("clinical fraction") + ylim(c(0.33,0.92)) + 
  ggtitle(paste0(countryval," herd immunity threshold"))
ggsave(paste0("covid_simul_output/HIT_",countryval,".png"),width=20,height=12,units="cm")

# HIT heatmap
# ggplot(R0_scan,aes(x=susc,y=clinfract,fill=HIT)) + geom_tile() + facet_wrap(~asympt_inf_str) + scale_fill_distiller(palette="Spectral") +
#   geom_text(aes(label = round(HIT,2)),size=3) + theme_bw() + standard_theme + xlab("susceptibility") + ylab("clinical fraction") + 
#   ggtitle(paste0(countryval," herd immunity threshold"))


# delete repetitions in par table
# age_dep_paramtable_clean=fun_paramtable_norep(age_dep_paramtable,scan_parameter)
# # plot parameter values
# ggplot(age_dep_paramtable_clean,aes(x=age_group,y=value,group=interaction(name,par_scale_cnt),color=name)) + 
#   geom_line(aes(size=par_scale_cnt)) + scale_size(range=c(0.25,2)) + theme_bw() + standard_theme + 
#   labs(color="transmission parameter",size='age dependence') + scale_y_continuous(breaks=(0:20)/20) +
#   ggtitle(paste(scan_param_fullname[grepl(gsub("_fract","",scan_parameter),scan_param_fullname)],"age dependence"))
# # SAVE
# ggsave(paste0("covid_simul_output/",countryval,"_simul/agedep_",scan_parameter,".png"),width=30,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### algebraic solution ----------------------------
### SEIR with 1 age group 
# for one age group: set initconds for other age groups to 0
# R0=N_tot*u_val*(y_val*((1/d_c + 1/d_p) - f_val/d_s) + f_val/d_s)
# N_tot*u_val*(y_val*(1/d_c - f_val/d_S) + f_val/d_S)
# INITCOND
inf_initval=10; initvals_seir_model=fcn_set_init_conds(inf_initval,init_inf_age_groups=1,init_inf_var="E",n_age,N_tot,
                                       vartype_list,ind_all_suscept)
# RUN
ptm<-proc.time(); ode_solution<-lsoda(initvals_seir_model,timesteps,func=seir_agestr_ode,parms=params); proc.time() - ptm
# PROCESS OUTPUT
df_ode_solution=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars)[[1]]
df_ode_solution_tidy=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars)[[2]]

### Analytical calc stationary solution (1 age group) ------------------
df_ode_solution[nrow(df_ode_solution),'S_1']/N_tot[1]
# part of the popul infected (HIT)
1 - df_ode_solution[nrow(df_ode_solution),'S_1']/N_tot[1]
# R0: R0
# analytical solution to S_eq
library("lamW") # install.packages("lamW"); 
Is0=initvals_seir_model[full_varname_list %in% "Is_1"]/N_tot[1]; E0=initvals_seir_model[full_varname_list %in% "E_1"]/N_tot[1]
Ic0=initvals_seir_model[full_varname_list %in% "Ic_1"]/N_tot[1]; Ip0=initvals_seir_model[full_varname_list %in% "Ip_1"]/N_tot[1]
S0=S0_val/N_tot[1]
sigmaval_renorm=u_val[1]*N_tot[1]
S_stat_sol=-lambertW0(-R0*S0*exp((-R0*S0*d_c*d_p*d_s + d_c*d_p*f_val[1]*sigmaval_renorm*(-E0*y_val[1] + E0 + Is0) + 
              d_c*d_s*sigmaval_renorm*(E0*y_val[1] + Ip0) + d_p*d_s*sigmaval_renorm*(E0*y_val[1] + Ic0 + Ip0))/(d_c*d_p*d_s)))/R0
S_stat_sol_simplif=-lambertW0(-R0*S0*exp(-R0*S0 ))/R0

# scan in the 3 params
gamma_range=seq(5,75,5)/100; f_range=seq(10,80,5)/100; sigma_range=seq(8,36,2)/10
param_scan_table=expand.grid(gamma_range,f_range,sigma_range)
S_stat_sol_vals=array()
# solve for S_ss
for (k in 1:nrow(param_scan_table)) {
  gamma=param_scan_table[k,1]; f_val=param_scan_table[k,2]; sigmaval_renorm=param_scan_table[k,3]
  R0=sigmaval_renorm*(gamma*((1/d_c + 1/d_p) - f_val/d_s) + f_val/d_s)
  S_stat_sol_vals[k]=-lambertW0(-R0*S0*exp((-R0*S0*d_c*d_p*d_s + d_c*d_p*f_val*sigmaval_renorm*(-E0*gamma + E0 + Is0) + 
      d_c*d_s*sigmaval_renorm*(E0*gamma + Ip0) + d_p*d_s*sigmaval_renorm*(E0*gamma + Ic0 + Ip0))/(d_c*d_p*d_s)))/R0 }
param_scan_sols_table=cbind(param_scan_table,S_stat_sol_vals); colnames(param_scan_sols_table)[1:3]=c('gamma','f_val','u_vals')
param_scan_sols_table[,"HIT"]=1-param_scan_sols_table$S_stat_sol_vals

# plot HIT
# facet by infectiousness
ttl_str=expression("attack rate ~ susceptibility ("~sigma~"), 
                   asympt. infectiousness (f), clinical fraction ("~gamma~")")
# facet by infectiousness
ggplot(param_scan_sols_table,aes(x=u_vals,y=HIT,group=factor(gamma),color=factor(gamma))) + geom_line() + 
  facet_wrap(~f_val) +  theme_bw() + standard_theme + 
  geom_text(data=data.frame(f_val=f_range),aes(x=0.98,y=0.98,label=paste0("f=",f_val)),size=4,color="black") +
  xlab(expression(sigma~"(susceptibility)")) + labs(color=expression(gamma)) + ggtitle(ttl_str)
# SAVE
ggsave(paste0("plots/",paste0(varname_list,collapse=""),"_analyt_sol_facetcontag.png"),width=30,height=20,units="cm")

# facet by clinical fraction
# ggplot(param_scan_sols_table,aes(x=u_vals,y=HIT,group=factor(f_val),color=factor(f_val))) + geom_line() + 
#   facet_wrap(~gamma) +  theme_bw() + standard_theme + 
#   geom_text(data=data.frame(gamma=gamma_range),aes(x=1.1,y=0.98,label=paste0("g=",gamma_range)),size=4,color="black") +
#   xlab(expression(sigma~"(susceptibility)")) + labs(color="f (infectiousness)") + ggtitle(ttl_str)
# # SAVE
# ggsave(paste0("plots/",paste0(varname_list,collapse=""),"_analyt_sol_facetclinfract.png"),width=30,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### R0 toy model
# r0vals=seq(6,30,2)/10; S_vals=(1:100)/100
# R0_scan=data.frame(melt(data.frame(S_ss=S_vals,sapply(r0vals,function(x) {x*(S_vals-1)-log(S_vals)})),id.vars="S_ss"),
#                    R0=unlist(lapply(r0vals,function(x) {rep(x,length(S_vals))})))
# ggplot(R0_scan,aes(x=S_ss,y=value,group=R0,color=R0)) + geom_line() + geom_hline(yintercept=0,color='red',linetype='dashed') + 
#   theme_bw() + standard_theme

### covidm parameters ------------------------------------------------------------------
# sample_size=1e4
# # d_e (pre-infectious days)
# mean_g=4; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# # d_p (presympt inf.)
# mean_g=1.5; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# # dC (Duration of symptomatic infectiousness in days)
# mean_g=3.5; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# # dS (Duration of asymptomatic infectiousness in days)
# mean_g=5; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# # R0
# mean(rnorm(sample_size,2.6,0.5))
# # Delay from symptom onset to becoming a severe case in days (mean: 7 days)
# mean_g=7; shape_g=8; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# # Duration of severe, non-critical disease in days ~ Gamma(μ=8,k=8)
# # Proportion of severe cases that become critical: 30%
# # Duration of severe, critical disease in days ~ Gamma(μ = 10, k = 10)
# # Delay from symptom onset to death in days ~ Gamma(μ=22,k=22)
# 
# lmic_symptom_fract=c(rep(0.2973718,2), rep(0.2230287,2), rep(0.4191036,2), rep(0.4445867,2), rep(0.5635720,2), rep(0.8169443,6))
# # library(data.table)
# age_groups <- data.table(age_group=c(1:16), age_low=c(seq(0,75,5) ), age_high=c( seq(4, 74, 5), 100) )
# covid_parameters<-list("clinical_fraction"=list(values=list("age_y"=19,"age_m"=50,"age_o"=68,
#                                                             "symp_y"=0.037,"symp_m"=0.3,"symp_o"=0.65)) )
# plot(age_groups$mid,getClinicalFraction(age_groups),type="b")
