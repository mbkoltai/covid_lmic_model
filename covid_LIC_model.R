# simulation study to compare potential mechanisms explaining the lower than expected case and death numbers in LIC

# to render as html: rmarkdown::render("RSV_kenya_SA_calculation.R",output_dir='output/cea_plots/')
# library(tidyverse); library(reshape2); library(matrixStats); library(rstudioapi); # library(fitdistrplus)
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
library(tidyverse); library(qs); library(wpp2019); library(countrycode); library(deSolve)
# functions and plotting theme
source("covid_LIC_fcns.R")
### ### ### ### ### ### ### ### ### ### ### ###
### SEIR age structured (covidm) --------------------------

# SEIR age structured ODEs
seir_agestr_ode <- function(t,X,parms){ # params=list(K_m,C_m,b_m_full,u_val,a_m,ind_all_susceptibles,ind_all_inf_vars)
  K_m=parms[[1]];C_m=parms[[2]];b_m_full=parms[[3]]; u_val=parms[[4]];a_m=parms[[5]]; 
  susc_ind=parms[[6]]; inf_vars_inds=parms[[7]]
  dXdt=b_m_full %*% diag(X[susc_ind]) %*% C_m %*% diag(u_val) %*% a_m %*% X[inf_vars_inds] + K_m%*%X; list(dXdt) }
# SYSTEM DIMENSIONS
n_age=16; vartype_list=c('S','E','Ip','Ic','Is'); infect_vartype=c('E','Ip','Ic','Is'); dim_sys=n_age*length(vartype_list)
# full list of vars
full_varname_list=fun_seir_agestr_varnames(vartype_list,n_age)
# convert name+age index to linear index: fun_sub2ind_seir_agestr(j_age=2,varname='Is',vartype_list,length(vartype_list),n_age=2)
# TOTAL POPULATION
# age struct in covidm
age_groups <- data.table(age_group=c(1:16), age_low=c(seq(0,75,5) ), age_high=c(seq(4,74,5),100))
# population data from wpp2019
i_cntr="Sudan";cntr_agestr=data.frame(agegroups=popF[popF$name %in% i_cntr,"age"],values=popF[popF$name %in% i_cntr,"2020"] +
                                        popM[popM$name %in% i_cntr,"2020"])
agegr_truthvals=sapply(strsplit(as.character(cntr_agestr$agegroups),"-"),"[[",1) %in% age_groups$age_low
N_tot=cntr_agestr$values[agegr_truthvals];N_tot[length(N_tot)]=N_tot[length(N_tot)]+sum(cntr_agestr$values[!agegr_truthvals])
# wpp2019 data is in thousands
if(all(N_tot<1e6)) {N_tot=N_tot*1e3}
# N_tot=rep(1e6,n_age)
# age-dependent parameters
u_val=rep(0.3,n_age)/N_tot; y_val=0.04+(1:n_age)*0.6/n_age; f_val=rep(0.5,n_age)
# constant transmission parameters
d_e=4; d_p=1.5; d_c=3.5; d_s=5; infect_first_ord_pars=c(d_e,d_p,d_c,d_s)
# KINETIC MATRIX (LINEAR TERMS)
K_m=fun_seir_agestr_kinmatr(infect_vartype,vartype_list,n_age,infect_first_ord_pars,y_val)
# vectors to scale FORCE of INFECTION terms
list_bm_am=fun_force_inf_vects(vartype_list,forceinf_vars=c("S","E"),n_age,f_val); b_m_full=list_bm_am[[1]]; a_m=list_bm_am[[2]]
# CONTACT MATRIX
# call covidm for Ethiopia contact matrix
if (!exists("covid_params")){cm_path="~/Desktop/research/models/epid_models/covid_model/lmic_model/covidm/"
source(file.path(cm_path,"R","covidm.R"));cm_force_rebuild=F;cm_build_verbose=T;cm_version=2
covid_params=cm_parameters_SEI3R("Ethiopia")}
C_m_covidm=covid_params$pop[[1]]$matrices$home + covid_params$pop[[1]]$matrices$work + covid_params$pop[[1]]$matrices$school
C_m=C_m_covidm # matrix(1,n_age,n_age)
# inf vector is: b_m_full %*% diag(c(1,2)) %*% C_m %*% diag(u_val) %*% a_m %*% matrix(1,8,1)
l=fun_inds_vartypes(n_age,vartype_list,infect_vartype); ind_all_inf_vars=l[[1]]; ind_all_susceptibles=l[[2]]
# INITIAL CONDITIONS
# seed epidemic by "E">0 in a given (or multiple) age groups
initvals_seir_model=fcn_set_init_conds(inf_initval=10,init_inf_age_groups=7,init_inf_vartype="E",
                                       n_age,N_tot,vartype_list,ind_all_susceptibles)

### run ODEs ----------------------------
# time resolution and duration
n_days_year=365; n_years=1; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=0.25)
# define inputs for ODE fcn
params=list(K_m, C_m, b_m_full, u_val, a_m, ind_all_susceptibles, ind_all_inf_vars)
# RUN # 
ptm<-proc.time(); ode_solution<-lsoda(initvals_seir_model,timesteps,func=seir_agestr_ode,parms=params); proc.time() - ptm
# PROCESS OUTPUT
# df_ode_solution = ode_solution %>% as.data.frame() %>% setNames(c("t",full_varname_list))
l_proc_sol=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars)
df_ode_solution=l_proc_sol[[1]]; df_ode_solution_tidy=l_proc_sol[[2]]

# PLOT
xval_breaks=seq(timesteps[1],timesteps[length(timesteps)],10); xlim_val=90; standard_theme$axis.text.y$size=6
df_age_groups=df_ode_solution_tidy %>% group_by(agegroup,vartype) %>% summarise(max_val=max(value),
        min_val=min(value),opt_val=(max(value)+min(value))/2,opt_fract_val=(max(fract_value)+min(fract_value))/2)
df_age_groups[,"agegroup_str"]=paste("age:",df_age_groups$agegroup)
# absolute or % value?
abs_or_fract=1
if (abs_or_fract==1){value_col="value";label_col="opt_val";savetag="absval"} else {
  value_col="fract_value"; label_col="opt_fract_val";savetag="fractval"}
ggplot(df_ode_solution_tidy,aes_string(x="t",y=value_col,group="variable",color="compartm")) + geom_line(size=1.05) + 
  facet_wrap(~agegroup+vartype,scales='free_y',ncol=4) + 
  geom_text(data=df_age_groups,aes_string(x=7,y=label_col,label="agegroup_str",group=NULL),size=3,color="black") +
  theme_bw() + standard_theme + theme(strip.background=element_blank(),strip.text=element_blank()) + 
  labs(linetype='vars',color='vars') + scale_x_continuous(breaks=xval_breaks,limits=c(0,xlim_val)) + 
  xlab('days') + ylab('') + ggtitle(paste0(paste0(vartype_list,collapse="-"),'-R simulation'))
# SAVE
ggsave(paste0("plots/SEIR_age_str_",n_age,"agegroups_",savetag,".png"),width=30,height=20,units="cm")
# ,minor_breaks=; seq(timesteps[1],timesteps[length(timesteps)],10)
# facet_grid(agegroup~vartype,rows=vars(agegroup),cols=vars(vartype),scales="free") + # 
# scale_linetype_manual("vars",values=c("dotdash","solid")) + # data=data.frame(agegroup=age_range),

# infected fraction of popul (HIT)
1 - round(df_ode_solution[nrow(df_ode_solution),grepl('S_',colnames(df_ode_solution))]/N_tot,2)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### algebraic solution

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### SEIR with 1 age group ----------------------------
# for one age group: set initconds for other age groups to 0
# R0=N_tot*u_val*(y_val*((1/d_c + 1/d_p) - f_val/d_s) + f_val/d_s)
# N_tot*u_val*(y_val*(1/d_c - f_val/d_S) + f_val/d_S)
# INITCOND
inf_initval=10; initvals_seir_model=fcn_set_init_conds(inf_initval,init_inf_age_groups=1,init_inf_vartype="E",n_age,N_tot,
                                       vartype_list,ind_all_susceptibles)
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
ttl_str=expression("Herd immunity threshold ~ susceptibility ("~sigma~"), asympt. infectiousness (f), clinical fraction ("~gamma~")")
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
sample_size=1e4
# d_e (pre-infectious days)
mean_g=4; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# d_p (presympt inf.)
mean_g=1.5; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# dC (Duration of symptomatic infectiousness in days)
mean_g=3.5; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# dS (Duration of asymptomatic infectiousness in days)
mean_g=5; shape_g=4; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# R0
mean(rnorm(sample_size,2.6,0.5))
# Delay from symptom onset to becoming a severe case in days (mean: 7 days)
mean_g=7; shape_g=8; mean(rgamma(sample_size,shape=shape_g,scale=mean_g/shape_g))
# Duration of severe, non-critical disease in days ~ Gamma(μ=8,k=8)
# Proportion of severe cases that become critical: 30%
# Duration of severe, critical disease in days ~ Gamma(μ = 10, k = 10)
# Delay from symptom onset to death in days ~ Gamma(μ=22,k=22)

lmic_symptom_fract=c(rep(0.2973718,2), rep(0.2230287,2), rep(0.4191036,2), rep(0.4445867,2), rep(0.5635720,2), rep(0.8169443,6))
# library(data.table)
age_groups <- data.table(age_group=c(1:16), age_low=c(seq(0,75,5) ), age_high=c( seq(4, 74, 5), 100) )
covid_parameters<-list("clinical_fraction"=list(values=list("age_y"=19,"age_m"=50,"age_o"=68,
                                                            "symp_y"=0.037,"symp_m"=0.3,"symp_o"=0.65)) )

# susceptibility from "Age-dependent effects...", 10-age bands
suscept_vals=c(0.4,0.38,0.79,0.86,0.8,0.82,0.88,0.74)
suscept_mean_age=4.5+(0:6)*10