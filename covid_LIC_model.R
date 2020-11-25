# simulation study to compare potential mechanisms explaining the lower than expected case and death numbers in LIC

# to render as html: rmarkdown::render("RSV_kenya_SA_calculation.R",output_dir='output/cea_plots/')
# library(tidyverse); library(reshape2); library(matrixStats); library(rstudioapi); # library(fitdistrplus)
library(tidyverse); library(qs); library(wpp2019); library(countrycode); library(deSolve)
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# plotting theme
source("covid_LIC_fcns.R")
### ### ### ### ### ### ### ### ### ### ### ###
### SEIR age structured (covidm) --------------------------

# SEIR age structured ODEs
seir_agestr_ode <- function(t,X,parms){ # params=list(K_m,C_m,b_m_full,u_val,a_m,ind_all_susceptibles,ind_all_inf_vars)
  K_m=parms[[1]];C_m=parms[[2]];b_m_full=parms[[3]]; u_val=parms[[4]];a_m=parms[[5]]; 
  susc_ind=parms[[6]]; inf_vars_inds=parms[[7]]
  dXdt=b_m_full %*% diag(X[susc_ind]) %*% C_m %*% diag(u_val) %*% a_m %*% X[inf_vars_inds] + K_m%*%X; list(dXdt) }

n_age=5; vartype_list=c('S','E','Ip','Ic','Is'); dim_sys=n_age*length(vartype_list)
infect_vartype=c('E','Ip','Ic','Is')
# full list of vars
full_varname_list=fun_seir_agestr_varnames(vartype_list,n_age)
# convert name & age index into linear index: fun_sub2ind_seir_agestr(j_age=2,varname='Is',vartype_list,length(vartype_list),n_age=2)
# total popul
N_tot=rep(1e6,n_age)
# params that could be age-dependent
u_val=rep(0.25,n_age)/N_tot; y_val=0.2+(1:n_age)*0.1; f_val=rep(0.5,n_age)
# constant transmission parameters
d_e=1; d_p=1; d_c=1; d_s=1; infect_first_ord_pars=c(d_e,d_p,d_c,d_s)
# set up kinetic matrix
K_m=fun_seir_agestr_kinmatr(infect_vartype,vartype_list,n_age,infect_first_ord_pars,y_val)
# infection terms
forceinf_vars=c("S","E"); n_lin_terms=sum(!vartype_list %in% forceinf_vars)
# vectors to scale force of infection terms
list_bm_am=fun_force_inf_vects(vartype_list,forceinf_vars,n_age,f_val); b_m_full=list_bm_am[[1]]; a_m=list_bm_am[[2]]
# contact matrix
C_m=matrix(1,n_age,n_age)
# inf vector is: b_m_full %*% diag(c(1,2)) %*% C_m %*% diag(u_val) %*% a_m %*% matrix(1,8,1)
ind_all_inf_vars=fun_inds_vartypes(n_age,vartype_list,infect_vartype)[[1]]
ind_all_susceptibles=fun_inds_vartypes(n_age,vartype_list,infect_vartype)[[2]]

# INITIAL CONDITIONS (full_varname_list=fun_seir_agestr_varnames(vartype_list,n_age))
# seed epidemic by "E">0 in a given (or multiple) age groups
ind_inf_initvals=fun_sub2ind_seir_agestr(j_age=1:2,varname="E",varname_list=vartype_list,n_var=length(vartype_list),n_age=n_age)
initvals_seir_model=matrix(0,dim_sys,1); initvals_seir_model[ind_all_susceptibles]=N_tot
inf_initval=10; initvals_seir_model[ind_inf_initvals]=inf_initval; initvals_m=matrix(initvals_seir_model,ncol=n_age)
initvals_seir_model[ind_all_susceptibles]=initvals_seir_model[ind_all_susceptibles]-colSums(initvals_m[2:nrow(initvals_m),])

### run ODEs ----------------------------
# time resolution and duration
n_days_year=365; n_years=1; max_time=n_years*n_days_year; timesteps <- seq(0,max_time,by=0.25)
# define inputs for ODE fcn
params=list(K_m, C_m, b_m_full, u_val, a_m, ind_all_susceptibles, ind_all_inf_vars)
# RUN # 
ptm<-proc.time(); ode_solution<-lsoda(initvals_seir_model,timesteps,func=seir_agestr_ode,parms=params); proc.time() - ptm
# PROCESS OUTPUT
# df_ode_solution = ode_solution %>% as.data.frame() %>% setNames(c("t",full_varname_list))
df_ode_solution=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars)[[1]]
df_ode_solution_tidy=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars)[[2]]

# PLOT
xval_breaks=seq(timesteps[1],timesteps[length(timesteps)],10); xlim_val=160
ggplot(df_ode_solution_tidy, aes(x=t,y=value,group=variable,color=compartm)) + # ,linetype=vartype
  geom_line(size=1.05) + # scale_linetype_manual("vars",values=c("dotdash","solid")) +
  facet_wrap(~agegroup+vartype,scales='free_y',nrow=n_age) + theme_bw() + standard_theme + labs(linetype='vars',color='vars') + #
  scale_x_continuous(breaks=xval_breaks,limits=c(0,xlim_val)) + # ,minor_breaks=; seq(timesteps[1],timesteps[length(timesteps)],10)
  xlab('days') + ylab('') + ggtitle(paste0(paste0(vartype_list,collapse="-"),'-R simulation'))

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
ind_inf_initvals=fun_sub2ind_seir_agestr(j_age=1,varname="E",varname_list=vartype_list,n_var=length(vartype_list),n_age=n_age)
inf_initval=10; S0_val=N_tot[1]-inf_initval; initvals_seir_model=matrix(0,dim_sys,1); initvals_seir_model[ind_inf_initvals]=inf_initval
initvals_seir_model[ind_all_susceptibles[1]]=S0_val
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

# from sympy (in python)
# dominant eigenvalue (2-age group SEIR model):
# (beta11*gamma2 + beta22*gamma1)/(2*gamma1*gamma2) + sqrt(beta11**2*gamma2**2 - 2*beta11*beta22*gamma1*gamma2 + 
# 4*beta12*beta21*gamma1*gamma2 + beta22**2*gamma1**2)/(2*gamma1*gamma2)

#### loglikelihood with negbinom and dirichlet distribs ---------------
# loglikelihood ~ prod_k(negbinom(incid_obs_day_k,pred_obs_day_k))*prod_m(dirmultinomial(obs_age_distrib,pred_age_distrib))
# num_k=1e3; sse_lglkl=matrix(NA,nrow=num_k,ncol=3)
# for (k in 1:num_k) { n_days=44; incid_pred=round(exp((0:n_days)/5.25)); 
# errorscale=round(incid_pred*0.001); errorscale[errorscale<1]=1; errorsvals=errorscale*floor(runif(n_days+1,min=-10,max=10))
# incid_obs=incid_pred+errorsvals; incid_obs[incid_obs<0]=0; mean_sse=mean((incid_obs-incid_pred)^2)
# mean_perc_err=mean(abs((incid_pred- incid_obs)/incid_pred))
# true_age_distrib=matrix(rep(c(1,20,50,30)/1e2,n_days),nrow=n_days,byrow=T)
# pred_age_distrib=matrix(rep(c(1,20,45,35)/1e2,n_days),nrow=n_days,byrow=T)
# # install.packages("extraDistr") # library("extraDistr")
# binm_probs=dnbinom(x=incid_obs,size=200,mu=incid_pred)
# dirichlet_probs=sapply(1:nrow(true_age_distrib),function(x) {ddirichlet(x=true_age_distrib[x,],alpha=pred_age_distrib[x,])})
# sse_lglkl[k,]=c(mean_perc_err,mean_sse,log10(prod(binm_probs)*prod(dirichlet_probs)))
# }

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

lmic_symptom_fract=c(rep(0.2973718, 2), rep(0.2230287, 2), rep(0.4191036, 2), rep(0.4445867, 2), rep(0.5635720, 2), rep(0.8169443, 6) )
# library(data.table)
age_groups <- data.table( age_group = c(1:16),   age_low = c( seq(0, 75, 5) ),   age_high = c( seq(4, 74, 5), 100) )
covid_parameters <- list( "clinical_fraction" = list( values = list( "age_y" = 19, "age_m" = 50, "age_o" = 68,
      "symp_y" = 0.037, "symp_m" = 0.3, "symp_o" = 0.65 )  ) )
#' calculate probability of clinical disease by age group
#' following posterior estimates by Davies et al
getClinicalFraction <- function(age_groups){
  #' smoothly interpolate between points (x0, y0) and (x1, y1) using cosine interpolation.
  #' for x < x0, returns y0; for x > x1, returns y1; for x0 < x < x1, returns the cosine interpolation between y0 and y1
  interpolate_cos = function(x, x0, y0, x1, y1)
  {
    ifelse(x < x0, y0, ifelse(x > x1, y1, y0 + (y1 - y0) * (0.5 - 0.5 * cos(pi * (x - x0) / (x1 - x0)))))
  }
  age_groups[, mid := mean(c(age_low, age_high)), by=seq_len(nrow(age_groups))]
  age_y = covid_parameters[["clinical_fraction"]][["values"]][["age_y"]]
  age_m = covid_parameters[["clinical_fraction"]][["values"]][["age_m"]]
  age_o = covid_parameters[["clinical_fraction"]][["values"]][["age_o"]]
  # definition of "young", "middle", and "old"
  young  = interpolate_cos(age_groups[, mid], age_y, 1, age_m, 0);
  old    = interpolate_cos(age_groups[, mid], age_m, 0, age_o, 1);
  middle = 1 - young - old;
  symp_y = covid_parameters[["clinical_fraction"]][["values"]][["symp_y"]]
  symp_m = covid_parameters[["clinical_fraction"]][["values"]][["symp_m"]]
  symp_o = covid_parameters[["clinical_fraction"]][["values"]][["symp_o"]]
  return(young * symp_y + middle * symp_m + old * symp_o)
}


### YYG model param fits ------------------------------------
# read in json file of param fits by youyanggu model
# install.packages("rjson"); library(rjson); library(jsonlite)
# hungary_params <- jsonlite::fromJSON( "../yyg-seir-simulator/best_params/latest/global/Hungary_ALL.json" )
# # fromJSON(file="yyg-seir-simulator/best_params/latest/global/Hungary_ALL.json")  # %>% as.data.frame
# hu_paramfits=data.frame(rbind(cbind(hungary_params$top10_params,rep('best',nrow(hungary_params$top10_params))), 
#       cbind(hungary_params$mean_params,rep('mean',nrow(hungary_params$top10_params))), 
#       cbind(hungary_params$median_params,rep('median',nrow(hungary_params$top10_params)))),country="HU")
# france_params <- jsonlite::fromJSON( "../yyg-seir-simulator/best_params/latest/global/France_ALL.json" )
# # fromJSON(file="yyg-seir-simulator/best_params/latest/global/Hungary_ALL.json")  # %>% as.data.frame
# fr_paramfits=data.frame(rbind(cbind(france_params$top10_params,rep('best',nrow(france_params$top10_params))), 
#                               cbind(france_params$mean_params,rep('mean',nrow(france_params$top10_params))), 
#                               cbind(france_params$median_params,rep('median',nrow(france_params$top10_params)))),country="FR")
# View(merge(fr_paramfits,hu_paramfits,by = c("X1","X3")))
######################
# # MSE vs MLE
# x<-rnorm(100,mean=20); b0<-10; b1<-20; y <- b1*x + b0 + rnorm(100);  df <- data.frame(x = x , y = y)
# #
# model <- lm(data = df , formula = y ~ x)
# # loglikelihood
# loglikelihood <- function(b0,b1){-sum(dnorm(df$y-(df$x*b1)-b0 , log=TRUE)) }    #-sum(log(R))
# # library(stats4)
# mle(loglikelihood, start = list(b0=1,b1=1))