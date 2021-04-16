# from name and age group to linear index ----------------------------------------------------------
fun_sub2ind_seir_agestr=function(j_age,varname,varname_list,n_var,n_age){
  varnum=which(varname_list %in% varname); k=(j_age-1)*n_var + varnum; k }

# generate all model names ----------------------------------------------------------
fun_seir_agestr_varnames=function(varname_list,n_age){
  array( sapply(1:n_age, function(x_age) {sapply(varname_list, function(x) {paste0(x,'_',x_age)}) } ) )
}

### function to have COVIDM incidence outputs
cm_multinom_process <- function(
  src, outcomes, delays,
  report = ""
) {
  if ("null" %in% names(outcomes)) {
    if (length(report) != length(outcomes)) report <- rep(report, length(outcomes))
    report[which(names(outcomes)=="null")] <- ""
    if (!("null" %in% names(delays))) {
      delays$null <- c(1, rep(0, length(delays[[1]])-1))
    }
  } else if (!all(rowSums(outcomes)==1)) {
    report <- c(rep(report, length(outcomes)), "")
    outcomes$null <- 1-rowSums(outcomes)
    delays$null <- c(1, rep(0, length(delays[[1]])-1))
  }
  nrow <- length(outcomes)
  list(
    source = src, type="multinomial", names=names(outcomes), report = report,
    prob = t(as.matrix(outcomes)), delays = t(as.matrix(delays))
  )
}

### generate kinetic matrix ----------------------------------------------------------
fun_seir_agestr_kinmatr=function(infect_vartype,vartype_list,n_age,infect_first_ord_pars,y_val){
  dim_sys=n_age*length(vartype_list); K_m=matrix(0,dim_sys,dim_sys)
for (i_age in 1:n_age){
  for (j_varname in infect_vartype) {
    row_ind=fun_sub2ind_seir_agestr(j_age=i_age,varname=j_varname,vartype_list,length(vartype_list),n_age=n_age)
    # outflow terms in diagonal
    outflow_par=infect_first_ord_pars[infect_vartype %in% j_varname]
    K_m[row_ind,row_ind]=-outflow_par
    # inflows: only Ip, Ic, Is have inflows, E--->Ip, Ip->Ic, E->Is
    # row is target, column is source variable
    if (j_varname %in% "E") {K_m[row_ind+1,row_ind]=outflow_par*y_val[i_age]; K_m[row_ind+3,row_ind]=outflow_par*(1-y_val[i_age])}
    if (j_varname %in% "Ip") {K_m[row_ind+1,row_ind]=outflow_par}
  }
}
  K_m
}

### create df with cumul rates
fcn_covidm_df <- function(lmic_simul,sel_vars,pars){   # c("S","E","Ip","Is","R")
  agegr_pop=lmic_simul[lmic_simul$t==0 & lmic_simul$compartment %in% sel_vars] %>% group_by(population,group) %>% summarise(pop=sum(value))
    agegr_pop <- agegr_pop %>% mutate(pop_perc=pop/sum(pop)); colnames(lmic_simul)[colnames(lmic_simul) %in% "value"]="number"
    lmic_simul=left_join(lmic_simul,agegr_pop,by=c("group","population")); 
    lmic_simul[,"percent_agegr_popul"]=round(lmic_simul$number/lmic_simul$pop,4)
    # fraction of the compartment at each timepoint
    lmic_simul=lmic_simul %>% group_by(t,compartment,population) %>% mutate(percent_cases_at_t=number/sum(number))
    # cumulative symptomatic attack rate. distribution of E between Is and Ia set by age-dependent parameter 
    age_sympt_fract=data.frame(group=unique(lmic_simul$group),sympt_fract=pars$pop[[1]]$y)
    lmic_simul_cumulcases=lmic_simul[lmic_simul$compartment %in% "R",]
    lmic_simul_cumulcases=left_join(lmic_simul_cumulcases,age_sympt_fract,by="group")
    for (k in c("number","percent_agegr_popul","percent_cases_at_t")) {
      lmic_simul_cumulcases[,k]=lmic_simul_cumulcases[,k]*lmic_simul_cumulcases[,"sympt_fract"]}
    lmic_simul_cumulcases$compartment="cumul_sympt_cases"
    if (!any(grepl("cumul",colnames(lmic_simul)))) {
      lmic_simul=rbind(lmic_simul,lmic_simul_cumulcases[,!colnames(lmic_simul_cumulcases) %in% "sympt_fract"]) }
    lmic_simul <- lmic_simul %>% pivot_longer(cols=c(number,percent_agegr_popul,percent_cases_at_t),names_to="variable")
    subset(lmic_simul,compartment %in% c(sel_vars,"cumul_sympt_cases"))
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### function to process covidm simul output -----------------
fcn_covidm_process_output <- function(run_dynamics,filter_vars,compartm_types_val,dynamics_type_val,populval,modelparams){
df_sim=subset(run_dynamics,!compartment %in% filter_vars)%>% select(!run) %>% group_by(t,compartment) %>% summarise(value=sum(value))
bind_rows(df_sim, subset(df_sim,compartment %in% "death_o") %>% 
        group_by(compartment) %>% mutate(value=cumsum(value),compartment="D") ) %>% mutate(value_perc_pop=value/populval,
        compartm_type=sapply(compartment,function(k) names(compartm_types_val)[sapply(compartm_types_val,function(x) sum(x %in% k))==1]),
        dynam_type=sapply(compartment,function(k) names(dynamics_type_val)[sapply(dynamics_type_val,function(x) sum(x %in% k))==1]) ) %>%
  mutate(date=seq(as.Date(modelparams$date0),as.Date(modelparams$time1),1)[t+1]) %>%  ungroup()
}

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### plot single simul and calc error -----------------
fcn_covidm_singlesim_error <- function(covidmsimul,df_data,fitting_dates){
  df_data_simul_deaths=left_join(subset(covidmsimul %>% rename(simul=value),compartment %in% "death_o")%>% select(date,simul),
        df_data %>% select(date,daily_baseline_subtr) %>% rename(data=daily_baseline_subtr)) %>% pivot_longer(!date)
abs_err=subset(df_data_simul_deaths, date >= fitting_dates[1] & date <= fitting_dates[2]) %>% group_by(date) %>%
  summarise(data=value[name=="data"],diff=value[name=="data"]-value[name=="simul"],sumval=sum(value),mae=abs(diff))
ggplot(df_data_simul_deaths,aes(x=date)) + geom_line(aes(y=value,group=name,color=name)) + theme_bw() + standard_theme + 
  ggtitle(paste0("MAE=",round(mean(abs_err$mae),2))) + scale_x_date(limits = fitting_dates)
}

### generate matrices to scale force of infection terms ---------
fun_force_inf_vects=function(vartype_list,forceinf_vars,n_age,f_val,N_tot){
b_m_full=matrix(0,length(vartype_list)*n_age,n_age)
for (i_age in 1:n_age) {for (k_var in 1:length(forceinf_vars)){
  b_m_full[(i_age-1)*length(vartype_list) + 1,i_age]=-1;b_m_full[(i_age-1)*length(vartype_list) + 2,i_age]=1
}}
# for summing infectious vars
a_m=matrix(0,n_age,length(infect_vartype)*n_age)
for (k in 1:n_age){  col_lims=((k-1)*length(infect_vartype)+1):(k*length(infect_vartype)); a_m[k,col_lims]=c(0,1,1,f_val[k])/N_tot[k] }
list(b_m_full,a_m)
}

### generate indices of inf and noninf variables
fun_inds_vartypes=function(n_age,vartype_list,infect_vartype){
  ind_all_inf_vars=unlist(lapply(1:n_age, 
      function(x){fun_sub2ind_seir_agestr(j_age=x,varname=infect_vartype,varname_list=vartype_list,
                                          n_var=length(vartype_list),n_age=n_age)}))
ind_all_susceptibles=unlist(lapply(1:n_age, 
      function(x){fun_sub2ind_seir_agestr(j_age=x,varname="S",varname_list=vartype_list,n_var=length(vartype_list),n_age=n_age)}))
list(ind_all_inf_vars,ind_all_susceptibles)
}

### set up plotting theme ------------
standard_theme=theme(# panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
    plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                     legend.title=element_text(size=14),legend.text=element_text(size=12),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))
### colors ----
library(wesanderson) # install.packages("wesanderson")
pal <- wes_palette("Zissou1", 100, type = "continuous")

### seir model without age struct
seir_ode <- function(t,X,parms){ 
  # params=list(K_m,b_vect,c_vect,susc_ind,inf_vars_inds)
  K_m=parms[[1]]; b=parms[[2]]; c_vect=parms[[3]]; susc_ind=parms[[4]]; inf_vars_inds=parms[[5]]
  dXdt=b%*%c_vect%*%X[inf_vars_inds]*X[susc_ind] + K_m%*%X; list(dXdt) }

### process deSolve output
fcn_proc_ode_output<-function(ode_solution,full_varname_list,ind_all_inf_vars){
df_ode_solution=ode_solution %>% as.data.frame() %>% setNames(c("t",full_varname_list))
# tidy format
df_ode_solution_tidy=df_ode_solution[,colSums(df_ode_solution)>0] %>% pivot_longer(!t)
df_ode_solution_tidy[,"vartype"]='noninf'; df_ode_solution_tidy[!grepl("S_",df_ode_solution_tidy$name),"vartype"]='inf'  
df_ode_solution_tidy[,"compartm"]=sapply(strsplit(as.character(df_ode_solution_tidy$name),"_"),"[[",1)
df_ode_solution_tidy[,"agegroup"]=sapply(strsplit(as.character(df_ode_solution_tidy$name),"_"),"[[",2)
df_ode_solution_tidy$agegroup=as.numeric(df_ode_solution_tidy$agegroup)
# fractional values as susceptible popul of age group
df_ode_solution_tidy[,"fract_value"]=df_ode_solution_tidy$value
for (k in 1:n_age) {df_ode_solution_tidy$fract_value[df_ode_solution_tidy$agegroup==k]=
    df_ode_solution_tidy$value[df_ode_solution_tidy$agegroup==k]/N_tot[k]}
# output
list(df_ode_solution,df_ode_solution_tidy) }

### set init conds
fcn_set_init_conds<-function(inf_initval,init_inf_age_groups,init_inf_vartype,n_age,N_tot,vartype_list,ind_all_susceptibles){
ind_inf_initvals=fun_sub2ind_seir_agestr(j_age=init_inf_age_groups,varname=init_inf_vartype,
                                         varname_list=vartype_list,n_var=length(vartype_list),n_age=n_age)
initvals_seir_model=matrix(0,length(vartype_list)*n_age,1); initvals_seir_model[ind_all_susceptibles]=N_tot
# inf_initval=10; 
initvals_seir_model[ind_inf_initvals]=inf_initval; initvals_m=matrix(initvals_seir_model,ncol=n_age)
initvals_seir_model[ind_all_susceptibles]=initvals_seir_model[ind_all_susceptibles]-colSums(initvals_m[2:nrow(initvals_m),])
initvals_seir_model
}

### clinical fraction interpolation
#' calculate probability of clinical disease by age group
#' following posterior estimates by Davies et al
getClinicalFraction <- function(age_groups){
  #' smoothly interpolate between points (x0, y0) and (x1, y1) using cosine interpolation.
  #' for x < x0, returns y0; for x > x1, returns y1; for x0 < x < x1, returns the cosine interpolation between y0 and y1
  interpolate_cos = function(x, x0, y0, x1, y1)
  {     ifelse(x < x0, y0, ifelse(x > x1, y1, y0 + (y1 - y0) * (0.5 - 0.5 * cos(pi * (x - x0) / (x1 - x0)))))   }
  age_groups[, mid := mean(c(age_low, age_high)), by=seq_len(nrow(age_groups))]
  age_y = covid_parameters[["clinical_fraction"]][["values"]][["age_y"]]
  age_m = covid_parameters[["clinical_fraction"]][["values"]][["age_m"]]
  age_o = covid_parameters[["clinical_fraction"]][["values"]][["age_o"]]
  # definition of "young", "middle", and "old"
  young=interpolate_cos(age_groups[, mid], age_y, 1, age_m, 0); old=interpolate_cos(age_groups[, mid], age_m, 0, age_o, 1);
  middle=1-young-old; 
  symp_y=covid_parameters[["clinical_fraction"]][["values"]][["symp_y"]]
  symp_m=covid_parameters[["clinical_fraction"]][["values"]][["symp_m"]]
  symp_o=covid_parameters[["clinical_fraction"]][["values"]][["symp_o"]]
  return(young * symp_y + middle * symp_m + old * symp_o)
}

### age structure of a country --------------------
fun_cntr_agestr=function(i_cntr,i_year,age_groups){
  if (!any((.packages()) %in% "wpp2019")) {library(wpp2019)}; if (!exists("popF")) {data("pop")}
  if (!any(colnames(popF) %in% i_year)) {i_year=colnames(popF)[ncol(popF)]; message("age groups from year ",i_year)}
  
  cntr_agestr=data.frame(agegroups=popF[popF$name %in% i_cntr,"age"],values=popF[popF$name %in% i_cntr,i_year] +
                                        popM[popM$name %in% i_cntr,i_year])
agegr_truthvals=sapply(strsplit(as.character(cntr_agestr$agegroups),"-"),"[[",1) %in% age_groups$age_low
N_tot=cntr_agestr$values[agegr_truthvals]
N_tot[length(N_tot)]=N_tot[length(N_tot)]+sum(cntr_agestr$values[!agegr_truthvals])
N_tot=N_tot*1e3; N_tot
}

### labels for SEIR dynamic plot ---------------------
fun_labels_table=function(df_ode_solution_tidy,age_groups){
df_age_groups=df_ode_solution_tidy %>% group_by(agegroup,vartype) %>% summarise(max_val=max(value),
          min_val=min(value),opt_val=(max(value)+min(value))/2,opt_fract_val=(max(fract_value)+min(fract_value))/2)
df_age_groups[,"agegroup_str"]=paste(age_groups$age_low[df_age_groups$agegroup],'-',
                                     age_groups$age_high[df_age_groups$agegroup],sep="")
df_age_groups}

### interpolate suscept as fcn of age ----------------
fun_interp_suscept=function(suscept_mean_age,suscept_vals,age_groups,N_tot){
interpol_susc=data.frame(approx(suscept_mean_age,suscept_vals))
interpol_susc$y[sapply(age_groups$mid, function(n_midage) {which.min(abs(interpol_susc$x - n_midage))})]/N_tot}

### approximate clinical fraction by "S" shape ------
# min_val=0.04; max_val=0.65; rep_min=5; rep_max=3; 
fun_lin_approx_agedep_par=function(agegroups,min_val,max_val,rep_min,rep_max){
n_intersteps=nrow(agegroups)-(rep_max+rep_min)+1; slope=(max_val-min_val)/n_intersteps
c(rep(min_val,rep_min),min_val+(1:(n_intersteps-1))*slope,rep(max_val,rep_max))
# lines(age_groups$mid,clin_fract_approx)
}

### plot age-dependent params ----------------
fun_plot_agedep_params=function(age_groups,u_val,y_val,f_val,N_tot,plot_flag){
df_agedep_pars=data.frame(age_group=paste(age_groups$age_low,'-',age_groups$age_high,sep=""),suscept=u_val,
  clin_fract=y_val,infectness=f_val) %>% pivot_longer(!age_group)
df_agedep_pars$age_group=factor(df_agedep_pars$age_group,levels=unique(df_agedep_pars$age_group))
if (nchar(plot_flag)>0){
  df_agedep_pars=data.frame(age_group=paste(age_groups$age_low,'-',age_groups$age_high,sep=""),suscept=u_val*sum(N_tot),
                            clin_fract=y_val,infectness=f_val) %>% pivot_longer(!age_group)
g<-ggplot(df_agedep_pars,aes(x=age_group,y=value,group=name,color=name)) + geom_line() + theme_bw() + standard_theme
g} else {
  df_agedep_pars }
}

### HIT vals calcul ----------------
fun_attackrate_df_calc=function(df_ode_solution,N_tot,f_scale,y_val){
attackrate_vals=data.frame(t(1 - round(df_ode_solution[nrow(df_ode_solution),grepl('S_',colnames(df_ode_solution))]/N_tot,3)))
attackrate_vals[,"age_group"]=factor(rownames(attackrate_vals),levels=rownames(attackrate_vals)); rownames(attackrate_vals)=c(); 
attackrate_vals=attackrate_vals[,c(2,1)]; colnames(attackrate_vals)[2]="sum_attackrate"
attackrate_vals[,"clinical_attackrate"]=attackrate_vals$sum_attackrate*y_val; attackrate_vals[,"par_scale"]=f_scale
attackrate_vals}

### function final values of a simul -----------------
fun_final_vals=function(df_ode_solution,N_tot,abs_or_fract){
  if (grepl("fract",abs_or_fract)){
  round(df_ode_solution[nrow(df_ode_solution),setdiff(which(round(df_ode_solution[nrow(df_ode_solution),],2)>0),1)]/N_tot,2)} else {
    round(df_ode_solution[nrow(df_ode_solution),setdiff(which(round(df_ode_solution[nrow(df_ode_solution),],2)>0),1)],2)} 
}

### function paramtable no repetitions
fun_paramtable_norep=function(age_dep_paramtable,scan_parameter){
age_dep_paramtable_clean=age_dep_paramtable[!(!grepl(scan_parameter,age_dep_paramtable$name) & age_dep_paramtable$par_scale!=1),]
susc_truthvals=grepl("suscept",age_dep_paramtable_clean$name)
age_dep_paramtable_clean$value[susc_truthvals]=age_dep_paramtable_clean$value[susc_truthvals]
age_dep_paramtable_clean}

### param scan table -------------------------
fun_paramtable_susc_clinfract_inf=function(scan_parameter,n_loop,midpoint_susc,delta_susc,midpoint_clinfract,delta_clinfr,
                                           mean_susc_exp,mean_clinfract,n_age){
for (par_scale_cnt in 0:n_loop) {
  # SUSCEPTIBILITY
  if (grepl("suscept",scan_parameter)) {par_scale=par_scale_cnt/n_loop} else {par_scale=1}
  min_val_susc=midpoint_susc-par_scale*delta_susc
  u_val=fun_lin_approx_agedep_par(min_val_susc,midpoint_susc+par_scale*delta_susc,rep_min=3,rep_max=10)/sum(N_tot)
  u_val=u_val*(mean_susc_exp/mean(u_val))
  # CLINICAL FRACTION
  # y_val=0.04+(1:n_age)*0.6/n_age # from "Age-dependent" getClinicalFraction(age_groups) # f_scale needs to be 0<=x<=1
  if (grepl("clin.*fract",scan_parameter)) {par_scale=par_scale_cnt/n_loop} else {par_scale=1}
  min_val_clinfr=midpoint_clinfract-par_scale*delta_clinfr
  y_val=fun_lin_approx_agedep_par(min_val_clinfr,midpoint_clinfract+par_scale*delta_clinfr,rep_min=5,rep_max=3)
  y_val=y_val*(mean_clinfract/mean(y_val))
  # INFECTIOUSNESS
  if (grepl("infect",scan_parameter)) {f_val=rep(0.5+(par_scale_cnt-n_loop/2)*(0.25/(n_loop/2)),n_age)} else {f_val=rep(0.5,n_age)}
  # plot params: fun_plot_agedep_params(age_groups,u_val,y_val,f_val,N_tot,'plot')
  # store param table
  if (par_scale_cnt==0) {
   age_dep_paramtable=cbind(fun_plot_agedep_params(age_groups,u_val,y_val,f_val,N_tot,''),par_scale_cnt)} else {
   age_dep_paramtable=rbind(age_dep_paramtable,cbind(fun_plot_agedep_params(age_groups,u_val,y_val,f_val,N_tot,''),par_scale_cnt))}
  # LOOP progress
  # print(paste("u (suscept)=[",u_val[1],u_val[length(y_val)],"] mean=",mean(u_val),"|", 
  #             "y (clin fract) =[",round(y_val[1],3),round(y_val[length(y_val)],3),"] mean=",mean(y_val),"|",
  #             "f (infect asympt) =[",round(f_val[1],3),round(f_val[length(f_val)],3),"] mean=",mean(f_val), sep=" "))
}
  age_dep_paramtable
}

### function mean age of country  ---------------------
fun_meanage_cntr=function(popM,popF,countryval,year_str){
full_agestr=data.frame(agegroups=popF[popF$name %in% countryval,"age"],values=popF[popF$name %in% countryval,year_str] +
popM[popM$name %in% countryval,year_str])
mean_age_groups=c(colMeans(sapply(strsplit(as.character(full_agestr$agegroups[1:(nrow(full_agestr)-1)]),"-"),as.numeric)),100)
sum((full_agestr$values/sum(full_agestr$values))*mean_age_groups)}

### function plot dyn by age groups ---------------------
fun_seir_agegroups_dyn=function(age_groups,df_ode_solution_tidy,xlim_val,abs_or_fract,timesteps,standard_theme){
# df_ode_solution_tidy=l_proc_sol[[2]]; xlim_val=365; 
df_age_groups=fun_labels_table(df_ode_solution_tidy,age_groups)
# absolute or % value?
# abs_or_fract=2; 
if (abs_or_fract==1){value_col="value";label_col="opt_val";savetag="absval"} else {
  value_col="fract_value"; label_col="opt_fract_val";savetag="fractval"}
g<-ggplot(df_ode_solution_tidy,aes_string(x="t",y=value_col,group="name",color="compartm")) + geom_line(size=1.05) + 
  facet_wrap(~agegroup+vartype,scales='free_y',ncol=4) +
  geom_text(data=df_age_groups,aes_string(x=8,y=label_col,label="agegroup_str",group=NULL),size=3,color="black") +
  theme_bw()+standard_theme+theme(strip.background=element_blank(),strip.text=element_blank(),axis.text.y=element_text(size=6)) +
  labs(linetype='vars',color='vars') + scale_x_continuous(breaks=seq(timesteps[1],xlim_val,10),limits=c(0,xlim_val)) + 
  xlab('days') + ylab('') + ggtitle(paste0(paste0(vartype_list,collapse="-"),'-R simulation'))
g
# SAVE
}

#### parameter scan for one parameter & country ----------------------
# fun_paramscan_singlepar=function(age_dep_paramtable,scan_parameter,infect_first_ord_pars,infect_vartype,vartype_list,
#                                  C_m, inf_initval,init_inf_age_groups,init_inf_vartype,N_tot,n_years,age_groups){
#   n_age=nrow(C_m)
#   scan_param=c("u_val","y_val","f_val")[sapply(c("susc","clin","infect"), function(x) {grepl(x,scan_parameter)})]
#   midvalue_agedep=round(median(unique(age_dep_paramtable$dep_fact)))
#   age_dep_paramtable_red=age_dep_paramtable[(age_dep_paramtable$name %in% scan_param) | 
#                                     ((!age_dep_paramtable$name %in% scan_param) & age_dep_paramtable$dep_fact==midvalue_agedep),]
#   for (age_dep_val in unique(age_dep_paramtable$dep_fact)){
#     if (sum(unique(age_dep_paramtable_red$dep_fact[age_dep_paramtable_red$name %in% "u_val"]) %in% age_dep_val)){
#       k_u=age_dep_val} else {k_u=midvalue_agedep}
#     if (sum(unique(age_dep_paramtable_red$dep_fact[age_dep_paramtable_red$name %in% "y_val"]) %in% age_dep_val)){
#       k_y=age_dep_val} else {k_y=midvalue_agedep}
#     if (sum(unique(age_dep_paramtable_red$dep_fact[age_dep_paramtable_red$name %in% "f_val"]) %in% age_dep_val)){
#       k_f=age_dep_val} else {k_f=midvalue_agedep}
#   u_val=age_dep_paramtable$value[age_dep_paramtable$dep_fact==k_u & age_dep_paramtable$name %in% "u_val"]
#   y_val=age_dep_paramtable$value[age_dep_paramtable$dep_fact==k_y & age_dep_paramtable$name %in% "y_val"]
#   f_val=age_dep_paramtable$value[age_dep_paramtable$dep_fact==k_f & age_dep_paramtable$name %in% "f_val"]
#   # constant transmission parameters
#   # d_e=1/4; d_p=1/1.5; d_c=1/3.5; d_s=1/5; infect_first_ord_pars=c(d_e,d_p,d_c,d_s)
#   # KINETIC MATRIX (LINEAR TERMS)
#   K_m=fun_seir_agestr_kinmatr(infect_vartype,vartype_list,n_age,infect_first_ord_pars,y_val)
#   # vectors to scale FORCE of INFECTION terms
#   list_bm_am=fun_force_inf_vects(vartype_list,forceinf_vars=c("S","E"),n_age,f_val,N_tot); b_m_full=list_bm_am[[1]]; a_m=list_bm_am[[2]]
#   # inf vector is: b_m_full %*% diag(c(1,2)) %*% C_m %*% diag(u_val) %*% a_m %*% matrix(1,8,1)
#   l=fun_inds_vartypes(n_age,vartype_list,infect_vartype); ind_all_inf_vars=l[[1]]; ind_all_susceptibles=l[[2]]
#   # INITIAL CONDITIONS
#   # seed epidemic by "E">0 in a given (or multiple) age groups
#   initvals_seir_model=fcn_set_init_conds(inf_initval,init_inf_age_groups,init_inf_vartype,n_age,N_tot,
#                                          vartype_list,ind_all_susceptibles)
#   ### run ODEs ----------------------------
#   # time resolution and duration
#   # n_days_year=365; n_years=1; 
#   max_time=n_years*365; timesteps <- seq(0,max_time,by=0.25)
#   # define inputs for ODE fcn
#   params=list(K_m, C_m, b_m_full, u_val, a_m, ind_all_susceptibles, ind_all_inf_vars)
#   # RUN # 
#   ode_solution<-lsoda(initvals_seir_model,timesteps,func=seir_agestr_ode,parms=params)
#   print(paste("u (susc)=[",round(u_val[1]/(10^floor(log10(u_val[1]))),2),"e",floor(log10(u_val[1])), 
#    round(u_val[length(u_val)]/(10^floor(log10(u_val[length(u_val)]))),2),"e",floor(log10(u_val[length(u_val)])),"] mean=",
#    round(mean(u_val)/(10^floor(log10(mean(u_val)))),2),"e",floor(log10(mean(u_val))), "|",
#               "y (clinfrac) =[",round(y_val[1],3),round(y_val[length(y_val)],3),"] mean=",round(mean(y_val),2),"|",
#               "f (inf asym) =[",round(f_val[1],3),round(f_val[length(f_val)],3),"] mean=",round(mean(f_val),2), sep=" "))
#   # PROCESS OUTPUT
#   l_proc_sol=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars); df_ode_solution=l_proc_sol[[1]]
#   # infected fraction of popul (HIT)
#   attack_rate=fun_attackrate_df_calc(df_ode_solution,N_tot,age_dep_val,y_val)
#   # print(attack_rate)
#   if (age_dep_val==0) {scan_table=attack_rate} else {scan_table=rbind(scan_table,attack_rate)}
#   } # end of for loop
#   scan_table[,"agegroup_names"]=factor(paste(age_groups$age_low,age_groups$age_high,sep = "-")[as.numeric(sapply(strsplit(as.character(
#     scan_table$age_group),"_"),"[[",2))],levels=paste(age_groups$age_low,age_groups$age_high,sep = "-"))
#   colorvar="col_scale"; scan_table[,colorvar]=factor(scan_table$par_scale)
#   scan_table
# }

### param scan for single parameter with COVIDM

fun_paramscan_singlepar=function(countryval,age_dep_paramtable,scan_parameter,sel_vars,n_years,age_groups){
  params=cm_parameters_SEI3R(gsub(countryval,"Sudan|Somalia","Ethiopia"))
  N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups)
  params$pop[[1]]$name=countryval; params$pop[[1]]$size=N_tot
  scan_param=c("u_val","y_val","f_val")[sapply(c("susc","clin","infect"), function(x) {grepl(x,scan_parameter)})]
  midvalue_agedep=round(median(unique(age_dep_paramtable$dep_fact)))
  age_dep_paramtable_red=age_dep_paramtable[(age_dep_paramtable$name %in% scan_param) | 
                                              ((!age_dep_paramtable$name %in% scan_param) & age_dep_paramtable$dep_fact==midvalue_agedep),]
  for (age_dep_val in unique(age_dep_paramtable$dep_fact)){
    if (sum(unique(age_dep_paramtable_red$dep_fact[age_dep_paramtable_red$name %in% "u_val"]) %in% age_dep_val)){
      k_u=age_dep_val} else {k_u=midvalue_agedep}
    if (sum(unique(age_dep_paramtable_red$dep_fact[age_dep_paramtable_red$name %in% "y_val"]) %in% age_dep_val)){
      k_y=age_dep_val} else {k_y=midvalue_agedep}
    if (sum(unique(age_dep_paramtable_red$dep_fact[age_dep_paramtable_red$name %in% "f_val"]) %in% age_dep_val)){
      k_f=age_dep_val} else {k_f=midvalue_agedep}
    params$pop[[1]]$u=age_dep_paramtable$value[age_dep_paramtable$dep_fact==k_u & age_dep_paramtable$name %in% "u_val"]
    params$pop[[1]]$y=age_dep_paramtable$value[age_dep_paramtable$dep_fact==k_y & age_dep_paramtable$name %in% "y_val"]
    params$pop[[1]]$fIa=age_dep_paramtable$value[age_dep_paramtable$dep_fact==k_f & age_dep_paramtable$name %in% "f_val"]
    ### run ODEs ----------------------------
    # time resolution and duration
    # n_days_year=365; n_years=1; 
    run=cm_simulate(params,1); lmic_simul=run$dynamics; lmic_simul=fcn_covidm_df(lmic_simul,sel_vars,params)
    print(paste("u (susc)=[",round(u_val[1]/(10^floor(log10(u_val[1]))),2),"e",floor(log10(u_val[1])), 
                round(u_val[length(u_val)]/(10^floor(log10(u_val[length(u_val)]))),2),"e",floor(log10(u_val[length(u_val)])),"] mean=",
                round(mean(u_val)/(10^floor(log10(mean(u_val)))),2),"e",floor(log10(mean(u_val))), "|",
                "y (clinfrac) =[",round(y_val[1],3),round(y_val[length(y_val)],3),"] mean=",round(mean(y_val),2),"|",
                "f (inf asym) =[",round(f_val[1],3),round(f_val[length(f_val)],3),"] mean=",round(mean(f_val),2), sep=" "))
    # PROCESS OUTPUT
    attack_rate=(subset(lmic_simul,t==max(lmic_simul$t) & compartment %in% c("cumul_sympt_cases","R") & variable %in% "percent_agegr_popul") %>%
                   pivot_wider(names_from=compartment,values_from=value))[,c("R","cumul_sympt_cases","group")]
    colnames(attack_rate) = c("sum_attackrate","clinical_attackrate","agegroup_names")
    attack_rate[,"par_scale"]=age_dep_val; attack_rate[,"age_group"]=paste0("S_",as.numeric(attack_rate$agegroup_names))
    attack_rate[,"col_scale"]=factor(attack_rate$par_scale)
    attack_rate=attack_rate[,c("age_group","sum_attackrate","clinical_attackrate","par_scale","agegroup_names","col_scale")]
    if (age_dep_val==0) {scan_table=attack_rate} else {scan_table=rbind(scan_table,attack_rate)}
  } # end of for loop
  # scan_table[,"agegroup_names"]=factor(paste(age_groups$age_low,age_groups$age_high,sep = "-")[as.numeric(sapply(strsplit(as.character(
  #   scan_table$age_group),"_"),"[[",2))],levels=paste(age_groups$age_low,age_groups$age_high,sep = "-"))
  # colorvar="col_scale"; scan_table[,colorvar]=factor(scan_table$par_scale)
  scan_table
}

### parameter scan for multiple parameters and countries with COVIDM --------------------------------
fun_onedim_multipar_multicntr_scan=function(cntr_list,year_str,age_groups,scan_param_fullname,n_loop,
                                                   agedep_param_lims,sel_vars,n_years){
  g(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,
    delta_clinfr,rep_min_clinfr,rep_max_clinfr) %=% agedep_param_lims
  age_dep_paramtable=fun_agedep_partable(agedep_param_lims,n_loop)
  for (k_cntr in 1:length(cntr_list)){
    countryval=cntr_list[k_cntr]; N_tot=fun_cntr_agestr(countryval,i_year=year_str,age_groups) 
    # uniform: N_tot=rep(1e6,n_age)
    # CONTACT MATRIX
    covid_params=cm_parameters_SEI3R(gsub("Sudan|Somalia","Ethiopia",countryval))
    params$pop[[1]]$name=countryval; params$pop[[1]]$size=N_tot
    C_m=Reduce('+',covid_params$pop[[1]]$matrices); C_m=fun_contmatr_recipr(C_m,N_tot); n_age=nrow(C_m)
    numchars=paste0(as.character(round(c(sum(N_tot)/1e6,median(C_m),max(Re(eigen(C_m)$values))),1)),collapse = "," )
    print(paste0(countryval," [popul,median(C_m),max(eigval(C_m))] = [",numchars,"]"))
    # scan in parameter
    # scan_param_fullname=c("susceptibility","clinical fraction","asymptomatic infectiousness")
    # c("u_val","y_val","f_val")
    for (k_par in 1:length(scan_param_fullname)){ # n_loop=6
      attack_rates_scan=fun_paramscan_singlepar(countryval,age_dep_paramtable,scan_param_fullname[k_par],sel_vars,n_years,age_groups)
      attack_rates_scan[,"scanpar"]=scan_param_fullname[k_par]; 
      if (k_par==1){attack_rates_scan_multipar=attack_rates_scan} else {
        attack_rates_scan_multipar=rbind(attack_rates_scan_multipar,attack_rates_scan)} }
    # long format
    attack_rates_scan_multipar[,"country"]=countryval; att_colns=colnames(attack_rates_scan_multipar)
    attack_rates_scan_multipar=attack_rates_scan_multipar %>% pivot_longer(!att_colns[!grepl("attack",att_colns)],values_to="fract_agegroup")
    # case numbers
    attack_rates_scan_multipar[,"agegr_pop"]=rep(as.vector(t(matrix(rep(N_tot,2),ncol=2))),
                                                 nrow(attack_rates_scan_multipar)/length(as.vector(t(matrix(rep(N_tot,2),ncol=2)))))
    attack_rates_scan_multipar[,"n_case"]=round(attack_rates_scan_multipar$agegr_pop*attack_rates_scan_multipar$fract_agegroup,1)
    attack_rates_scan_multipar[,"per_1000_popul"]=round(1e3*attack_rates_scan_multipar$n_case/sum(N_tot),1)
    
    if (k_cntr==1) {attack_rates_multipar_multicntr=attack_rates_scan_multipar} else {
      attack_rates_multipar_multicntr=rbind(attack_rates_multipar_multicntr,attack_rates_scan_multipar)}
  } # end of for loop for cntrs
  attack_rates_multipar_multicntr=attack_rates_multipar_multicntr %>% group_by(par_scale,scanpar,name,country) %>% 
    mutate(fraction_cases=n_case/sum(n_case))
  attack_rates_multipar_multicntr
}

### extract process of optimisation
fcn_extract_optimresults <- function(optim_proc,parnames){
  optim_outputs=data.frame(gsub("\\s+", " ", optim_proc[which(grepl("BUILD",optim_proc)):(which(grepl("\\$par",optim_proc))-3)]))
  colnames(optim_outputs)="onecol" # str_replace_all(optim_proc[grepl('\\$',optim_proc)],'\\$','')
  optim_outputs=optim_outputs %>% separate(col=onecol,into=c("message","count","sse0","sse"),sep=' ') # 
  optim_outputs[,2:ncol(optim_outputs)]=sapply(optim_outputs[2:ncol(optim_outputs)],as.numeric)
  optim_outputs
}

### parameter scan for multiple parameters and countries (own model) --------------------------------
# fun_onedim_multipar_multicntr_scan=function(cntr_list,year_str,age_groups,cm_path,scan_param_fullname,n_loop,
#                                      agedep_param_lims,infect_first_ord_pars,infect_vartype,vartype_list,
#                                      inf_initval,init_inf_age_groups,init_inf_var,n_years){
#   g(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,
#     delta_clinfr,rep_min_clinfr,rep_max_clinfr) %=% agedep_param_lims
#   age_dep_paramtable=fun_agedep_partable(agedep_param_lims,n_loop)
# for (k_cntr in 1:length(cntr_list)){countryval=cntr_list[k_cntr]; N_tot=fun_cntr_agestr(countryval,i_year=year_str,age_groups) 
# # uniform: N_tot=rep(1e6,n_age)
# # CONTACT MATRIX
# if (!exists("covid_params")){cm_force_rebuild=F;cm_build_verbose=T;cm_version=2; source(file.path(cm_path,"R","covidm.R"))}
# covid_params=cm_parameters_SEI3R(gsub("Sudan","Ethiopia",countryval))
# C_m=Reduce('+',covid_params$pop[[1]]$matrices); C_m=fun_contmatr_recipr(C_m,N_tot); n_age=nrow(C_m)
# numchars=paste0(as.character(round(c(sum(N_tot)/1e6,median(C_m),max(Re(eigen(C_m)$values))),1)),collapse = "," )
# print(paste0(countryval," [popul,median(C_m),max(eigval(C_m))] = [",numchars,"]"))
# # scan in parameter
# # scan_param_fullname=c("susceptibility","clinical fraction","asymptomatic infectiousness")
# # c("u_val","y_val","f_val")
# for (k_par in 1:length(scan_param_fullname)){ # n_loop=6
# attack_rates_scan=fun_paramscan_singlepar(age_dep_paramtable,scan_parameter=scan_param_fullname[k_par],infect_first_ord_pars,
#                   infect_vartype,vartype_list,C_m,inf_initval,init_inf_age_groups,init_inf_var,N_tot,n_years,age_groups)
# attack_rates_scan[,"scanpar"]=scan_param_fullname[k_par]; if (k_par==1){attack_rates_scan_multipar=attack_rates_scan} else {
#   attack_rates_scan_multipar=rbind(attack_rates_scan_multipar,attack_rates_scan)} }
# # long format
# attack_rates_scan_multipar[,"country"]=countryval; att_colns=colnames(attack_rates_scan_multipar)
# attack_rates_scan_multipar=attack_rates_scan_multipar %>% pivot_longer(!att_colns[!grepl("attack",att_colns)],values_to="fract_agegroup")
# # case numbers
# attack_rates_scan_multipar[,"agegr_pop"]=rep(as.vector(t(matrix(rep(N_tot,2),ncol=2))),
#                                     nrow(attack_rates_scan_multipar)/length(as.vector(t(matrix(rep(N_tot,2),ncol=2)))))
# attack_rates_scan_multipar[,"n_case"]=round(attack_rates_scan_multipar$agegr_pop*attack_rates_scan_multipar$fract_agegroup,1)
# attack_rates_scan_multipar[,"per_1000_popul"]=round(1e3*attack_rates_scan_multipar$n_case/sum(N_tot),1)
# 
# if (k_cntr==1) {attack_rates_multipar_multicntr=attack_rates_scan_multipar} else {
#   attack_rates_multipar_multicntr=rbind(attack_rates_multipar_multicntr,attack_rates_scan_multipar)}
# } # end of for loop for cntrs
#   attack_rates_multipar_multicntr=attack_rates_multipar_multicntr %>% group_by(par_scale,scanpar,name,country) %>% 
#     mutate(fraction_cases=n_case/sum(n_case))
#   attack_rates_multipar_multicntr
# }

### calculate NGM and R0 -----------------
fun_NGM_R0=function(C_m,N_tot,y_val,u_val,f_val,lin_rates){
NGM=matrix(0,nrow=length(N_tot),ncol=length(N_tot)); g(d_p,d_c,d_s) %=% lin_rates
for (i_row in 1:length(y_val)){ for (j_col in 1:length(y_val)){
  # agegr_frac=N_tot[i_row]/sum(N_tot) # 1
  NGM[i_row,j_col]=u_val[i_row]*C_m[i_row,j_col]*(y_val[j_col]*(1/d_p+1/d_c) + (1-y_val[j_col])*unique(f_val)/d_s) } }
list(NGM,max(Re(eigen(NGM)$values)))
}

### make contact matrix symmetric
fun_contmatr_recipr=function(C_m_orig,N_tot){ C_m=C_m_orig # matrix(0,nrow=nrow(C_m_orig),ncol=ncol(C_m_orig))
matr_combs=permutations(n=nrow(C_m_orig),r=2,repeats.allowed=T)
for (k_matr in 1:nrow(matr_combs)) {
  C_m[matr_combs[k_matr,1],matr_combs[k_matr,2]]=(C_m_orig[matr_combs[k_matr,1],matr_combs[k_matr,2]] + 
 C_m_orig[matr_combs[k_matr,2],matr_combs[k_matr,1]]*(N_tot[matr_combs[k_matr,2]]/N_tot[matr_combs[k_matr,1]]) )/2 }
C_m
}

### agedep paramtable for multidim scan -----------------
fun_agedep_partable=function(agedep_param_lims,k_max){
  g(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,
    delta_clinfr,rep_min_clinfr,rep_max_clinfr) %=% agedep_param_lims
for (k_par in 0:k_max){ u_val=fun_lin_approx_agedep_par(min_val=midpoint_susc-(k_par/k_max)*delta_susc,
                                                        max_val=midpoint_susc+(k_par/k_max)*delta_susc,rep_min_susc,rep_max_susc)
  y_val=fun_lin_approx_agedep_par(midpoint_clinfract-(k_par/k_max)*delta_clinfr,
                                midpoint_clinfract+(k_par/k_max)*delta_clinfr,rep_min_clinfr,rep_max_clinfr)
  f_val=rep(k_par/k_max,16)
if (k_par==0){agedep_partable=data.frame(dep_fact=k_par,cbind(1:16,u_val,y_val,f_val))} else{
  agedep_partable=rbind(agedep_partable,data.frame(dep_fact=k_par,cbind(1:16,u_val,y_val,f_val)))} } 
  colnames(agedep_partable)[2]="agegroup"
  agedep_partable %>% pivot_longer(!c(agegroup,dep_fact))
}

### multidimensional scan for R0 -----------------
fun_multidim_scan=function(countryval,C_m,N_tot,k_max,agedep_param_lims,lin_rates){
R0_scan=data.frame(matrix(0,ncol=5,nrow=k_max^3)); k_count=0; colnames(R0_scan)=c("susc","clinfract","asympt_inf","R0","HIT")
g(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,
  delta_clinfr,rep_min_clinfr,rep_max_clinfr) %=% agedep_param_lims
for (k_susc in 0:k_max){
  for (k_clinfrac in 0:k_max){
    for (k_infect in 0:k_max){
      u_val=fun_lin_approx_agedep_par(min_val=midpoint_susc-(k_susc/k_max)*delta_susc,
                                      max_val=midpoint_susc+(k_susc/k_max)*delta_susc,rep_min_susc,rep_max_susc)
      y_val=fun_lin_approx_agedep_par(midpoint_clinfract-(k_clinfrac/k_max)*delta_clinfr,
                                      midpoint_clinfract+(k_clinfrac/k_max)*delta_clinfr,rep_min_clinfr,rep_max_clinfr); 
      f_val=k_infect/k_max; k_count=k_count+1; 
      R0=fun_NGM_R0(C_m,N_tot,y_val,u_val,f_val,lin_rates)[[2]]; R0_scan[k_count,]=c(k_susc,k_clinfrac,k_infect,R0,1-1/R0)
    }
  } 
}
R0_scan[,"asympt_inf_str"]=factor(paste0("asympt infectness=",round(R0_scan$asympt_inf/k_max,2)),
                                  levels=unique(paste0("asympt infectness=",round(R0_scan$asympt_inf/k_max,2))))
R0_scan[,"clinfract_str"]=factor(paste0("clin fract agedep=",round(R0_scan$clinfract/k_max,2)),
                                  levels=unique(paste0("clin fract agedep=",round(R0_scan$clinfract/k_max,2) )))
R0_scan[,"susc_str"]=factor(paste0("suscept agedep=",round(R0_scan$susc/k_max,2)),
                                 levels=unique(paste0("suscept agedep=",round(R0_scan$susc/k_max,2) )))
R0_scan[,"country"]=countryval
R0_scan
}

### ODE fcns for somalia model
# simple_sir_somalia_interpol <- function(t,x,parms) {
#   S<-x[1]; E<-x[2]; I_C<-x[3]; I_CR<-x[4]; I_CD=x[5]; R=x[6]; D=x[7]
#   with(as.list(c(x,parms)), {
#     beta_variable=beta_eff*input_npi(t); force_inf=beta_variable*(I_C + I_CR + sever_inf*I_CD)
#     force_inf=beta_variable*(I_C + I_CR + sever_inf*I_CD)
#     dS = -S*force_inf; dE = S*force_inf - d_e*E + input_seeding(t)*inv_timestep
#     dI_C=d_e*E - k_fast*I_C; dI_CR=(1-IFR_estim)*k_fast*I_C - d_cr*I_CR; dI_CD=IFR_estim*k_fast*I_C - d_cd*I_CD
#     dR=d_cr*I_CR; dD=d_cd*I_CD; cumulinf=d_e*E  # infect who'll recover or die
#     return(list(c(dS,dE,dI_C,dI_CR,dI_CD,dR,dD,cumulinf)))
#   }) }
#
#
# simple_sir_somalia_interpol <- function(t,x,parms) {
#   S<-x[1]; E<-x[2]; I_C<-x[3]; I_CR<-x[4]; I_CD=x[5]; I_S=x[6]; R=x[7]; D=x[8]
#   with(as.list(c(x,parms)), {
#     beta_variable=beta_eff*input_npi(t); force_inf=beta_variable*(I_C + I_CR + sever_inf*I_CD + f_asympt*I_S); d_s=d_cr
#     dS = -S*force_inf; dE = S*force_inf - d_e*E + input_seeding(t)*inv_timestep
#     dI_C=symptshare*d_e*E - k_fast*I_C
#     dI_CR=(1-IFR_estim)*k_fast*I_C - d_cr*I_CR
#     dI_CD=IFR_estim*k_fast*I_C - d_cd*I_CD
#     dI_S=(1-symptshare)*d_e*E - d_s*I_S
#     dR=d_cr*I_CR + d_s*I_S; dD=d_cd*I_CD; cumul_sympt_inf=symptshare*d_e*E; cumul_asympt_inf=(1-symptshare)*d_e*E  
#     return(list(c(dS,dE,dI_C,dI_CR,dI_CD,dI_S,dR,dD,cumul_sympt_inf,cumul_asympt_inf)))
#   }) }

# ODE fcn 
simple_sir_somalia_interpol <- function(t,x,parms) {
  S<-x[1]; E<-x[2]; I_P<-x[3]; I_C<-x[4]; I_S=x[5]; R=x[6]
  with(as.list(c(x,parms)), {
    beta_variable=beta_eff*input_npi(t); force_inf=beta_variable*(I_P + I_C + f_asympt*I_S)
    dS = -S*force_inf; dE = S*force_inf - d_e*E + input_seeding(t)*inv_timestep
    dI_P=symptshare*d_e*E - d_p*I_P; dI_C=d_p*I_P - d_c*I_C
    dI_S=(1-symptshare)*d_e*E - d_s*I_S
    dR=d_c*I_C + d_s*I_S
    cumul_sympt_inf=d_p*I_P; cumul_asympt_inf=(1-symptshare)*d_e*E
    return(list(c(dS,dE,dI_P,dI_C,dI_S,dR,cumul_sympt_inf,cumul_asympt_inf)))
  }) }

### function to run and plot single simul for somalia SIR model -----------------------------
fcn_somal_sir_singlesimul <- function(varcateglist,time_steps,num_params,day0,NPI_input,fcn_flag){
  popul_tot=num_params["popul_tot"]; timespan_dates=NPI_input$date
  initvals_S_I=c(popul_tot,rep(0,length(varcateglist$name_vars)-1)); names(initvals_S_I)=c(varcateglist$name_vars)
  day0_num=as.numeric(day0-timespan_dates[1]); 

  timesteps<-seq(0,length(NPI_input$OxCGRT_scaled),by=unique(diff(time_steps)))
  timevar_signal<-data.frame(t=(1:nrow(NPI_input))-1,
                             npi_index=1-(1-NPI_input$OxCGRT_scaled)*NPI_input$NPI_on*num_params["compliance"],seeding=0)
  timevar_signal$seeding[timevar_signal$t %in% day0_num:(day0_num+num_params["seeding_duration"])]=num_params["seed_size_val"]
  # print(timevar_signal)
input_npi <- approxfun(timevar_signal[,c("t","npi_index")]); input_seeding<-approxfun(timevar_signal[,c("t","seeding")])
### create ODE object
symptshare=num_params["sympt_share"]
params=list(beta_eff=num_params['beta']/popul_tot,d_e=1/num_params['d_e'],d_p=1/num_params['d_p'],d_c=1/num_params['d_c'],
    d_s=1/num_params['d_s'],IFR_estim=num_params['IFR'],seed_size=num_params["seed_size_val"],
    symptshare=num_params["sympt_share"],f_asympt=num_params["asympt_infness"],inv_timestep=1/unique(diff(timesteps)))
# ODE fcn 
ode_fcn <- function(t,x,parms) {
  S<-x[1]; E<-x[2]; I_P<-x[3]; I_C<-x[4]; I_S=x[5]; R=x[6]
  with(as.list(c(x,parms)), {
    beta_variable=beta_eff*input_npi(t); force_inf=beta_variable*(I_P + I_C + f_asympt*I_S)
    dS = -S*force_inf; dE = S*force_inf - d_e*E + input_seeding(t)*inv_timestep
    dI_P=symptshare*d_e*E - d_p*I_P;  dI_C=d_p*I_P - d_c*I_C;  dI_S=(1-symptshare)*d_e*E - d_s*I_S
    dR=d_c*I_C + d_s*I_S
    cumul_sympt_inf=d_p*I_P; cumul_asympt_inf=(1-symptshare)*d_e*E; cumul_presympt_inf=symptshare*d_e*E; incid_E=S*force_inf
    return(list(c(dS,dE,dI_P,dI_C,dI_S,dR,cumul_sympt_inf,cumul_asympt_inf,cumul_presympt_inf,incid_E)))
  }) }
# RUN
  df_ode_solution <- ode(y=initvals_S_I,times=timesteps,func=ode_fcn,parms=params,method="euler") %>% # euler
    as.data.frame() %>% setNames(c("t",varcateglist$name_vars)) # %>% filter(t %% 1 ==0) 
# process output
  if (nchar(fcn_flag)>0) {
  df_ode_solution = df_ode_solution %>% 
    mutate(new_E=incid_E-lag(incid_E,default=incid_E[1]),
      new_sympt_inf=cumul_sympt_inf-lag(cumul_sympt_inf,default=cumul_sympt_inf[1]),
      new_presympt_inf=cumul_presympt_inf-lag(cumul_presympt_inf,default=cumul_presympt_inf[1]),
    new_asympt_inf=cumul_asympt_inf-lag(cumul_asympt_inf,default=cumul_asympt_inf[1]),new_recov=R-lag(R,default=R[1]),
    new_deaths=symptshare*num_params['IFR']*lag(new_recov,n=num_params["d_death"]/unique(diff(timesteps))),
    new_deaths_poiss=rpois(n=1:nrow(df_ode_solution),replace_na(new_deaths,0)),
    D=cumsum(replace_na(new_deaths,0)),R_recov=R-D) %>% select(varcateglist$sel_vars) %>% pivot_longer(cols=!t) %>% 
    mutate(case_death_var=ifelse(name %in% varcateglist$case_vars,"case","fatal"),
      cumul_trans_var=ifelse(name %in% varcateglist$cumul_var,"cumul","transient"),
      state_delta_var=ifelse(name %in% varcateglist$delta_var,"delta","state"),
      date=timespan_dates[t+1],name=factor(name,levels=unique(name)),
      introd_date=day0,seed_size=num_params["seed_size_val"],IFR_sympt=num_params['IFR'], # ,IFR_all_inf=num_params['IFR']*symptshare
      npi_compliance=num_params["compliance"],beta=num_params['beta'],
      R0=num_params['beta']*(symptshare*(num_params['d_c'] + num_params['d_p']) + 
                               (1-symptshare)*num_params["asympt_infness"]*num_params['d_s']),
      seeding_duration=num_params["seeding_duration"]) %>% filter(!(is.na(date) | is.na(value)))
    }
df_ode_solution
}

### load IFRs from imperial estimates
fcn_load_ifr <- function(ifr_filepath){
read_csv(ifr_filepath) %>% group_by(agegroup) %>% 
  summarise(agegroup=unique(agegroup),mean=mean(mean)/1e2,CI95_lower=mean(CI95_lower)/1e2,CI95_upper=mean(CI95_upper)/1e2) %>% 
  mutate(agegroup_min=ifelse(grepl("-|\\+",agegroup),as.numeric(sapply(strsplit(agegroup,"-|\\+"),function(x) {x[1]})),NA ),
         agegroup_max=ifelse(grepl("-",agegroup),as.numeric(sapply(strsplit(agegroup,"-"),function(x) {x[2]})),NA ) ) %>% 
  filter(!is.na(agegroup_min))
}

### load cntr age struct, merge age groups above a cutoff ---------------------------
fcn_load_age_str <- function(cntrval,cutoff_age){
data.frame(agegroup=popF[popF$name %in% cntrval,"age"],
  value=(popF[popF$name %in% cntrval,"2020"] + popM[popM$name %in% cntrval,"2020"])*1e3) %>%  
  mutate(agegroup_min=ifelse(grepl("-|\\+",as.character(agegroup)),
                             as.numeric(sapply(strsplit(as.character(agegroup),"-|\\+"),function(x) {x[1]})),NA),
         agegroup_max=ifelse(grepl("-|\\+",as.character(agegroup)),
                             as.numeric(sapply(strsplit(as.character(agegroup),"-|\\+"),function(x) {x[2]})),NA),
  agegroup=ifelse(agegroup_min<cutoff_age,as.character(agegroup),paste0(cutoff_age,"+"))) %>% group_by(agegroup) %>% 
  summarise(agegroup=unique(agegroup),value=sum(value),agegroup_min=unique(agegroup_min)[1],agegroup_max=max(agegroup_max)) %>%
    mutate(agegroup_max=ifelse(is.na(agegroup_max),agegroup_min+4,agegroup_max),agegroup_mean=(agegroup_min+agegroup_max+1)/2) %>%
    rename(agegroupsize=value) %>% arrange(agegroup_min)
}

### merge IFRs above a given age (typically 75) ---------------------------
fcn_merge_ifr_above_age <- function(age_ifr_data,cutoff_age){
  age_ifr_data %>% mutate(agegroup=ifelse(agegroup_min<cutoff_age,as.character(agegroup),paste0(cutoff_age,"+"))) %>% 
    group_by(agegroup) %>% 
    summarise(ifr_mean=sum(mean*agegroupsize/sum(agegroupsize)),agegroup=unique(agegroup)[1],
          agegroup_mean=sum(agegroup_mean*agegroupsize/sum(agegroupsize)),agegroupsize=sum(agegroupsize)) %>%
    mutate(n=sum(agegroupsize),agegroup_perc=agegroupsize/n) %>% select(!n) %>% arrange(agegroup_mean) }

### get_OxCGRT -----------------------------

fcn_get_OxCGRT <- function(OxCGRT_url,cntr_name){
  OxCGRT_sel=subset(read_csv(OxCGRT_url),CountryName %in% cntr_name)
  # PLOT
  sel_cols=(unlist(lapply(OxCGRT_sel,class)) %in% "numeric") & !grepl("Confirmed|_|Legacy|ForDisplay",colnames(OxCGRT_sel))
  #   theme(legend.text=element_text(size=6),legend.position = "bottom")
  #
  OxCGRT_sel[,"date"]=as.Date(paste(sapply(strsplit(as.character(OxCGRT_sel$Date), "(?<=[0-9]{4})", perl=TRUE), "[[",1),
                                    sapply(strsplit(sapply(strsplit(as.character(OxCGRT_sel$Date), "(?<=[0-9]{4})", perl=TRUE), "[[",2), 
                                                    "(?<=[0-9]{2})", perl=TRUE),function(x){paste(x,collapse="-")}),sep="-"))
  OxCGRT_sel[,"NPI_on"]=0; OxCGRT_sel$NPI_on[min(which(OxCGRT_sel$StringencyIndex>0)):nrow(OxCGRT_sel)]=1
  # timespan of model
  OxCGRT_sel=merge(data.frame(date=seq(as.Date("2019-11-01"),max(OxCGRT_sel$date),1)),
                   OxCGRT_sel[,c("date","StringencyIndex","NPI_on")],by="date",all=TRUE)
  OxCGRT_sel$StringencyIndex[1:which.min(is.na(OxCGRT_sel$StringencyIndex))-1]=0
  OxCGRT_sel=OxCGRT_sel[!is.na(OxCGRT_sel$StringencyIndex),]; OxCGRT_sel$NPI_on[is.na(OxCGRT_sel$NPI_on)]=0
  # need to convert it into [0,1] to scale susceptibility (assume: pre-data period had not restrictions)
  OxCGRT_sel[,"OxCGRT_scaled"]=1-(OxCGRT_sel$StringencyIndex)/100 # timespan_dates=OxCGRT_sel$date
  OxCGRT_sel
}

### plot single simulation
fcn_plot_singlesim_separ_vars <- function(dfplot,NPI_input,popultot){
df_plot = dfplot %>% mutate(case_death_var_spec=ifelse((case_death_var %in% "case")|(cumul_trans_var %in% "cumul"),
          paste0(case_death_var," (%)"),paste0(case_death_var," (number)")),
          value_abs_fract_cond=ifelse((case_death_var %in% "case")|(cumul_trans_var %in% "cumul"),100*value/popultot,value))
ggplot(df_plot, aes(x=date,y=value_abs_fract_cond,group=name,color=name)) + geom_line(aes(linetype=state_delta_var),size=1.25) +
  scale_linetype_manual(values=c("twodash","solid")) + facet_wrap(cumul_trans_var~case_death_var_spec,scales="free",ncol=2) +
  theme_bw() + standard_theme + theme(legend.position="top") + 
  geom_vline(xintercept=as.Date(unique(df_plot$introd_date)),color="red") + 
  geom_vline(xintercept=as.Date(NPI_input$date[min(which(NPI_input$NPI_on>0))]),color="green") + 
  scale_x_date(date_breaks="2 weeks",expand=expansion(0.01,0)) + xlab('time (days)') + ylab(('')) + 
  labs(linetype="variables",color="variables") + labs(title=paste('introduction date:',unique(df_plot$introd_date), 
      ", initial import (E):",df_plot$seed_size,"/day for ",unique(df_plot$seeding_duration)," days"),
      subtitle=paste0("R0=",round(unique(df_plot$R0),2),", compliance=",df_plot$npi_compliance,
                                                                 ", IFR (sympt)=",round(1e2*df_plot$IFR_sympt,2),"%") )
}

### calc R0 ----------------------
# fcn_calc_R0 <- function(beta,ddeath,drecov,IFR){ beta*((1-IFR)/drecov + IFR/ddeath) }

### calc attack rate and cumul death rate ----------------------
fcn_attack_death_rate <- function(df,recov_var,plot_flag,y_breaks_lims){
  if (recov_var=="R"){
  df_output = subset(df,t %in% c(min(t),max(t)) & cumul_trans_var=="cumul")[,1:3] %>% 
    pivot_wider(names_from=c(name),values_from=value) %>% mutate(attack_rate=c(0,sum(c(R[t==max(t)],D[t==max(t)]))/S[t==min(t)]),
                                                                 cumul_death_rate=c(0,D[t==max(t)]/S[t==min(t)]))
  } else if (recov_var=="R_recov") {
    df_output = subset(df,t %in% c(min(t),max(t)) & cumul_trans_var=="cumul")[,1:3] %>% 
      pivot_wider(names_from=c(name),values_from=value) %>% mutate(attack_rate=c(0,sum(c(R_recov[t==max(t)],D[t==max(t)]))/S[t==min(t)]),
                                                  cumul_death_rate=c(0,D[t==max(t)]/S[t==min(t)]))}
  if (plot_flag=="plot"){
  ggplot(subset(df_output,t==max(t)) %>% select(t,attack_rate,cumul_death_rate) %>% 
    pivot_longer(cols=c(attack_rate,cumul_death_rate))) + geom_segment(aes(x=0,xend=0+1,y=value*1e2,yend=value*1e2,color=name),size=2) +
    theme_bw() + standard_theme + scale_x_continuous(breaks=NULL) + xlab("") + ylab("%") + labs(color="") +
    scale_y_log10(breaks=y_breaks_lims$breaks,limits=y_breaks_lims$limits) + 
      geom_text(aes(y=1.1*value*1e2,label=round(value*1e2,2)),x=0.5) }
  else if (grepl("calc",plot_flag)){df_output} else {print("choose <calc> or <plot>")}
}

### plot single simul vs data -----------------
fcn_plot_singlesim_data_error_deaths <- function(df,outbdrdailyestimates,selvars,NPI_input,baselineburialrate,
                                                death_rate_percap,n_per_persday,popul,datelims,plot_flag){
scale_factor=baselineburialrate/(death_rate_percap*popul/n_per_persday)
# rmse
singlesim_comp=left_join(subset(df,name %in% selvars)[,c("t","date","value")] %>% filter(t %% 1==0) %>% 
    rename(value_simul=value) %>% mutate(value_simul=value_simul) %>% mutate(value_simul_scaled=value_simul*scale_factor),
    outbdrdailyestimates[,c("date","rollmeanweek")] %>% rename(value_data=rollmeanweek) %>% mutate(value_data=value_data),by="date") %>% 
    mutate(sq_error=(value_data-value_simul_scaled)^2,mae=abs(value_data-value_simul_scaled)) %>% 
    filter(!is.na(value_data)) # rmse=sqrt(mean(singlesim_comp$sq_error))
  
  if (nchar(plot_flag)==0){rmse } else {
  ggplot(singlesim_comp[,c("date","value_simul_scaled","value_data")] %>% pivot_longer(col=!date),aes(x=date,y=value,color=name)) + 
    geom_line() + geom_point() + scale_x_date(date_breaks="1 week",expand=expansion(0.01,0),limits=as.Date(datelims)) + 
    scale_y_continuous(expand=expansion(0.01,0)) + theme_bw() + standard_theme + ggtitle(paste0("MAE=",round(mean(singlesim_comp$mae),2))) + 
      labs(caption=paste('introduction date:',unique(df$introd_date), ", initial import (E):",df$seed_size,"/day for ",
      unique(df$seeding_duration)," days"),subtitle=paste0("R0=",round(unique(df$R0),2),", compliance=",df$npi_compliance,
        ", IFR (all inf.)=",round(1e2*df$IFR_sympt,2),"%") ) + geom_vline(xintercept=as.Date(unique(df$introd_date)),color="red") + 
      geom_vline(xintercept=as.Date(NPI_input$date[min(which(NPI_input$NPI_on>0))]),color="green")
    }
}

### plot peak dates cases vs deaths of single simul --------------
fcn_plot_peakdates <- function(df_plot,datelims,datebreakval,sel_vars,yminval,logflag){
  peak_dates=subset(df_plot,name %in% sel_vars) %>% group_by(name) %>% summarise(maxval_date=date[value==max(value,na.rm=T)])
  peakplot <- ggplot(subset(df_plot, name %in% sel_vars) %>% filter(t %% 1 ==0),aes(x=date,y=value,group=name,color=name)) + geom_line() +
    geom_point() + scale_x_date(limits=as.Date(datelims),date_breaks=datebreakval,expand=expansion(0.01,0)) +
    geom_vline(data=peak_dates,aes(xintercept=maxval_date,color=name),size=1.5) + labs(color="peak of infections/deaths") + 
    theme_bw() + standard_theme + theme(axis.text.x=element_text(vjust=0.5)) 
   if (nchar(logflag)>0)  {
     peakplot <- peakplot+ scale_y_log10(limits=c(yminval,max(subset(df_plot,name %in% sel_vars)$value)),expand=expansion(0.01,0))}
  peakplot
}

### best paramsets  --------------
fcn_best_paramsets <- function(dferrors,paramtable,threshold,partype)  {
  if (grepl("numer",partype)){
  best_parsets=left_join(dferrors %>% rename(err_type=name,error=value) %>% group_by(CDR,err_type) %>% 
        filter(error<quantile(error,threshold)), paramtable %>% rowid_to_column(var="parset_ID") %>% 
        select(-c(seedsize_scanvals,introd_date_scanvals)),by="parset_ID") %>% 
  mutate(R0=betaval*((1-IFR)*d_c+IFR*d_death),IFR=1e2*IFR,compliance=1e2*compliance) %>% pivot_longer(cols=-c(parset_ID,err_type,error,CDR))
    list(best_parsets,best_parsets  %>% group_by(err_type,CDR,name) %>% 
           summarise(meanval=mean(value),medianval=median(value),minval=min(value),maxval=max(value)) )
    } else {
  # dates # best_paramsets_dates=
  best_parsets= left_join(error_pred_data_comp %>% rename(err_type=name,error=value) %>% group_by(CDR,err_type) %>% 
      filter(error<quantile(error,threshold)),paramtable %>% rowid_to_column(var="parset_ID") %>% 
      select(-c(seedsize_scanvals,betaval,IFR,compliance)),by="parset_ID") %>% pivot_longer(cols=-c(parset_ID,err_type,error,CDR)) 
  list(best_parsets,
  best_parsets %>% group_by(err_type,CDR,name) %>% 
    summarise(meanval=mean(value),medianval=median(value),minval=min(value),maxval=max(value)) %>% mutate(name="date of introduction"))
      }
}

### plot best trajs together with data
fcn_plot_best_fits <- function(dfplot,df_data,OxCGRTinput,seeddur,baselineburial,CDRval,CDR_metric,plotdatelims,popul,geomtext_vals,
                               ylimval,sizelimval,dc,ddeath,title_str,colorpal){
  scalefactor=baselineburial/(popul*CDRval/CDR_metric) # err_type
  min_error_vals=dfplot %>% group_by(introd_date_scanvals) %>% filter(value_error==min(value_error,na.rm=T)) %>% 
    summarise(compliance=unique(compliance),R0=unique(R0),min_error=min(value_error,na.rm=T),IFR=unique(IFR))

ggplot(subset(dfplot,date>=plotdatelims[1] & date<plotdatelims[2])) + 
  geom_line(aes(x=date,y=value*scalefactor,group=R0,color=R0,size=best_score)) + # ,alpha=best_score/max(best_score)
  geom_line(data=subset(df_data,date<=max(dfplot$date)),aes(x=date,y=new_graves_best_ipol-baselineburial),
            color="black",linetype="dashed",size=0.25) +
  geom_line(data=subset(df_data,date<=max(dfplot$date)),aes(x=date,y=rollmeanweek),color="black") +
  facet_wrap(~introd_date_scanvals,scales="free") + 
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=8)) + 
  geom_vline(aes(xintercept=introd_date_scanvals),color="red") +
  geom_vline(xintercept=as.Date(min(subset(OxCGRTinput,NPI_on>0)$date)),color="green") +
  scale_size_continuous(limits=c(min(df_plot$best_score),sizelimval)) + scale_color_gradientn(colours=colorpal) +
  scale_x_date(date_breaks="2 weeks",expand=expansion(0,0),limits=plotdatelims) +
  scale_y_continuous(expand=expansion(0.01,0),limits=c(0,ylimval)) + xlab('time (days)') + 
  ylab('deaths/day above baseline') + labs(title=title_str,subtitle=paste0("CDR=",CDRval,", seeding duration=",seeddur,
        ", error=",unique(dfplot$name_error)), color="R0") + guides(size=FALSE) + 
  geom_text(data=min_error_vals,aes(label=paste0("error (best)=",round(min_error,max(2-round(log10(min_error)),0)),", compl=",compliance,
                      ",\nR0=",R0,",IFR=",IFR*1e2,"%"),x=geomtext_vals[[1]],y=ylimval*geomtext_vals[[2]]),size=3,hjust=0)
}

### param scan parallel or sequential --------------
fcn_paramscan_parallel_seq <- function(paral_seq,k_lim,sir_varnames,varcateglist,timesteps,paramtable,transm_pars,clinfract,
                                       OxCGRT_input,N_tot,subpopul,filterdates,parnames_number){
  g(d_c,d_e,d_death,seed_dur) %=% transm_pars
  if (k_lim=="full"){k_lim=nrow(paramtable); print(paste0(nrow(paramtable)," parsets"))}
  
  if (grepl("par",paral_seq)){
    fcn_somal_sir_singlesimul<-.GlobalEnv$fcn_somal_sir_singlesimul
    df_output <- foreach(k=1:k_lim,.combine=rbind,.packages=c("tidyr","deSolve","dplyr","RcppRoll")) %dopar% {
      # RUN & process output
      temp<-fcn_somal_sir_singlesimul(sir_varnames,varcateglist,timesteps,paramtable$introd_date_scanvals[k],
                             c("beta"=paramtable$betaval[k],"gamma"=1/d_c,"sigma"=1/d_e,"d_death"=1/d_death,"IFR"=paramtable$IFR[k]), 
                             "clin_fract"=clinfract,"seed_size_val"=paramtable$seedsize_scanvals[k],"seeding_duration"=seed_dur,
                             OxCGRT_input,compliance=paramtable$compliance[k], subpopul*N_tot/sum(N_tot),c(),"")[[1]] %>% 
        filter(name %in% "new_deaths") %>% select(-any_of(c("case_death_var","cumul_trans_var","state_delta_var","introd_date",
        "seed_size","IFR","npi_compliance","beta","seeding_duration"))) %>% filter(date %in% filterdates) %>% mutate(parset_ID=k)
    }
  } else {
    df_output=list()
    for (k in 1:k_lim) { # nrow(paramtable)
      # RUN & process output
  df_ode_solution=fcn_somal_sir_singlesimul(sir_varnames,varcateglist,timesteps,day0=paramtable$introd_date_scanvals[k],
    num_params=c(beta=paramtable$betaval[k],gamma=1/d_c,sigma=1/d_e,d_death=1/d_death,IFR=paramtable$IFR[k]),
    clin_fract=clinfract,seed_size_val=paramtable$seedsize_scanvals[k],seeding_duration=transm_pars["seeding_duration"],
    OxCGRT_input,compliance=paramtable$compliance[k], N_tot*mogadishu_popul/sum(N_tot),c(),"")[[1]] %>% filter(name %in% "new_deaths") %>%
    select(-any_of(c("case_death_var","cumul_trans_var","state_delta_var","introd_date","seed_size","IFR","npi_compliance","beta",
                     "seeding_duration"))) %>% filter(date %in% filterdates) %>% mutate(parset_ID=k)
    df_output[[k]]=df_ode_solution 
      # print progress
      if (k %% 10 == 0) {print(paste0(k,",",round(k/k_lim,3)*1e2,"%, ",
                                      paste0(c("params:",as.character(paramtable[k,])),collapse=" ")))}
    }
    df_output=bind_rows(df_output, .id="column_label") 
    # %>% mutate(parset_ID=group_indices(df_output,.dots=parnames_number)) # datasource="model",
  }
  n_parcol=which(colnames(df_output) %in% "introd_date") 
  df_output # %>% mutate(parset_ID=group_indices(df_output,.dots=parnames_number))
}

### ### ### ### ### ### ### ### ### ### ### ###
### heatmap of MSE ------------------------
# PLOT
# err_type="negbinom"; n_text=2; n_round=0
# for (k in 1:length(scan_params$IFR)){
# df_heatmap=subset(error_pred_data_comp,IFR==scan_params$IFR[k] & grepl(err_type,name)) #  & beta==scan_params$betaval[2]
# colnames(df_heatmap)[grepl("npi",colnames(df_heatmap))]="npi_c"; colnames(df_heatmap)[grepl("durat",colnames(df_heatmap))]="seed_t"
# # best 1% param sets
# # minval_mse=min(subset(error_pred_data_comp,grepl(err_type,name))[,"value"]) 
# min_error_limit=10^(floor(log10(sort(subset(error_pred_data_comp,grepl(err_type,name))$value)[round(0.01*nrow(error_pred_data_comp))])*10)/10 )
# hmap_lims=summary(subset(error_pred_data_comp,grepl(err_type,name) & !is.infinite(value))$value)[c("Min.","Median","Max.")]*c(0.9,1,1.1)
# highl_tiles=subset(df_heatmap, value<min_error_limit) # %>% select(introd_date,seed_size) # midpoint_val=2*min(df_heatmap$value)
# # PLOT
# ggplot(df_heatmap, aes(x=seed_size,y=introd_date,fill=value)) + geom_tile(color="black") + 
#  geom_text(aes(label=round(value,n_round)),color="black",size=n_text) + geom_rect(data=highl_tiles,size=1,fill=NA,colour="green",
#   aes(xmin=seed_size-unique(diff(scan_params$seedsize_scanvals))/2,xmax=seed_size+unique(diff(scan_params$seedsize_scanvals))/2,
#   ymin=introd_date-unique(diff(scan_params$introd_date_scanvals))/2,ymax=introd_date+unique(diff(scan_params$introd_date_scanvals))/2)) +
#   theme_bw() + standard_theme + theme(axis.text.y=element_text(size=6)) + 
#   scale_fill_gradientn(colors=c("blue","white","red"),values=scales::rescale(array(hmap_lims))) + #,limits=c(0,hmap_lims[3]) 
#   facet_grid(npi_c~seed_t~beta,labeller=labeller(npi_c=label_both,seed_t=label_both,beta=label_both)) + 
#   scale_y_date(breaks=unique(df_heatmap$introd_date),expand=c(0,0)) + scale_x_continuous(breaks=unique(df_heatmap$seed_size),expand=c(0,0)) +
#   labs(title="Error for grid search",subtitle=paste0("IFR=",unique(df_heatmap$IFR)*100,"%, error: ",
#   gsub(gsub(err_type,"poiss","NLL, poisson"),"negbinom","NLL, neg.binom."))) + xlab("seed size") + ylab("introd. date")
# # SAVE
# filetag=paste0("IFR_",1e2*unique(df_heatmap$IFR),"_error_",err_type) 
# ggsave(paste0("simul_output/somalia_output/errors/error_parscan_SIR_heatmap_",filetag,".png"),width=32,height=40,units="cm")
# print(k)
# }

### object size -------------- 
fcn_objs_mem_use <- function(min_size){
  mem_use_df=round(data.frame(unlist(sapply(ls(envir=.GlobalEnv), function(n) object.size(get(n)), simplify = FALSE)))/1e6,1)
  colnames(mem_use_df)[1]<-"size (Mb)"; mem_use_df[,"objs"]=rownames(mem_use_df)
  mem_use_df<-mem_use_df[order(mem_use_df$size,decreasing=T),]; rownames(mem_use_df)<-c()
  mem_use_df[mem_use_df$size>min_size,c(2,1)]
}

### progress bar for parallel processing --------------
f_progress <- function(){
  pb <- txtProgressBar(min=1, max=n-1,style=3); count <- 0
  function(...) {
    count <<- count + length(list(...)) - 1
    setTxtProgressBar(pb,count);     Sys.sleep(0.01);     flush.console() ;    c(...)
  }
}

############
### assign multiple variables -----------------
'%=%' = function(l, r, ...) UseMethod('%=%')
# Binary Operator
'%=%.lbunch' = function(l, r, ...) {
  Envir = as.environment(-1)
  if (length(r) > length(l))
    warning("RHS has more args than LHS. Only first", length(l), "used.")
  if (length(l) > length(r))  {
    warning("LHS has more args than RHS. RHS will be repeated.")
    r <- extendToMatch(r, l)
  }
  for (II in 1:length(l)) {
    do.call('<-', list(l[[II]], r[[II]]), envir=Envir)
  }
}
extendToMatch <- function(source, destin) {
  s <- length(source)
  d <- length(destin)
  # Assume that destin is a length when it is a single number and source is not
  if(d==1 && s>1 && !is.null(as.numeric(destin)))
    d <- destin
  dif <- d - s
  if (dif > 0) {
    source <- rep(source, ceiling(d/s))[1:d]
  }
  return (source)
}
# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}
