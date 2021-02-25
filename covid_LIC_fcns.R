# from name and age group to linear index ----------------------------------------------------------
fun_sub2ind_seir_agestr=function(j_age,varname,varname_list,n_var,n_age){
  varnum=which(varname_list %in% varname); k=(j_age-1)*n_var + varnum; k }

# generate all model names ----------------------------------------------------------
fun_seir_agestr_varnames=function(varname_list,n_age){
  array( sapply(1:n_age, function(x_age) {sapply(varname_list, function(x) {paste0(x,'_',x_age)}) } ) )
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
fcn_covidm_df <- function(lmic_simul,sel_vars,params){   # c("S","E","Ip","Is","R")
agegr_pop=lmic_simul[lmic_simul$t==0 & lmic_simul$compartment %in% sel_vars] %>% group_by(population,group) %>% 
  summarise(pop=sum(value)); agegr_pop <- agegr_pop %>% mutate(pop_perc=pop/sum(pop))
colnames(lmic_simul)[colnames(lmic_simul) %in% "value"]="number"
lmic_simul=left_join(lmic_simul,agegr_pop,by=c("group","population")); 
lmic_simul[,"percent_agegr_popul"]=round(lmic_simul$number/lmic_simul$pop,4)
# fraction of the compartment at each timepoint
lmic_simul=lmic_simul %>% group_by(t,compartment,population) %>% mutate(percent_cases_at_t=number/sum(number))
# cumulative symptomatic attack rate. distribution of E between Is and Ia set by age-dependent parameter 
age_sympt_fract=data.frame(group=unique(lmic_simul$group),sympt_fract=params$pop[[1]]$y)
lmic_simul_cumulcases=lmic_simul[lmic_simul$compartment %in% "R",]
lmic_simul_cumulcases=left_join(lmic_simul_cumulcases,age_sympt_fract,by="group")
for (k in c("number","percent_agegr_popul","percent_cases_at_t")) {
  lmic_simul_cumulcases[,k]=lmic_simul_cumulcases[,k]*lmic_simul_cumulcases[,"sympt_fract"]}
lmic_simul_cumulcases$compartment="cumul_sympt_cases"
if (!any(grepl("cumul",colnames(lmic_simul)))) {
  lmic_simul=rbind(lmic_simul,lmic_simul_cumulcases[,!colnames(lmic_simul_cumulcases) %in% "sympt_fract"]) }
lmic_simul <- lmic_simul %>% pivot_longer(cols=c(number,percent_agegr_popul,percent_cases_at_t),names_to="variable")
lmic_simul
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
fun_lin_approx_agedep_par=function(min_val,max_val,rep_min,rep_max){
n_intersteps=nrow(age_groups)-(rep_max+rep_min)+1; slope=(max_val-min_val)/n_intersteps
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
# simple_sir_somalia <- function(t,x,parms) {
#   S<-x[1]; E<-x[2]; I_R<-x[3]; I_D=x[4]; R=x[5]; D=x[6]; newinfvar=x[7]
#   with(as.list(c(x,parms)), {
#     dS=-beta_dyn*S*(I_R+I_D); dE=beta_dyn*S*(I_R+I_D)-sigma*E;
#     dI_R=(1-IFR_estim)*sigma*E - gamma*I_R; dI_D=IFR_estim*sigma*E - d_death*I_D # infect who'll recover or die
#     dR=gamma*I_R; dD=d_death*I_D; d_newinfvar=sigma*E; d_beta=0; # if (t<100) {print(c(t,beta_dyn))}
#     return(list(c(dS,dE,dI_R,dI_D,dR,dD,d_newinfvar,d_beta))) }) }
# # event function with importation and stringency index
# eventfun <- function(t,y,parms){with (as.list(c(y,parms)),{
#   E <- E + seed_size*round(t>day0_num & t<(day0_num+seed_t)); beta_dyn <- beta*timevar_input[floor(t)+1]
#   return(c(S,E,I_R,I_D,R,D,newinfvar,beta_dyn)) }) }
# # by interpolation (defined within the fcn)
# simple_sir_somalia_interpol <- function(t,x,parms) {
#   S<-x[1]; E<-x[2]; I_R<-x[3]; I_D=x[4]; R=x[5]; D=x[6]; newinfvar=x[7]
#   with(as.list(c(x,parms)), {
#     beta_variable=beta*input_npi(t);
#     dS=-beta_variable*S*(I_R+I_D); dE=beta_variable*S*(I_R+I_D)-sigma*E + input_seeding(t)*4
#     dI_R=(1-IFR_estim)*sigma*E - gamma*I_R; dI_D=IFR_estim*sigma*E - d_death*I_D # infect who'll recover or die
#     dR=gamma*I_R; dD=d_death*I_D; d_newinfvar=beta_variable*S*(I_R+I_D); d_beta=0; # if (t<100) {print(c(t,beta_dyn))}
#     return(list(c(dS,dE,dI_R,dI_D,dR,dD,d_newinfvar,d_beta))) }) }

### function to run and plot single simul for somalia SIR model -------------
fcn_somal_sir_singlesimul <- function(sir_varnames,var_categ_list,timesteps,day0,num_params,clin_fract,
                                      seed_size_val,seeding_duration,OxCGRT_input,compliance,N_tot,xlimvals,plot_flag){
popul_tot=sum(N_tot); timespan_dates=OxCGRT_input$date
betaval=num_params['beta']; gamma=num_params['gamma']; d_death=num_params['d_death']; IFR_estim=num_params['IFR']
sigma=num_params['sigma']
initvals_S_I=c(popul_tot,rep(0,length(sir_varnames)-1),betaval/popul_tot); names(initvals_S_I)=c(sir_varnames,"beta_dyn")
day0_num=as.numeric(day0-timespan_dates[1]); suscept_reduct=1-OxCGRT_input$OxCGRT_scaled

timesteps<-seq(0,length(OxCGRT_input$OxCGRT_scaled),by=0.25)
timevar_signal<-data.frame(t=(1:nrow(OxCGRT_input))-1,npi_index=1-suscept_reduct*OxCGRT_input$NPI_on*compliance,seeding=0)
timevar_signal$seeding[timevar_signal$t %in% day0_num:(day0_num+seeding_duration)]=seed_size_val
# print(timevar_signal)
input_npi <- approxfun(timevar_signal[,c("t","npi_index")]) # ,method="constant"
input_seeding <- approxfun(timevar_signal[,c("t","seeding")]) # ,method="constant"
### create ODE object
params=list(beta=betaval/popul_tot,gamma=gamma,sigma=sigma,d_death=d_death,IFR_estim=IFR_estim,
            seed_size=seed_size_val,inv_timestep=1/unique(diff(timesteps)))
simple_sir_somalia_interpol <- function(t,x,parms) {
  S<-x[1]; E<-x[2]; I_R<-x[3]; I_D=x[4]; R=x[5]; D=x[6]; newinfvar=x[7]
  with(as.list(c(x,parms)), {
    beta_variable=beta*input_npi(t);
    dS=-beta_variable*S*(I_R+I_D); dE=beta_variable*S*(I_R+I_D)-sigma*E + input_seeding(t)*inv_timestep
    dI_R=(1-IFR_estim)*sigma*E - gamma*I_R; dI_D=IFR_estim*sigma*E - d_death*I_D # infect who'll recover or die
    dR=gamma*I_R; dD=d_death*I_D; d_newinfvar=beta_variable*S*(I_R+I_D); d_beta=0; # if (t<100) {print(c(t,beta_dyn))}
    return(list(c(dS,dE,dI_R,dI_D,dR,dD,d_newinfvar,d_beta))) }) }
####
# RUN
df_ode_solution <- euler(y=initvals_S_I,times=timesteps,func=simple_sir_somalia_interpol,parms=params) %>%
  # lsoda(y=initvals_S_I,times=timesteps,func=simple_sir_somalia,parms=params,events=list(func=eventfun,time=timesteps)) %>%
  as.data.frame() %>% setNames(var_categ_list$name_vars) %>% filter(t %% 1 ==0) %>% mutate(new_deaths=D-lag(D, default=D[1]),
  new_infections=newinfvar-lag(newinfvar,default=newinfvar[1]),symptom_cases=(I_R+I_D)*sum((N_tot/sum(N_tot))*clin_fract)) %>% 
  select(var_categ_list$sel_vars) %>% pivot_longer(cols=!t) %>% mutate(
  case_death_var=ifelse(name %in% var_categ_list$case_vars,"case","fatal"),
  cumul_trans_var=ifelse(name %in% var_categ_list$cumul_var,"cumul","transient"), 
 state_delta_var=ifelse(name %in% var_categ_list$delta_var,"delta","state"),
 date=timespan_dates[t+1],name=factor(name,levels=unique(name)),
 introd_date=day0,seed_size=seed_size_val,IFR=IFR_estim,
 npi_compliance=compliance, beta=betaval,seeding_duration=seeding_duration) %>% filter(!(is.na(date) | is.na(value))) #
#######
if (nchar(plot_flag)>0){
p<-ggplot(df_ode_solution, aes(x=date,y=value,group=name,color=name)) + geom_line(aes(linetype=state_delta_var),size=1.25) +
   scale_linetype_manual(values=c("twodash","solid")) + facet_wrap(cumul_trans_var~case_death_var,scales="free",ncol=2) +
   theme_bw() + standard_theme + theme(legend.position="top") + 
  geom_vline(xintercept=as.Date(day0),color="red") + 
  geom_vline(xintercept=as.Date(OxCGRT_input$date[min(which(OxCGRT_input$NPI_on>0))]),color="green") + 
   scale_x_date(date_breaks="2 weeks",limits=c(timespan_dates[xlimvals[1]],timespan_dates[xlimvals[2]])) + 
   xlab('time (days)') + ylab(('# cases')) + labs(linetype="variables",color="variables") +
    labs(title=paste('SIR model, introduction date:',unique(df_ode_solution$introd_date), ", initial import (E):",seed_size_val),
      subtitle=paste0("beta=",betaval,", compliance=",compliance) )} else {
     p<-c()
   }

list(df_ode_solution,p)
}

### param scan parallel or sequential --------------
fcn_paramscan_parallel_seq <- function(paral_seq,k_lim,sir_varnames,var_categ_list,timesteps,param_table,transm_pars,y_val,OxCGRT_input,
                                       N_tot,subpopul,newsurface_weekly){
  g(d_c,d_e,d_death) %=% transm_pars
  if (k_lim=="full"){k_lim=nrow(param_table); print(paste0(nrow(param_table)," parsets"))}
  
if (grepl("par",paral_seq)){
  fcn_somal_sir_singlesimul<-.GlobalEnv$fcn_somal_sir_singlesimul
  df_ode_solution_scan <- foreach(k=1:k_lim,.combine=rbind,.packages=c("tidyr","deSolve","dplyr","RcppRoll")) %dopar% {
  # RUN & process output
  temp<-fcn_somal_sir_singlesimul(sir_varnames,var_categ_list,timesteps,param_table$introd_date_scanvals[k],
  c("beta"=param_table$betaval[k],"gamma"=1/transm_pars[1],"sigma"=1/transm_pars[2],"d_death"=1/transm_pars[3],"IFR"=param_table$IFR[k]),
  clin_fract=y_val,seed_size_val=param_table$seedsize_scanvals[k],seeding_duration=param_table$seeding_duration[k],
     OxCGRT_input,compliance=param_table$compliance[k], N_tot*subpopul/sum(N_tot),c(),"")[[1]] %>% 
     filter(name %in% "new_deaths") %>% select(-any_of(c("case_death_var","cumul_trans_var","state_delta_var"))) %>%
     mutate(sum_weekly=roll_sum(value,7,fill=NA,align="right")) %>% filter(date %in% newsurface_weekly$date) } 
} else {
  df_ode_solution_scan=list()
  for (k in 1:k_lim) { # nrow(param_table)
  # RUN & process output
  df_ode_solution=fcn_somal_sir_singlesimul(sir_varnames,var_categ_list,timesteps,
            day0=param_table$introd_date_scanvals[k],
            num_params=c(beta=param_table$betaval[k],gamma=1/d_c,sigma=1/d_e,d_death=1/d_death,IFR=param_table$IFR[k]),clin_fract=y_val,
            seed_size_val=param_table$seedsize_scanvals[k], seeding_duration=param_table$seeding_duration[k],
            OxCGRT_input,compliance=param_table$compliance[k], N_tot*mogadishu_popul/sum(N_tot),c(),"")[[1]] %>% 
    filter(name %in% "new_deaths") %>% select(-any_of(c("case_death_var","cumul_trans_var","state_delta_var"))) %>%
    mutate(sum_weekly=roll_sum(value,7,fill=NA,align="right")) %>% filter(date %in% newsurface_weekly$date)
  df_ode_solution_scan[[k]]=df_ode_solution 
  # print progress
  if (k %% 10 == 0) {print(paste0(k,",",round(k/k_lim,3)*1e2,"%, ",
                           paste0(c("params:",as.character(param_table[k,])),collapse=" ")))}
  }
  df_ode_solution_scan=bind_rows(df_ode_solution_scan, .id = "column_label")  
}
  df_ode_solution_scan
}

### SEIR model with dynamic deaths ----------------------
seir_deaths_npi_interpol <- function(t,x,parms) {
  S<-x[1]; E<-x[2]; I_R<-x[3]; I_D=x[4]; D=x[5] # ; npi_val=x[6]
  gamma=0.2866; sigma=0.2503; d_death=1/21; inv_timestep=1/0.25; popul=2.2e6
  with(as.list(c(x,parms)), {
    beta_variable=(beta/popul)*(1-input_npi(t)*npi_compl)
    dS=-beta_variable*S*(I_R+I_D); dE=beta_variable*S*(I_R+I_D)-sigma*E + input_seeding(t)*seed_size*inv_timestep
    dI_R=(1-IFR_estim)*sigma*E - gamma*I_R; dI_D=IFR_estim*sigma*E-d_death*I_D; dD=d_death*I_D # infect who'll recover or die
    # print(c(t,1-input_npi(t)*npi_compl))
    return(list(c(dS,dE,dI_R,dI_D,dD))) }) } # dR=gamma*I_R


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
