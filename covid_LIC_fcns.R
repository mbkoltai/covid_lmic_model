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
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
    plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                     legend.title=element_text(size=14),legend.text=element_text(size=12),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))

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

########### parameter scan for one parameter & country ----------------------
fun_paramscan_singlepar=function(age_dep_paramtable,infect_first_ord_pars,infect_vartype,vartype_list,n_age,
                                 C_m, inf_initval,init_inf_age_groups,init_inf_vartype,N_tot,n_years,age_groups){
  for (age_dep_val in unique(age_dep_paramtable$par_scale_cnt)){
    parnames=unique(age_dep_paramtable$name)
  u_val=age_dep_paramtable$value[age_dep_paramtable$par_scale_cnt==age_dep_val & age_dep_paramtable$name %in% parnames[1]]
  y_val=age_dep_paramtable$value[age_dep_paramtable$par_scale_cnt==age_dep_val & age_dep_paramtable$name %in% parnames[2]]
  f_val=age_dep_paramtable$value[age_dep_paramtable$par_scale_cnt==age_dep_val & age_dep_paramtable$name %in% parnames[3]]
  # constant transmission parameters
  # d_e=1/4; d_p=1/1.5; d_c=1/3.5; d_s=1/5; infect_first_ord_pars=c(d_e,d_p,d_c,d_s)
  # KINETIC MATRIX (LINEAR TERMS)
  K_m=fun_seir_agestr_kinmatr(infect_vartype,vartype_list,n_age,infect_first_ord_pars,y_val)
  # vectors to scale FORCE of INFECTION terms
  list_bm_am=fun_force_inf_vects(vartype_list,forceinf_vars=c("S","E"),n_age,f_val,N_tot); b_m_full=list_bm_am[[1]]; a_m=list_bm_am[[2]]
  # inf vector is: b_m_full %*% diag(c(1,2)) %*% C_m %*% diag(u_val) %*% a_m %*% matrix(1,8,1)
  l=fun_inds_vartypes(n_age,vartype_list,infect_vartype); ind_all_inf_vars=l[[1]]; ind_all_susceptibles=l[[2]]
  # INITIAL CONDITIONS
  # seed epidemic by "E">0 in a given (or multiple) age groups
  initvals_seir_model=fcn_set_init_conds(inf_initval,init_inf_age_groups,init_inf_vartype,n_age,N_tot,vartype_list,ind_all_susceptibles)
  ### run ODEs ----------------------------
  # time resolution and duration
  # n_days_year=365; n_years=1; 
  max_time=n_years*365; timesteps <- seq(0,max_time,by=0.25)
  # define inputs for ODE fcn
  params=list(K_m, C_m, b_m_full, u_val, a_m, ind_all_susceptibles, ind_all_inf_vars)
  # RUN # 
  ode_solution<-lsoda(initvals_seir_model,timesteps,func=seir_agestr_ode,parms=params)
  print(paste("u (susc)=[",round(u_val[1]/(10^floor(log10(u_val[1]))),2),"e",floor(log10(u_val[1])), 
   round(u_val[length(u_val)]/(10^floor(log10(u_val[length(u_val)]))),2),"e",floor(log10(u_val[length(u_val)])),"] mean=",
   round(mean(u_val)/(10^floor(log10(mean(u_val)))),2),"e",floor(log10(mean(u_val))), "|",
              "y (clinfrac) =[",round(y_val[1],3),round(y_val[length(y_val)],3),"] mean=",round(mean(y_val),2),"|",
              "f (inf asym) =[",round(f_val[1],3),round(f_val[length(f_val)],3),"] mean=",round(mean(f_val),2), sep=" "))
  # PROCESS OUTPUT
  l_proc_sol=fcn_proc_ode_output(ode_solution,full_varname_list,ind_all_inf_vars); df_ode_solution=l_proc_sol[[1]]
  # infected fraction of popul (HIT)
  attack_rate=fun_attackrate_df_calc(df_ode_solution,N_tot,age_dep_val,y_val)
  # print(attack_rate)
  if (age_dep_val==0) {scan_table=attack_rate} else {scan_table=rbind(scan_table,attack_rate)}
  } # end of for loop
  scan_table[,"agegroup_names"]=factor(paste(age_groups$age_low,age_groups$age_high,sep = "-")[as.numeric(sapply(strsplit(as.character(
    scan_table$age_group),"_"),"[[",2))],levels=paste(age_groups$age_low,age_groups$age_high,sep = "-"))
  colorvar="col_scale"; scan_table[,colorvar]=factor(scan_table$par_scale)
  scan_table
}

### parameter scan for multiple parameters and countries --------------------------------
fun_multipar_multicntr_scan=function(cntr_list,year_str,age_groups,covidm_abs_path,scan_param_fullname,n_loop,
                                     midpoint_susc,delta_susc,midpoint_clinfract,delta_clinfr,mean_susc_exp,mean_clinfract,n_age,
                                     infect_first_ord_pars,infect_vartype,vartype_list,
                                     inf_initval,init_inf_age_groups,init_inf_var,n_years){
for (k_cntr in 1:length(cntr_list)){countryval=cntr_list[k_cntr]; N_tot=fun_cntr_agestr(countryval,i_year=year_str,age_groups) 
# uniform: N_tot=rep(1e6,n_age)
# CONTACT MATRIX
if (!exists("covid_params")){cm_force_rebuild=F; cm_build_verbose=T; cm_version=2; source(file.path(covidm_abs_path,"R","covidm.R"))}
covid_params=cm_parameters_SEI3R(gsub("Sudan","Ethiopia",countryval))
C_m=Reduce('+',covid_params$pop[[1]]$matrices)
numchars=paste0(as.character(round(c(sum(N_tot)/1e6,median(C_m),max(Re(eigen(C_m)$values))),1)),collapse = "," )
print(paste0(countryval," [popul,median(C_m),max(eigval(C_m))] = [",numchars,"]"))
# scan in parameter
# scan_param_fullname=c("susceptibility","clinical fraction","asymptomatic infectiousness")
for (k_par in 1:length(scan_param_fullname)){ # n_loop=6
age_dep_paramtable=fun_paramtable_susc_clinfract_inf(scan_param_fullname[k_par],n_loop,midpoint_susc,delta_susc,
                                                     midpoint_clinfract,delta_clinfr,mean_susc_exp,mean_clinfract,n_age)
attack_rates_scan=fun_paramscan_singlepar(age_dep_paramtable,infect_first_ord_pars,infect_vartype,vartype_list,n_age,C_m,inf_initval,
                                          init_inf_age_groups,init_inf_var,N_tot,n_years,age_groups)
attack_rates_scan[,"scanpar"]=scan_param_fullname[k_par]; if (k_par==1){attack_rates_scan_multipar=attack_rates_scan} else {
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

### calculate NGM and R0 -----------------
fun_NGM_R0=function(N_tot,y_val,u_val,f_val,C_m,lin_rates){
NGM=matrix(0,nrow=length(N_tot),ncol=length(N_tot)); g(d_p,d_c,d_s) %=% lin_rates
for (i_row in 1:length(y_val)){ for (j_col in 1:length(y_val)){
  # agegr_frac=N_tot[i_row]/sum(N_tot) # 1
  NGM[i_row,j_col]=u_val[i_row]*C_m[i_row,j_col]*(y_val[j_col]*(1/d_p+1/d_c) + (1-y_val[j_col])*unique(f_val)/d_s) } }
list(NGM,max(Re(eigen(NGM)$values)))
}

### multidimensional scan for R0 -----------------
fun_multidim_scan=function(N_tot,k_max,agedep_param_lims,lin_rates){
R0_scan=data.frame(matrix(0,ncol=5,nrow=k_max^3)); k_count=0; colnames(R0_scan)=c("susc","clinfract","asympt_inf","R0","HIT")
g(midpoint_susc,delta_susc,rep_min_susc,rep_max_susc,midpoint_clinfract,delta_clinfr,rep_min_clinfr,rep_max_clinfr) %=% agedep_param_lims
for (k_susc in 1:k_max){
  for (k_clinfrac in 1:k_max){
    for (k_infect in 1:k_max){
      u_val=fun_lin_approx_agedep_par(min_val=midpoint_susc-(k_susc/k_max)*delta_susc,
                                      max_val=midpoint_susc+(k_susc/k_max)*delta_susc,rep_min_susc,rep_max_susc)
      y_val=fun_lin_approx_agedep_par(midpoint_clinfract-(k_clinfrac/k_max)*delta_clinfr,
                                      midpoint_clinfract+(k_clinfrac/k_max)*delta_clinfr,rep_min_clinfr,rep_max_clinfr); 
      f_val=k_infect/k_max; k_count=k_count+1; 
      R0=fun_NGM_R0(N_tot,y_val,u_val,f_val,C_m,lin_rates)[[2]]; R0_scan[k_count,]=c(k_susc,k_clinfrac,k_infect,R0,1-1/R0)
    }
  } 
}
R0_scan[,"asympt_inf_str"]=factor(paste0("asympt inf=",R0_scan$asympt_inf/k_max),
                                  levels=unique(paste0("asympt inf=",R0_scan$asympt_inf/k_max)))
R0_scan
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
