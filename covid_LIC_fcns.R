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
fun_force_inf_vects=function(vartype_list,forceinf_vars,n_age,f_val){
b_m_full=matrix(0,length(vartype_list)*n_age,n_age)
for (i_age in 1:n_age) {for (k_var in 1:length(forceinf_vars)){
  b_m_full[(i_age-1)*length(vartype_list) + 1,i_age]=-1;b_m_full[(i_age-1)*length(vartype_list) + 2,i_age]=1
}}
# for summing infectious vars
a_m=matrix(0,n_age,length(infect_vartype)*n_age)
for (k in 1:n_age){  col_lims=((k-1)*length(infect_vartype)+1):(k*length(infect_vartype)); a_m[k,col_lims]=c(0,1,1,f_val[k]) }
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
df_ode_solution_tidy=reshape2::melt(df_ode_solution[,colSums(df_ode_solution)>0],id.vars='t')
df_ode_solution_tidy[,"vartype"]='noninf'; df_ode_solution_tidy[!grepl("S_",df_ode_solution_tidy$variable),"vartype"]='inf'  
df_ode_solution_tidy[,"compartm"]=sapply(strsplit(as.character(df_ode_solution_tidy$variable),"_"),"[[",1)
df_ode_solution_tidy[,"agegroup"]=sapply(strsplit(as.character(df_ode_solution_tidy$variable),"_"),"[[",2)
list(df_ode_solution,df_ode_solution_tidy)
}

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