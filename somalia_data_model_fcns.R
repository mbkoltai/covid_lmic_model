### set up plotting theme ------------
standard_theme=theme(# panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
  plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
  legend.title=element_text(size=14),legend.text=element_text(size=12),
  axis.title=element_text(size=14), text=element_text(family="Calibri"))
### colors ----
# library(wesanderson) # install.packages("wesanderson")
# pal <- wes_palette("Zissou1", 100, type = "continuous")

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

### merge IFRs above a given age (typically 75) ---------------------------
fcn_merge_ifr_above_age <- function(age_ifr_data,cutoff_age){
  age_ifr_data %>% mutate(agegroup=ifelse(agegroup_min<cutoff_age,as.character(agegroup),paste0(cutoff_age,"+"))) %>% 
    group_by(agegroup) %>% summarise(ifr_mean=sum(mean*agegroupsize/sum(agegroupsize)),agegroup=unique(agegroup)[1],
              agegroup_mean=sum(agegroup_mean*agegroupsize/sum(agegroupsize)),agegroupsize=sum(agegroupsize)) %>%
    mutate(n=sum(agegroupsize),agegroup_perc=agegroupsize/n) %>% select(!n) %>% arrange(agegroup_mean) }

### load cntr age struct, merge age groups above a cutoff ---------------------------
fcn_load_age_str <- function(cntrval,n_year,cutoff_age){
  data.frame(agegroup=popF[popF$name %in% cntrval,"age"],
             value=(popF[popF$name %in% cntrval,n_year] + popM[popM$name %in% cntrval,n_year])*1e3) %>%  
    mutate(agegroup_min=ifelse(grepl("-|\\+",as.character(agegroup)),
                               as.numeric(sapply(strsplit(as.character(agegroup),"-|\\+"),function(x) {x[1]})),NA),
           agegroup_max=ifelse(grepl("-|\\+",as.character(agegroup)),
                               as.numeric(sapply(strsplit(as.character(agegroup),"-|\\+"),function(x) {x[2]})),NA),
           agegroup=ifelse(agegroup_min<cutoff_age,as.character(agegroup),paste0(cutoff_age,"+"))) %>% group_by(agegroup) %>% 
    summarise(agegroup=unique(agegroup),value=sum(value),agegroup_min=unique(agegroup_min)[1],agegroup_max=max(agegroup_max)) %>%
    mutate(agegroup_max=ifelse(is.na(agegroup_max),agegroup_min+4,agegroup_max),agegroup_mean=(agegroup_min+agegroup_max+1)/2) %>%
    rename(agegroupsize=value) %>% arrange(agegroup_min)
}

### load IFRs from imperial estimates
fcn_load_ifr <- function(ifr_filepath){
  read_csv(ifr_filepath) %>% group_by(agegroup) %>% 
    summarise(agegroup=unique(agegroup),mean=mean(mean)/1e2,CI95_lower=mean(CI95_lower)/1e2,CI95_upper=mean(CI95_upper)/1e2) %>% 
    mutate(agegroup_min=ifelse(grepl("-|\\+",agegroup),as.numeric(sapply(strsplit(agegroup,"-|\\+"),function(x) {x[1]})),NA ),
           agegroup_max=ifelse(grepl("-",agegroup),as.numeric(sapply(strsplit(agegroup,"-"),function(x) {x[2]})),NA ) ) %>% 
    filter(!is.na(agegroup_min))
}

### approximate clinical fraction by "S" shape ------
# min_val=0.04; max_val=0.65; rep_min=5; rep_max=3; 
fun_lin_approx_agedep_par=function(agegroups,min_val,max_val,rep_min,rep_max){
  n_intersteps=nrow(agegroups)-(rep_max+rep_min)+1; slope=(max_val-min_val)/n_intersteps
  c(rep(min_val,rep_min),min_val+(1:(n_intersteps-1))*slope,rep(max_val,rep_max))
  # lines(age_groups$mid,clin_fract_approx)
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
fcn_covidm_singlesim_error <- function(covidmsimul,intr_date,seedsize,df_data,fitting_dates){
  df_data_simul_deaths=left_join(subset(covidmsimul %>% rename(simul=value),compartment %in% "death_o")%>% select(date,simul),
                                 df_data %>% select(date,daily_baseline_subtr) %>% rename(data=daily_baseline_subtr)) %>% pivot_longer(!date)
  abs_err=subset(df_data_simul_deaths, date >= fitting_dates[1] & date <= fitting_dates[2]) %>% group_by(date) %>%
    summarise(data=value[name=="data"],diff=value[name=="data"]-value[name=="simul"],sumval=sum(value),mae=abs(diff))
  ggplot(df_data_simul_deaths,aes(x=date)) + geom_line(aes(y=value,group=name,color=name)) + theme_bw() + standard_theme + 
    ggtitle(paste0("MAE=",round(mean(abs_err$mae),2), ", introd date: ", intr_date, ", seedsize=", sum(seedsize))) + 
    scale_x_date(limits = fitting_dates,date_breaks="2 weeks",expand=expansion(0.01,0)) + scale_y_continuous(expand=expansion(0.01,0))
}

### function to have COVIDM incidence outputs ---------------------------------
cm_multinom_process <- function(src, outcomes, delays, report="") {
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


### ### ### ### 
cm_plot_posterior_mod <- function(fit,plot_params)
{  if (!any("cm.fit" %in% class(fit))) {
  stop("cm_plot_posterior requires a cm.fit object.")  }
  
  # get prior distributions
  priors = NULL;
  for (pr in 1:length(fit$priors)) {
    ep = cm_evaluate_distribution(fit$priors[[pr]]);
    setDT(ep);
    ep[, p := p / max(p)]
    ep[, parameter := names(fit$priors)[pr]];
    priors = rbind(priors, ep);
  }
  posterior = melt(fit$posterior, id.vars = NULL, measure.vars = 5:ncol(fit$posterior), variable.name = "parameter");
  
  # put in correct order
  priors[, parameter := factor(parameter, levels = names(fit$priors))]
  posterior[, parameter := factor(parameter, levels = names(fit$priors))]
  
  # plot
  ggplot() + geom_line(data=priors, aes(x=x, y=p), linetype="22", size = plot_params$line_width) + # prior
    geom_histogram(data = posterior, aes(value, y=stat(ndensity), fill=parameter), alpha=0.75,bins=plot_params$n_bin) +
    facet_wrap(~parameter, scales = "free") +     labs(x = NULL, y = NULL) + ggtitle(paste0("CDR=",round(plot_params$cdr_val,3))) +
    theme(legend.position = "none", strip.background = element_blank()) + theme_bw() + standard_theme
}

### cm plot pairwise -----------------
cm_plot_pairwise_mod <- function(fit, plot_theme)
{  if (!any("cm.fit" %in% class(fit))) {     stop("cm_plot_pairwise requires a cm.fit object.");   }
  
  # get prior distributions
  priors = NULL;
  for (pr in 1:length(fit$priors)) {
    ep = cm_evaluate_distribution(fit$priors[[pr]]);
    setDT(ep);
    ep[, p := p / max(p)]
    ep[, parameter := names(fit$priors)[pr]];
    priors = rbind(priors, ep);
  }
  posterior = melt(fit$posterior, id.vars = NULL, measure.vars = 5:ncol(fit$posterior), variable.name = "parameter");
  
  # put in correct order
  priors[, parameter := factor(parameter, levels = names(fit$priors))];
  posterior[, parameter := factor(parameter, levels = names(fit$priors))];
  
  # assemble plots
  plotlist = list();
  i = 1;
  for (r in 1:length(fit$priors)) {
    for (c in 1:r) {
      rp = names(fit$priors)[r];
      cp = names(fit$priors)[c];
      if (r == c) { # posterior element on its own
        plot = ggplot() +
          geom_line(data = priors[parameter == rp], aes(x = x, y = p), linetype = "22") + # colour = "grey", , size = 0.25
          geom_histogram(data = posterior[parameter == rp], aes(value, y = stat(ndensity)), fill = "#88bbff", bins = 30) +
          labs(x = NULL, y = NULL, title = rp) + theme_bw() + plot_theme
      } else { # pairwise plot
        data = cbind(posterior[parameter == cp, .(x = value)], posterior[parameter == rp, .(y = value)]);
        plot = ggplot(data) + geom_density2d(aes(x, y), colour = "#ffbb88") +   labs(x = NULL, y = NULL) + theme_bw() + plot_theme}
      plotlist[[i]] = plot;
      i = i + 1;
    }
    for (c in seq_len(length(fit$priors) - r) + r) {
      plotlist[[i]] = NULL;
      i = i + 1;     }
    cat("\n");
  }
  plot_grid(plotlist = plotlist, nrow = length(fit$priors), ncol = length(fit$priors))
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
  if (dif > 0) {    source <- rep(source, ceiling(d/s))[1:d]  }
  return (source)
}
# Grouping the left hand side
g = function(...) {
  List = as.list(substitute(list(...)))[-1L]
  class(List) = 'lbunch'
  return(List)
}