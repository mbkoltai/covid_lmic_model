# libraries for population data and country codes
# load popul data: "pop" has totals, "popF", "popM" by age groups
data("pop")
# sum(popF[popF$name %in% "Sudan","2020"] + popM[popM$name %in% "Sudan","2020"]) ==  pop[pop$name %in% "Sudan","2020"]
sudan_pop_2020=pop[pop$name %in% "Sudan","2020"]*1e3
# sessionInfo()
# path
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# plotting theme
standard_theme=theme(panel.grid=element_line(linetype="dashed",colour="black",size=0.1),
                     plot.title=element_text(hjust=0.5,size=16),axis.text.x=element_text(size=9,angle=90),axis.text.y=element_text(size=9),
                     legend.title=element_text(size=14),legend.text=element_text(size=12),
                     axis.title=element_text(size=14), text=element_text(family="Calibri"))
### read in modelling scenario analysis --------------------
SDN_alls=qread('SDN/alls.qs')
# cumulative number of cases
cumul_number=SDN_alls$med[SDN_alls$scen_id==1 & SDN_alls$compartment %in% "R" & SDN_alls$age %in% "all"]
# plot age groups separately for one scenario
ggplot(SDN_alls[SDN_alls$scen_id==1,],aes(x=t,y=med,group=age,color=age)) + geom_line() + scale_color_discrete() +
  facet_wrap(~compartment,scales='free') + standard_theme + theme_bw() + xlab('days') + ylab('median value') + xlim(c(0,250))
# plot total # for all scenarios
ggplot(SDN_alls[SDN_alls$age %in% "all",],aes(x=t,y=med,group=scen_id,color=scen_id)) + geom_line() +
  scale_colour_gradientn(colours=terrain.colors(length(unique(SDN_alls$scen_id)))) + facet_wrap(~compartment,scales='free') + 
  standard_theme + theme_bw() + xlab('days') + ylab('median value') + xlim(c(0,6e2))
# TOTAL CASES by intervention scenario
truthvals_final_totalcases=SDN_alls$age %in% "all" & SDN_alls$t==max(SDN_alls$t) & SDN_alls$compartment %in% "R"
norm_factor=1 # 1, sudan_pop_2020
sudan_reported_deaths=total_reported_deaths_lastdate[total_reported_deaths_lastdate$country %in% "Sudan",]$value
if (norm_factor>1) {fract_tag='fraction'} else {fract_tag='number'}
ggplot(SDN_alls[truthvals_final_totalcases,],aes(x=scen_id,y=med/norm_factor,ymin=lo/norm_factor,ymax=hi/norm_factor)) + 
  geom_pointrange(fatten=1.8,size=0.8) + scale_x_continuous(breaks=1:length(unique(SDN_alls$scen_id))) + 
  # insert data
  theme_bw() + standard_theme +  xlab('intervention scenarios') + ylab(paste0('cumul cases ',fract_tag,' [median value, lo-hi]'))
# fatten=1.8,size=0.8,position=position_dodge(width=0.6)
# SAVE
cumulcase_filename=paste0('plots/cumul_case_',fract_tag,'_scen_id.png')
ggsave(cumulcase_filename,width=30,height=20,units="cm")
# TOTAL DEATHS
start_date="1-MAR-2020";days_since_beginning=as.numeric(as.Date(format(Sys.Date(),"%d-%B-%Y"),"%d-%B-%Y")-as.Date(start_date,"%d-%B-%Y"))
total_deaths_seir=SDN_alls[SDN_alls$compartment %in% "death_o" & SDN_alls$t<days_since_beginning,] %>% group_by(scen_id) %>% 
  summarise(sum_lolo=sum(lo.lo),sum_med=sum(med),sum_hihi=sum(hi.hi),sum_lo=sum(lo),sum_hi=sum(hi))
# geom_hline(yintercept=sudan_reported_deaths,linetype='dashed',size=rep_line_size,color='red')
ggplot(total_deaths_seir,aes(x=scen_id,y=sum_med/norm_factor,ymin=sum_lolo/norm_factor,ymax=sum_hihi/norm_factor)) + 
  geom_pointrange(fatten=3,size=0.6) + 
  scale_x_continuous(breaks=1:length(unique(SDN_alls$scen_id))) + scale_y_log10(breaks=10^c(0:6)) + coord_cartesian(ylim=c(1e2,5e5)) +
  geom_hline(yintercept=sudan_reported_deaths,linetype='dashed',size=rep_line_size,color='red') +
  theme_bw() + standard_theme + xlab('intervention scenarios') + ylab(paste0('cumul deaths ',fract_tag)) +
  ggtitle('initial SEIR predictions [2.5%, median, 97.5%] and reported deaths')
# labs(color='expected',linetype='reported')
####
cumuldeaths_filename=paste0('plots/cumul_deaths_',fract_tag,'_scen_id.png')
ggsave(cumuldeaths_filename,width=30,height=20,units="cm")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### compare expected death toll to reported/modelled --------------------
library(COVIDCurve) # remotes::install_github("mrc-ide/COVIDCurve")
IFR_by_age_imperial=read_csv("IFR_by_age_imperial.csv")
# deaths = prevalence*IFR_by_age*popul_by_age
sudan_popul_age_groups=popF[popF$name %in% "Sudan","2020"] + popM[popM$name %in% "Sudan","2020"]
unique(popF$age)
num_agegroups=unique(IFR_by_age_imperial$agegroup)[!sapply(as.numeric(substr(unique(IFR_by_age_imperial$agegroup),1,1)),is.na)]
# ifr_all_age_groups=IFR_by_age_imperial$mean[!grepl("without",IFR_by_age_imperial$variable) & 
#                                               IFR_by_age_imperial$agegroup %in% num_agegroups]
ifr_all_age_groups=IFR_by_age_imperial[!grepl("without",IFR_by_age_imperial$variable) & 
                                         IFR_by_age_imperial$agegroup %in% num_agegroups,c("mean","CI95_lower","CI95_upper")]
ifr_all_age_groups=rbind(ifr_all_age_groups,ifr_all_age_groups[rep(nrow(ifr_all_age_groups),
                                                                   length(unique(popF$age))-nrow(ifr_all_age_groups)),])
# convert from percentage to fraction
if (any(ifr_all_age_groups>1)) {ifr_all_age_groups=ifr_all_age_groups/1e2}

# calculate expected fatality for subsah african cntrs
subsah_afr_cntrs_iso3=c("AGO","BDI","BEN","BFA","BWA","CAF","CIV","CMR","COD","COG","COM","DJI","ERI","ETH","GAB","GHA","GIN",
                        "GMB","GNB","GNQ","KEN","LBR","LSO","MDG","MLI","MOZ","MRT","MWI","NAM","NER","NGA","RWA","SDN","SEN","SLE","SOM","SSD",
                        "STP","SWZ","TCD","TGO","TZA","UGA","ZAF","ZMB","ZWE") # "SYC","MUS"
subsah_afr_cntrs_name=countrycode(subsah_afr_cntrs_iso3,origin='iso3c',destination='un.name.en')
subsah_afr_cntrs_iso3n=countrycode(subsah_afr_cntrs_iso3,origin='iso3c',destination='iso3n')
# load data from ourworldindata
# https://covid.ourworldindata.org/data/ecdc/total_cases.csv
# https://covid.ourworldindata.org/data/ecdc/total_reported_deaths.csv
total_reported_deaths <- read_csv("https://covid.ourworldindata.org/data/ecdc/total_deaths.csv")
# few cntr names are different
subsah_afr_cntrs_name[!subsah_afr_cntrs_name %in% colnames(total_reported_deaths)] 
subsah_afr_cntrs_name[!subsah_afr_cntrs_name %in% colnames(total_reported_deaths)]=
  colnames(total_reported_deaths)[unlist(lapply(c("Ivoire",glob2rx("Democratic*Congo"),"Gambia","Bissau","Swaziland","Tanzania"), 
                                                function(x){which(grepl(x,colnames(total_reported_deaths)))} ))]
if (sum(!subsah_afr_cntrs_name %in% colnames(total_reported_deaths))==0){print("all cntrs found")} else{
  print(subsah_afr_cntrs_name[!subsah_afr_cntrs_name %in% colnames(total_reported_deaths)]) }

# age struct of selected cntrs
agestruct_subsah_afr_cntrs=data.frame(sapply(subsah_afr_cntrs_iso3n, function(x){(popF[popF$country_code %in% x,"2020"] + 
                                                                                    popM[popM$country_code %in% x,"2020"]) }))
colnames(agestruct_subsah_afr_cntrs)=subsah_afr_cntrs_name; rownames(agestruct_subsah_afr_cntrs)=unique(popF$age)
# total expected deaths for 100% prevalence = IFR_by_agegroup*popul_by_agegroup
subsah_afr_cntrs_total_exp_death=data.frame(country=subsah_afr_cntrs_name,t(sapply(1:ncol(agestruct_subsah_afr_cntrs),
                                                                                   function(x) {round(sapply(1:3, function(n){1e3*ifr_all_age_groups[,n]%*%agestruct_subsah_afr_cntrs[,x]}))})))
colnames(subsah_afr_cntrs_total_exp_death)[2:ncol(subsah_afr_cntrs_total_exp_death)]=colnames(ifr_all_age_groups)
# IFR for a given country
pop_struct=popF[popF$name %in% "Colombia","2020"] + popM[popM$name %in% "Colombia","2020"]
# agestruct_subsah_afr_cntrs$Sudan
1e2*(ifr_all_age_groups$mean %*% pop_struct/sum(pop_struct) )
# 

prevalence_range=c((1:19)/10,seq(2,80,0.5))/100
total_exp_death_preval_varied=do.call(rbind,lapply(prevalence_range, 
                                                   function(x) {cbind(country=subsah_afr_cntrs_total_exp_death$country,prevalence=x,x*subsah_afr_cntrs_total_exp_death[,2:4])}) )
# total reported deaths: total_reported_deaths
total_reported_deaths_lastdate=data.frame(
  country=colnames(total_reported_deaths)[colnames(total_reported_deaths) %in% subsah_afr_cntrs_name],
  value=t(total_reported_deaths[nrow(total_reported_deaths),
                                colnames(total_reported_deaths) %in% subsah_afr_cntrs_name])); row.names(total_reported_deaths_lastdate)=c()
# have same order as predicted values
total_reported_deaths_lastdate=total_reported_deaths_lastdate[match(unique(total_exp_death_preval_varied$country),
                                                                    total_reported_deaths_lastdate$country),]
# find prevalence level closest to reported deaths
total_exp_death_preval_varied[,'reported']=rep(total_reported_deaths_lastdate$value,length(prevalence_range))
total_exp_death_preval_varied[,'diff_from_mean']=total_exp_death_preval_varied$mean - total_exp_death_preval_varied$reported
best_approx_preval=total_exp_death_preval_varied %>% group_by(country) %>% slice(which.min(abs(diff_from_mean)))
# PLOT: expected deaths ~ (prevalence*ascertainment)
truthvals=total_exp_death_preval_varied$prevalence<(10^1.5)/100
# total_exp_death_preval_varied$country %in% as.character(unique(total_exp_death_preval_varied$country)[1:5])
rep_line_size=0.5; insert_txt_size=3
ggplot(total_exp_death_preval_varied[truthvals,]) + geom_line(aes(x=prevalence*1e2,y=mean,linetype='expected (IFR)'),color='blue') +
  geom_ribbon(aes(x=prevalence*1e2,ymin=CI95_lower,ymax=CI95_upper),fill='blue',alpha=0.1) + 
  # horizontal line showing reported deaths
  geom_hline(data=total_reported_deaths_lastdate,aes(yintercept=value,linetype='reported'),size=rep_line_size,color='red') + 
  # prevalence value providing best approx
  geom_vline(data=best_approx_preval,aes(xintercept=prevalence*1e2,linetype='reported'),size=rep_line_size,color='red') +
  geom_text(data=best_approx_preval,aes(x=10^(log10(prevalence)+0.3)*1e2,y=mean/3.5,label=paste0(prevalence*1e2,"%")),
            size=insert_txt_size) + 
  facet_wrap(~country,scales='free') + theme_bw() + standard_theme + 
  theme(axis.text.x=element_text(angle=0,size=7,vjust=1,hjust=0.5),axis.text.y=element_text(size=6,angle=90,hjust=0.5,vjust=0)) +
  scale_linetype_manual(name='',values=c(2,1),guide=guide_legend(override.aes=list(color=c("blue","red"),size=c(0.6,rep_line_size)))) +
  scale_x_log10() + scale_y_log10(breaks=scales::trans_breaks("log10",function(x) 10^x),
                                  labels=scales::trans_format("log10",scales::math_format(10^.x))) +
  xlab('% [prevalence x ascertainment]') + ylab('deaths') + labs(color='expected',linetype='reported') + 
  ggtitle('Expected prevalence using age-dependent IFRs IF no under-reporting')
# labs(color='data source',linetype='mean') + guides(xintercept=FALSE,linetype=guide_legend(ncol=2))
ggsave('plots/cumul_deaths_expected_observed.png',width=36,height=24,units="cm")

# make a plot where mean is the prevalence level where predicted deaths closest to reported by mean IFR,
# lower CI95 is where it's closest to lower CI95 of IFR, upper ~ ~ ~
best_approx_preval_ci95_low=total_exp_death_preval_varied %>% group_by(country) %>% slice(which.min(abs(CI95_lower-reported))) %>%
  mutate(diff_from_CI_low=CI95_lower-reported,prevalence_CI_low=prevalence)
best_approx_preval_ci95_up=total_exp_death_preval_varied %>% group_by(country) %>% slice(which.min(abs(CI95_upper-reported))) %>%
  mutate(diff_from_CI_up=CI95_upper-reported,prevalence_CI_up=prevalence)
best_approx_preval_all=bind_cols(best_approx_preval[,c("country","prevalence","mean","reported","diff_from_mean")],
                                 best_approx_preval_ci95_low[,c("prevalence_CI_low","CI95_lower")],best_approx_preval_ci95_up[,c("prevalence_CI_up","CI95_upper")])
best_approx_preval_all[,"country_no"]=1:nrow(best_approx_preval_all)
####
# pointrange plot
ggplot(best_approx_preval_all,aes(x=country,y=prevalence*1e2,ymin=prevalence_CI_up*1e2,ymax=prevalence_CI_low*1e2)) +
  geom_pointrange(fatten=2,size=0.6,linetype='dashed',color='blue') + # aes(,shape=,),position=position_dodge(width=0.6)
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=12,angle=90,hjust=0.95,vjust=0.2)) + 
  scale_y_log10(breaks=scales::trans_breaks("log10",function(x) 10^x),labels=scales::trans_format("log10",scales::math_format(10^.x))) + 
  geom_text(aes(x=country_no-0.5,y=prevalence*1e2,label=paste0(prevalence*1e2,"%")),size=4,color='red',angle=90) + 
  geom_text(aes(x=country_no,y=0.9*prevalence_CI_up*1e2,label=paste0(prevalence_CI_up*1e2)),size=3) +
  geom_text(aes(x=country_no,y=1.1*prevalence_CI_low*1e2,label=paste0(prevalence_CI_low*1e2)),size=3) + 
  xlab('') + ylab('% [prevalence x reporting]') + ggtitle('Best fit prevalence to reported death data (mean IFR +/- CI95)') +
  scale_x_discrete(expand=c(0.03,0))
# labs(color='data source',shape='data source',linetype='data source') + # guides(linetype=FALSE) +
ggsave('plots/best_fit_prevalence_range.png',width=36,height=20,units="cm")  

#### expected IFRs ---------------
# x="South Africa";pop_struct=popF[popF$name %in% x,"2020"]+popM[popM$name %in% x,"2020"]; 
# expected IFRs
cntr_list=c("Brazil") # "South Africa","Egypt","Sudan","Nigeria")
sapply(cntr_list,function(x) { 1e2*(ifr_all_age_groups$mean %*% (popF[popF$name %in% x,"2020"]+popM[popM$name %in% x,"2020"])/
                                      sum(popF[popF$name %in% x,"2020"]+popM[popM$name %in% x,"2020"]))})
# # # expected deaths with x% prevalence
prev=0.5; cntr_list="Sudan"
prev*sapply(cntr_list,function(x) {1e3*(ifr_all_age_groups$mean %*% (popF[popF$name %in% x,"2020"]+popM[popM$name %in% x,"2020"]) )})
# # reported deaths
# total_reported_deaths_lastdate[total_reported_deaths_lastdate$country %in% c("South Africa","Egypt","Sudan","Nigeria"),]
