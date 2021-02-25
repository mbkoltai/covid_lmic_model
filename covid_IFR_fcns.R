library(COVIDCurve); library(wpp2019) # remotes::install_github("mrc-ide/COVIDCurve")
IFR_by_age_imperial=read_csv("IFR_by_age_imperial.csv"); data("popF"); data("popM")
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
total_reported_deaths <- read_csv("https://github.com/owid/covid-19-data/raw/master/public/data/jhu/total_deaths.csv") 
# https://covid.ourworldindata.org/data/ecdc/total_deaths.csv
# some cntr names are different
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
# cntr_ex="Somalia"; pop_struct=popF[popF$name %in% cntr_ex,"2020"] + popM[popM$name %in% cntr_ex,"2020"]
# agestruct_subsah_afr_cntrs$Sudan

### transm parameters age dependence -------------------
age_groups <- data.frame(age_group=c(1:16), age_low=c(seq(0,75,5) ), age_high=c(seq(4,74,5),100))
# population data from wpp2019
# countryval="Sudan"; N_tot=fun_cntr_agestr(countryval,i_year="2020",age_groups) # uniform: N_tot=rep(1e6,n_age)
min_val_susc=1e-2; maxval_susc=0.15; midpoint_susc=(min_val_susc+maxval_susc)/2; delta_susc=maxval_susc-midpoint_susc
midpoint_clinfract=0.4; delta_clinfr=0.35; rep_min_susc=3; rep_max_susc=5; rep_min_clinfr=5; rep_max_clinfr=3; 
depval=1
u_val=fun_lin_approx_agedep_par(min_val=midpoint_susc-depval*delta_susc,max_val=midpoint_susc+depval*delta_susc,
                                rep_min_susc,rep_max_susc)
y_val=fun_lin_approx_agedep_par(midpoint_clinfract-depval*delta_clinfr,midpoint_clinfract+depval*delta_clinfr,
                                rep_min_clinfr,rep_max_clinfr)
clinical_fraction=y_val; rm(y_val)