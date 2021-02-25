rm(list=ls()); currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
lapply(c("tidyverse","deSolve","gtools","rstudioapi","wpp2019","countrycode","coronavirus",
         "RcppRoll","scales","dttr2","wpp2019","foreach","parallel","doParallel"), library,character.only=TRUE)

### Linelist data ----------------
#' ## LINELIST DATA
#' ### Testing data on national level
national_charts_testing=read_csv("somalia_linelist/national_charts_testing.csv")
selcolnames=c("Total Tested","7 day moving average","Number of test for 1 confirmed case","TOTAL TESTED (%)",
              "Confirmed cases","TESTED POSITIVE from total (%)","CUMULATIVE CASES","Death","CUMULATIVE DEATHS","Recovery",
              "CUMULATIVE RECOVERY","Positivity rate (%)","CFR (%)")
national_charts_testing = national_charts_testing %>% pivot_longer(cols=all_of(selcolnames))
# plot
ggplot(subset(national_charts_testing,!Month %in% 'TOTAL' & 
                !grepl("Recovery|CUMUL|TOTAL TESTED|from total",name)),aes(x=`Week Number`,y=value,group=1)) + 
  geom_line() + geom_point() + facet_wrap(~name,scales="free",nrow=2) + theme_bw() + standard_theme + xlab("") +
  geom_rect(aes(xmin=4,xmax=6,ymin=0,ymax=Inf),fill="pink",color=NA,alpha=0.01) + 
  labs(title="National-level COVID19 metrics, weeks 12-32",caption="source: line list")
#' Testing, detected cases and deaths went up together in weeks 15-16
# save
# ggsave("simul_output/somalia_output/national_charts_testing.png",width=32,height=18,units="cm")

#' ### By state
state_charts=read_csv("somalia_linelist/state_charts.csv")
state_charts = state_charts %>% pivot_longer(cols = !c("Month","Week Number","State"))
state_charts$`Week Number`=gsub("Week-","w",state_charts$`Week Number`)
state_charts[,"week"]=as.numeric(gsub("w","",state_charts$`Week Number`))
# variables are
# c("TOTAL TESTED", "Weekly Average", "Confirmed Case", "CUMULATIVE CASES", "Death", "CUMULATIVE DEATHES",
# "CFR (%)", "Recovery", "CUMULATIVE RECOVERY", "Positivity (%)")

# plot
# ggplot(subset(state_charts, !grepl("Recovery|CUMUL|from total|Weekly",name) & week<33),
#        aes(x=`Week Number`,y=value,group=1)) + geom_line() + geom_point(size=0.5,fill=NA,shape=1) +
#   facet_grid(name~State,scales = "free",switch="y") + theme_bw() + standard_theme + theme(axis.text.x = element_text(size=6))
# # save
# ggsave("simul_output/somalia_output/state_charts.png",width=32,height=18,units="cm")

# plot with facet_wrap
ggplot(subset(state_charts, !grepl("Recovery|CUMUL|from total|Weekly",name) & week<33),
       aes(x=`Week Number`,y=value,group=1)) + geom_line() + geom_point(size=0.5,fill=NA,shape=1) +
  facet_wrap(name~State,scales="free",ncol=length(unique(state_charts$State))) +
  theme_bw() + standard_theme + theme(axis.text.x=element_text(size=5),axis.text.y=element_text(size=7)) + 
  labs(title="COVID19 metrics by state, weeks 12-32",caption = "source: line list")
#' most testing in Banadir and Somaliland
# save
# ggsave("simul_output/somalia_output/state_charts_facetwrap.png",width=30,height=20,units="cm")

#' ### Attack rate by district
attackrate_by_district=read_csv("somalia_linelist/attackrate_by_district.csv")
selcolnames=c("Total Pop","Positivity","Attack rate per 100K","Number of Cases","Total Tested",
              "Increase or decrease from mean","tests_per_100k","tests_per_case")
attackrate_by_district=attackrate_by_district %>% 
  mutate(tests_per_100k=(`Total Tested`/`Total Pop`)*1e5,tests_per_case=`Total Tested`/`Number of Cases`) %>% 
  pivot_longer(cols=all_of(selcolnames)) %>% 
  mutate(jointname=ifelse(name %in% c("Attack rate per 100K pop","tests_per_100k"),"attackrate_tests_per100k",name))
# plot
p <- ggplot(subset(attackrate_by_district,!grepl("Increase|Total Pop|Number of|Total Tested|tests_per_case",name)),
            aes(x=`DISRICT CODE`,y=value,group=1)) + geom_line() + geom_point(fill=NA,shape=1) + 
  facet_wrap(name~STATE,scales="free",ncol=7) + ylab("") + # facet_grid(name~STATE,scales="free",switch="y") + 
  theme_bw() + standard_theme + theme(axis.text.x = element_text(size=7),legend.position="top",legend.title=element_blank()) +
  labs(title="Cumulative attack rate, positivity and testing rate by state",caption = "source: line list")
p
#' Most data is from Banadir, some from Somaliland. Only in Banadir have most districts a non-zero/NA value.
#
# if (class(p$facet)[1] %in% "FacetWrap") {filenme="attackrate_by_district_facetwrap.png"} else {filenme="attackrate_by_district.png"}
# ggsave(paste0("simul_output/somalia_output/",filenme),width=32,height=16,units="cm")

#' ### Age distribution
# "summary_by_age_grp.csv" | "cases_by_age_gender.csv" | "death_by_age_gender.csv"
summary_by_age_grp=read_csv("somalia_linelist/summary_by_age_grp.csv")
# death_by_age_gender=read_csv("somalia_linelist/death_by_age_gender.csv")
somal_agestr=data.frame(`Age Group`=popF[popF$name %in% "Somalia","age"],
                        value=popF[popF$name %in% "Somalia","2020"]+popM[popM$name %in% "Somalia","2020"])
summary_by_age_grp$total_popul=c(c(sapply(seq(1,19,2),function(x) {sum(somal_agestr$value[c(x,x+1)])})*1e3,
                                   somal_agestr$value[nrow(somal_agestr)]*1e3),
                                 sum(c(sapply(seq(1,19,2),function(x) {sum(somal_agestr$value[c(x,x+1)])})*1e3,somal_agestr$value[nrow(somal_agestr)]*1e3)))
summary_by_age_grp$total_popul_perc=c(N_tot,sum(N_tot))/sum(N_tot)

summary_by_age_grp=summary_by_age_grp %>% mutate(total_tested_perc=`Total Tested`/`Total Tested`[`Age Group` %in% "National"]) %>% 
  pivot_longer(cols=!c("Age Group"))
summary_by_age_grp$`Age Group`=factor(summary_by_age_grp$`Age Group`,levels=unique(summary_by_age_grp$`Age Group`))

# plot
ggplot(subset(summary_by_age_grp,!`Age Group` %in% "National" & !name %in% "total_popul"),aes(x=`Age Group`,y=value)) + 
  geom_bar(stat="identity",color="black",width =0.7) + facet_wrap(~name,scales="free") + theme_bw() + standard_theme + 
  labs(title="Age distribution of cumulative COVID19 metrics (weeks 12-32)",caption = "source: line list") + xlab("") + ylab("")
# SAVE
# ggsave("simul_output/somalia_output/summary_by_age_grp.png",width=32,height=24,units="cm")

# rmarkdown::render("somalia_data_model.R",output_dir = "simul_output/somalia_output/")
