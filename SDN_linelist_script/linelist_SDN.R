# Covid in Sudan - Analysis of line list data
# author: Fabienne Krauer, CMMID LSHTM
# last updated: 2020-12-07

# Housekeeping ---------------------------------------------
rm(list=ls()); set.seed(42)
library(tidyverse); library(readxl); library(binom); library(lubridate); library(openxlsx); date <- Sys.Date()
theme_set(theme_minimal())

# Read helper data  -------------------------------------
currentdir_path=dirname(rstudioapi::getSourceEditorContext()$path); setwd(currentdir_path)
# Population sizes by state
pop_state <- readxl::read_xlsx("Sudan_datadict.xlsx", trim_ws = T, sheet="population")
# translation of state names
states <- readxl::read_xlsx("Sudan_datadict.xlsx", trim_ws = T, sheet="states")
# Translation to english for symptoms
clin_other <- readRDS("clin_other_translation.rds")

# UN WPP population data for Sudan
data(pop); popF <- popF[popF$name=="Sudan",c("name", "age", "2020")]; popF$gender <- "female"
popM <- popM[popM$name=="Sudan",c("name", "age", "2020")]; popM$gender <- "male"
pop_age_gender <- merge(popF, popM, by=intersect(names(popF), names(popM)), all=T)
rm(list=c("popF", "popM", "popMT", "popFT", "pop")); colnames(pop_age_gender)[3] <- "pop"
pop_age_gender$pop <- pop_age_gender$pop*1000; pop_age_gender$age <- factor(pop_age_gender$age, 
                       levels=c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49",
                                "50-54", "55-59", "60-64", "65-69", "70-74", "75-79", "80-84", "85-89", "90-94", "95-99", "100+"))
pop_age_gender$age_group <- ifelse(pop_age_gender$age %in% c("0-4", "5-9"), "0-9",
                            ifelse(pop_age_gender$age %in% c("10-14", "15-19"), "10-19",
                                   ifelse(pop_age_gender$age %in% c("20-24", "25-29"), "20-29",
                                          ifelse(pop_age_gender$age %in% c("30-34", "35-39"),"30-39",
                                                 ifelse(pop_age_gender$age %in% c("40-44", "45-49"),"40-49",
                                                        ifelse(pop_age_gender$age %in% c("50-54", "55-59"), "50-59",
                                                               ifelse(pop_age_gender$age %in% c("60-64", "65-69"), "60-69",
                                                                      ifelse(is.na(pop_age_gender$age), NA, "70+"))))))))
table(pop_age_gender$age, pop_age_gender$age_group)
pop_age_gender <- pop_age_gender %>% dplyr::group_by(gender, age_group) %>% dplyr::summarise(pop=sum(pop))
# Totals by age group (M & F)
pop_age <- pop_age_gender %>% dplyr::group_by(age_group) %>% dplyr::summarise(pop=sum(pop))

# Cleaning / prep -------------------------------------

# Problems:
# the dates are imported as excel numeric dates, which can be converted back but it is just a pain
# when I convert the excel to csv and import it, I lose the arabic text. 
# Read main data set
data <- readxl::read_xlsx("Sudan_Covid_19.xlsx", trim_ws = T)
length(unique(data$NO.))==nrow(data) # TRUE if observation numbers are unique

#rename columns to english, it is easier to work with in R
colnames(data) <- c("no", "state_rep_orig", "age", "gender", "age_group", "state_orig",
                    "locality", "address", "transmission", "date_onset", "date_hosp", "date_rep",
                    "week", "month", "outcome", "date_death", "sample", "date_sample",
                    "sample_type", "sample_no", "result", "clin_fever", "clin_cough", "clin_dyspnea",
                    "clin_throat", "clin_headache", "clin_anosmia", "clin_other", 
                    "job1", "job2", "place_isol", "outcome_det", "place_samp",
                    "sample1_no", "sample1_res", "sample2_no", "sample2_res",
                    "sample3_no", "sample3_res", "sample4_no", "sample4_res",
                    "notes", "nationality", "destination",
                    "date_depart", "countries_transit", "date_transit", "date_arrival",
                    "means_travel", "entry_point", "group")

# Translate states
data <- merge(data, states, by.x="state_rep_orig", by.y="orig", all.x=T)
colnames(data)[which(colnames(data)=="translation")] <- "state_rep"; sort(unique(data$state_rep))
data <- merge(data, states, by.x="state_orig", by.y="orig", all.x=T); colnames(data)[which(colnames(data)=="translation")] <- "state"
rm(states)

# If state is missing, impute state of reporting (all are from Kharthoum)
data$state <- ifelse((is.na(data$state) & !is.na(data$state_rep)) | (data$state=="unknown" & data$state_rep!="unknown"), 
                     data$state_rep, data$state)

# Age/Age group
sort(unique(data$age))
sort(unique(data$age_group))

data$age <- ifelse(data$age=="غير معروف", NA, data$age) # replace unknown with NA
data$age <- ifelse(data$age=="8.3333333333333329E-2", "0.083", data$age) # this seems to be an incorrect importation into R
data$age <- as.numeric(data$age)
data$age_group <- ifelse(data$age_group=="غير معروف", NA, data$age_group)
data$age_group <- recode(data$age_group, "01-09.9"="0-9", "10 - 19.9"="10-19", "20 -29.9" = "20-29", "30 -39.9"="30-39", 
                         "40-49.9"="40-49", "50-59.9" ="50-59" , "60-69.9"="60-69", "70 & More"="70+") 
# This could also be recoded in the excel sheet itself

# check agreement of age and agecats
data %>% dplyr::group_by(age_group) %>% dplyr::summarise(minage=min(age, na.rm=T), meanage=mean(age, na.rm=T), maxage=max(age, na.rm=T),
                   n=n()) # some ages are very high in the oldest category, not very likely
# Gender
data$gender <- ifelse(data$gender=="غير معروف", "unknown",ifelse(data$gender=="ذكر", "male", "female"))
  
# locality/address
sort(unique(data$locality))
data$locality <- ifelse(data$locality=="0" | data$locality=="غير معروف", NA, data$locality) # replace unknown with NA
data$address <- ifelse(data$address=="غير معروف", NA, data$address) # replace unknown with NA
data$address <- gsub("[0-9]{9}", "", data$address) # remove the phone numbers
data$address <- gsub("[0-9]{10}", "", data$address) # remove the phone numbers

# Unify and correct all date formats
sort(unique(data$date_onset))
data$date_onset <- ifelse(data$date_onset=="28/8/2020", "44071", data$date_onset)
data$date_onset <- ifelse(data$date_onset %in% c("لم يشعر", "غيرمعروف", "غير معروف"),NA, data$date_onset)
data$date_onset <- as.Date(as.numeric(data$date_onset), orig=as.Date("1899-12-31"))-1 # convert from Excel numeric to date
summary(data$date_onset)
data$date_hosp <- as.Date(data$date_hosp); summary(data$date_hosp)

# correct missing date of report for obs 28001, which has been removed during data import because of wrong date format
data$date_rep[data$no==28001] <- as.Date("2020-04-19"); data$date_rep <- as.Date(data$date_rep)
summary(data$date_rep)
sort(unique(data$date_death))
data$date_death <- ifelse(data$date_death %in% c(".", "مؤسسة صحية", "لاينطبق", "غير معروف"), NA, data$date_death)
data$date_death <- as.Date(as.numeric(data$date_death), orig=as.Date("1899-12-31"))-1 # convert from Excel numeric to date
summary(data$date_death)

sort(unique(data$date_sample))
data$date_sample <- ifelse(data$date_sample %in% c("لاينطبق", "لا ينطبق"), NA, data$date_sample)
data$date_sample <- as.Date(as.numeric(data$date_sample), orig=as.Date("1899-12-31"))-1 # convert from Excel numeric to date
summary(data$date_sample)
# if result is positive and date of sampling is missing, impute date of reporting
data$date_sample <- dplyr::if_else(is.na(data$date_sample) & data$result %in% c("positive", "Positive") & 
                                     !is.na(data$result) & !is.na(data$date_rep),
                           data$date_rep, data$date_sample)
sort(unique(data$date_depart))
data$date_depart <- ifelse(data$date_depart =="20/3", "43910", data$date_depart)
data$date_depart <- ifelse(data$date_depart %in% c("لاينطبق", "لا ينطبق", "بدون اسم"),NA, data$date_depart)
data$date_depart <- as.Date(as.numeric(data$date_depart), orig=as.Date("1899-12-31"))-1 # convert from Excel numeric to date
summary(data$date_depart)

# Outcome 
sort(unique(data$outcome))
data$outcome <- ifelse(data$outcome=="توفي", "deceased", ifelse(data$outcome=="حي", "alive", NA))
sort(unique(data$outcome_det))
data$outcome_det <- gsub("\\r\\n", "",data$outcome_det) # remove carriage returns
data$outcome_det <- ifelse(data$outcome_det=="هروب", "escaped", 
                       ifelse(data$outcome_det=="عزل", "isolated", 
                              ifelse(data$outcome_det=="شفاء", "recovering", 
                                     ifelse(data$outcome_det=="سالب/توفي", "negative, deceased", 
                                            ifelse(data$outcome_det=="سالب", "negative", 
                                                   ifelse(data$outcome_det=="توفي", "deceased", NA))))))
table(data$outcome_det, data$outcome)
# Sampling
data$sample <- ifelse(data$sample=="نعم", "yes", ifelse(data$sample=="لا", "no", NA))
# Result
sort(unique(data$result)); table(data$result)
data$result <- ifelse(data$result %in% c("negative", "Negative"), "negative",
                      ifelse(data$result %in% c("positive", "Positive"), "positive",
                             ifelse(data$result %in% c("pending", "Pending"), "pending",
                                    ifelse(data$result %in% c("rejected", "Rejected", "REJECTED"), "rejected", NA))))
table(data$result)

# clinical symptoms
sort(unique(data$clin_fever)); data$clin_fever <- ifelse(data$clin_fever=="نعم", "yes", ifelse(data$clin_fever=="لا", "no", NA))
sort(unique(data$clin_cough))
data$clin_cough <- ifelse(data$clin_cough=="نعم", "yes", ifelse(data$clin_cough=="لا", "no", NA)); sort(unique(data$clin_dyspnea))
data$clin_dyspnea <- ifelse(data$clin_dyspnea=="نعم", "yes", ifelse(data$clin_dyspnea=="لا", "no", NA))
sort(unique(data$clin_throat)); data$clin_throat <- ifelse(data$clin_throat=="نعم", "yes", ifelse(data$clin_throat=="لا", "no", NA))
sort(unique(data$clin_headache)); data$clin_headache<-ifelse(data$clin_headache=="نعم", "yes", ifelse(data$clin_headache=="لا", "no", NA))
sort(unique(data$clin_anosmia)); data$clin_anosmia <- ifelse(data$clin_anosmia=="نعم", "yes", ifelse(data$clin_anosmia=="لا", "no", NA))
sort(unique(data$clin_other))

# merge with google english translation for other symptoms
data <- merge(data, clin_other, by=intersect(names(data), names(clin_other)), all.x=T)
rm(clin_other)

# type of transmission (community or contact? Have to check with MOH again)
sort(unique(data$transmission))
data$transmission <- ifelse(data$transmission== "مخالط", "mixed",
                            ifelse(data$transmission== "غير معروف", "unknown",
                                   ifelse(data$transmission==  "غير مخالط", "no contact",
                                          ifelse(data$transmission=="N", "not applicable",
                                                 ifelse(data$transmission=="T", "contact", data$transmission))))) 
# What group does the tested individual belong to?
# travellers and institutions (expats) are different from the suspected cases because they are mostly
# tested pre-emptively
sort(unique(data$group)); data$group <- ifelse(is.na(data$group), "suspected case", data$group)
data$group <- ifelse(data$group=="مؤوسسات", "institution",ifelse(data$group=="نعم", "traveller", data$group))

# Calculate some new variables

# Isoweeks for the different dates
data$isoweek_onset <- isoweek(data$date_onset)
data$isoweek_hosp <- isoweek(data$date_hosp)
data$isoweek_rep <- isoweek(data$date_rep)
data$isoweek_sample <- isoweek(data$date_sample)
data$isoweek_death <- isoweek(data$date_death)

# Months
data$month_sample <- month(data$date_sample)
#data$month_sample <- factor(data$month_sample, 
#                            levels=month.abb[min(unique(data$month_sample), na.rm=T):max(unique(data$month_sample), na.rm=T)],
#                            ordered=T)
data$month_rep <- month(data$date_rep)
# Year
data$year_sample <- year(data$date_sample)
# Time between onset of symptoms and sampling
data$t_onset_sample <- as.numeric(data$date_sample - data$date_onset) 
# Time between onset of symptoms and reporting
data$t_onset_rep <- as.numeric(data$date_rep - data$date_onset) 
# Time between onset of symptoms and hospital admission
data$t_onset_hosp <- as.numeric(data$date_hosp - data$date_onset) 
# Time between onset of symptoms and death
data$t_onset_death <- as.numeric(data$date_death - data$date_onset) 
# Asymptomatic (based on absence of symptoms)
data$asymptomatic <- ifelse(data$clin_fever=="no" & !is.na(data$clin_fever) & data$clin_cough=="no" & !is.na(data$clin_cough) &
                      data$clin_dyspnea=="no" & !is.na(data$clin_dyspnea) & data$clin_throat=="no" & !is.na(data$clin_throat) &
                      data$clin_headache=="no" & !is.na(data$clin_headache) & 
                      (data$clin_other %in% c("لاتوجد" ,"لا") | is.na(data$clin_other)), "yes", "no")
saveRDS(data, paste0("linelist_clean_", date, ".rds"))
write_csv(data,paste0("linelist_clean_", date, ".csv"))

# Results (on maps) -------------------------------------
update <- max(data$date_sample, na.rm=T)
# load Map
map <- readRDS(url("https://biogeo.ucdavis.edu/data/gadm3.6/Rsp/gadm36_SDN_1_sp.rds"))
# if (!require(gpclib)) install.packages("gpclib", type="source")
map <- fortify(map, region="NAME_1")
data=data[,!is.na(names(data))]
## Epidemiological characteristics
# 1. Total number of tests conducted by state
test_total <- data %>% dplyr::filter(!is.na(data$result)) %>% dplyr::group_by(state) %>% 
 dplyr::summarise(ntest=n(), npos=length(no[result=="positive"]), ntest_res=length(no[result %in% c("positive", "negative")]))
test_total[,c("mean","low95","up95")]<-binom.confint(test_total$npos,test_total$ntest_res,method="exact")[,c("mean","lower","upper")]*100
test_total$ntest_cat<-cut(test_total$ntest,breaks=c(0,99,199,499,999,max(test_total$ntest)),labels=c("<100","<200","<500","<1000","1000+"))
test_total <- merge(test_total, pop_state, by="state", all=T)

test_total$ntest_pc <- test_total$ntest*1e5/test_total$population
test_total$ntest_pc_cat <- cut(test_total$ntest_pc, breaks=c(0,9,19,49,99, max(test_total$ntest_pc, na.rm=T)),
                               labels=c("<10","<20", "<50", "<100", "100+"))
test_total_map <- merge(test_total, map, by.x="state", by.y="id", all=T)

fig1a <- ggplot(test_total_map) + scale_fill_brewer(palette="YlOrRd") + coord_equal() + 
 geom_polygon(aes(x=long, y=lat, group=group, fill=ntest_cat), color="grey30") + guides(fill=guide_legend("N tested"))
fig1a

fig1b <- ggplot(test_total_map) + scale_fill_brewer(palette="YlOrRd") + coord_equal() + 
  geom_polygon(aes(x=long, y=lat, group=group, fill=ntest_pc_cat), color="grey30") +
  guides(fill=guide_legend("N tested per 100,000 population"))
fig1b

### ### ### ### ###
### age distribution by result
agedistrib_testresult = data %>% group_by(result,age_group) %>% summarise(n = n()) %>% mutate(freq = n / sum(n))# tally()
ggplot(agedistrib_testresult[grepl("tive",agedistrib_testresult$result),],aes(x=age_group,y=freq,group=result,color=result)) + 
  geom_line() + geom_point() + theme_bw() + standard_theme + scale_y_continuous(breaks=(0:20)/40)
ggsave("SDN_linelist_agedistrib.png",width=30,height=18,units="cm") 

# 2. Contribution of different groups among all tests conducted
group <- data.frame(round(prop.table(table(data$group)),2)*100); colnames(group) <- c("group", "percentage")

# 3. Positivity rate by state and month for suspected cases 
# (omit travellers and institutions)
test_monthly <- data %>% dplyr::filter(data$result %in% c("positive", "negative") & group=="suspected case") %>% 
  dplyr::group_by(state, month_sample) %>% dplyr::summarise(npos=length(no[result=="positive"]), ntested=n())
test_monthly[,c("mean","low95","up95")] <- binom.confint(test_monthly$npos, test_monthly$ntested, 
                                                         method="exact")[,c("mean","lower","upper")]*100
fig2 <- ggplot(test_monthly) +
  geom_line(aes(x=month_sample, y=mean), color="red") +
  geom_ribbon(aes(x=month_sample, ymin=low95, ymax=up95), alpha=0.3, fill="red") +
  scale_x_continuous(breaks=seq(1,10,by=2), labels=month.abb[seq(1,10,by=2)]) +
  facet_wrap(~ state) + 
  xlab(NULL) + ylab("Percentage positives among all tested") +
  theme(panel.grid.minor = element_blank())
fig2 


# 4. cumulative cases by state
cumulinc <- data %>% dplyr::filter(result=="positive") %>% 
  dplyr::group_by(state) %>% 
  dplyr::summarise(cases=n(), deaths=length(no[outcome=="deceased"]))
cumulinc$cases_cat <- cut(cumulinc$cases, 
                          breaks=c(0,99,199,499,999, 4999, max(cumulinc$cases)),
                          labels=c("<100", "<200", "<500", "<1000", "<5000", "5000+"))

cumulinc <- merge(cumulinc, pop_state, by="state", all=T)
cumulinc$cases_pc <- cumulinc$cases*100000/cumulinc$population
cumulinc$cases_pc_cat <- cut(cumulinc$cases_pc, 
                             breaks=c(0,9,19,49,99, max(cumulinc$cases_pc, na.rm=T)),
                             labels=c("<10", "<20", "<50", "<100", "100+"))

cumulinc_map <- merge(cumulinc, map, by.x="state", by.y="id", all=T)

fig3a <- ggplot(cumulinc_map) + 
  scale_fill_brewer(palette="YlOrRd") +
  geom_polygon(aes(x=long, y=lat, group=group, fill=cases_cat), color="grey30") +
  coord_equal() + guides(fill=guide_legend("cases"))

fig3a

fig3b <- ggplot(cumulinc_map) + 
  scale_fill_brewer(palette="YlOrRd") +
  geom_polygon(aes(x=long, y=lat, group=group, fill=cases_pc_cat), color="grey30") +
  coord_equal() + guides(fill=guide_legend("cases per 100,000 population"))

fig3b

# 5. Weekly incidence by state
weeklyinc <- data %>% dplyr::filter(result=="positive") %>% 
                  dplyr::group_by(isoweek_sample, state) %>% 
                  dplyr::summarise(year=year_sample[n()],
                                   cases=n(), deaths=length(no[outcome=="deceased"]))
# impute last date of the iso week for plotting
weeklyinc$date <- ISOweek::ISOweek2date(paste0(weeklyinc$year, "-W",weeklyinc$isoweek_sample, "-7"))

fig4a <- ggplot(weeklyinc[weeklyinc$state=="Khartoum",]) + 
  geom_line(aes(x=date, y=cases, color=state)) + xlab(NULL)

fig4a

fig4b <- ggplot(weeklyinc[weeklyinc$state!="Khartoum",]) + 
  geom_line(aes(x=date, y=cases, color=state)) + xlab(NULL)
fig4b

fig4c <- ggplot(weeklyinc[weeklyinc$state=="Khartoum",]) + 
  geom_line(aes(x=date, y=deaths, color=state)) + xlab(NULL)
fig4c

fig4d <- ggplot(weeklyinc[weeklyinc$state!="Khartoum",]) + 
  geom_line(aes(x=date, y=deaths, color=state)) + xlab(NULL)
fig4d

fig5 <- ggplot(data=weeklyinc) + 
  geom_line(aes(x=date, y=cases)) + 
  xlab(NULL) +
  facet_wrap(~state, scale="free_y", ncol=3)
fig5

# 6. age-specific attack rate
attackrate_age <- data %>% dplyr::filter(data$result %in% c("positive") & group=="suspected case" & !is.na(age_group)) %>% 
                  dplyr::group_by(age_group) %>% dplyr::summarise(ncases=n())
attackrate_age <- merge(attackrate_age, pop_age, by="age_group", all=T)
attackrate_age[,c("mean","low95","up95")] <- binom.confint(attackrate_age$ncases, attackrate_age$pop, method="exact")[,c("mean","lower","upper")]*100

fig6 <- ggplot(attackrate_age) + 
  geom_point(aes(x=as.factor(age_group), y=mean)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=low95, ymax=up95)) +
  ylab("age-specific attack rate (% of population infected)") + xlab("age (years)") 
fig6

# Delays
summary(data$t_onset_rep[data$result=="positive" & data$group=="suspected case"])
summary(data$t_onset_sample[data$result=="positive" & data$group=="suspected case"])
summary(data$t_onset_hosp[data$result=="positive" & data$group=="suspected case"])


## Clinical characteristics

# For the clinical characteristics, subset the data to only positive results and only suspected cases
# because travellers and institutions are sampled differently and the CFR may be biased
subset <- data[data$result=="positive" & !is.na(data$result) & data$group=="suspected case",]

# Overall CFR
prop.table(table(subset$outcome))

# age distribution
fig7 <- ggplot(subset) + geom_histogram(aes(age)) + xlab("age (years)")
fig7


# CFR
CFR_age <- subset %>% 
  dplyr::group_by(age_group) %>% 
  dplyr::summarise(ndead=length(outcome[outcome=="deceased"]),
                  ntotal=n())
CFR_age[,c("mean","low95","up95")] <- binom.confint(CFR_age$ndead, CFR_age$ntotal, method="exact")[,c("mean","lower","upper")]*100

fig8a <- ggplot(CFR_age) + 
  geom_point(aes(x=as.factor(age_group), y=mean)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=low95, ymax=up95)) +
  ylab("case fatality ratio (%)") + xlab("age (years)") 

fig8a

# CFR by month
CFR_month <- subset %>% 
  dplyr::group_by(month_sample) %>% 
  dplyr::summarise(ndead=length(outcome[outcome=="deceased"]),
                   ntotal=n())
CFR_month[,c("mean","low95","up95")] <- binom.confint(CFR_month$ndead, CFR_month$ntotal, method="exact")[,c("mean","lower","upper")]*100

fig8b <- ggplot(CFR_month) +
  geom_line(aes(x=month_sample, y=mean), color="red") +
  geom_ribbon(aes(x=month_sample, ymin=low95, ymax=up95), alpha=0.3, fill="red") +
  scale_x_continuous(breaks=seq(1,10,by=2), labels=month.abb[seq(1,10,by=2)]) +
  xlab(NULL) + ylab("case fatality ratio (%)") +
  theme(panel.grid.minor = element_blank())
fig8b

# CFR by state
CFR_state <- subset %>% 
  dplyr::group_by(state) %>% 
  dplyr::summarise(ndead=length(outcome[outcome=="deceased"]),
                   ntotal=n())
CFR_state[,c("mean","low95","up95")] <- binom.confint(CFR_state$ndead, CFR_state$ntotal, method="exact")[,c("mean","lower","upper")]*100

fig8c <- ggplot(CFR_state) + 
  geom_point(aes(x=as.factor(state), y=mean)) +
  geom_linerange(aes(x=as.factor(state), ymin=low95, ymax=up95)) +
  ylab("overall case fatality ratio (%)") + xlab(NULL) +
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
fig8c      

  
# CFR by gender
CFR_gender <- subset %>% 
  dplyr::group_by(gender) %>% 
  dplyr::summarise(ndead=length(outcome[outcome=="deceased"]),
                   ntotal=n())
CFR_gender[,c("mean","low95","up95")] <- binom.confint(CFR_gender$ndead, CFR_gender$ntotal, method="exact")[,c("mean","lower","upper")]*100

fig8d <- ggplot(CFR_gender) + 
  geom_point(aes(x=as.factor(gender), y=mean)) +
  geom_linerange(aes(x=as.factor(gender), ymin=low95, ymax=up95)) +
  ylab("overall case fatality ratio (%)") + xlab(NULL)
fig8d     


# Percentage men among all cases by age group
gender <- subset %>% filter(!is.na(age_group)) %>% 
  dplyr::group_by(age_group) %>% 
  dplyr::summarise(total=n(),
                   men=length(no[gender=="male" & !is.na(gender)]))
gender[,c("mean","low95","up95")] <- binom.confint(gender$men, gender$total, method="exact")[,c("mean","lower","upper")]*100

gender <- merge(gender, pop_age_gender[pop_age_gender$gender=="male",c("age_group", "pop")], by="age_group", all=T)
colnames(gender)[7] <- "pop_male"
gender <- merge(gender, pop_age, by="age_group", all=T)
gender[,c("popmale_mean","popmale_low95","popmale_up95")] <- binom.confint(gender$pop_male, gender$pop, method="exact")[,c("mean","lower","upper")]*100

fig9 <- ggplot(gender) + 
  geom_point(aes(x=as.factor(age_group), y=popmale_mean, color="general population")) +
  geom_linerange(aes(x=as.factor(age_group), ymin=popmale_low95, ymax=popmale_up95, color="general population")) +
  geom_point(aes(x=as.factor(age_group), y=mean,  color="positive cases")) +
  geom_linerange(aes(x=as.factor(age_group), ymin=low95, ymax=up95 ,color="positive cases")) +
  ylab("% males among total population in age group") + xlab("age (years)") +
  scale_color_manual(breaks=c("general population", "positive cases"), values=c("grey40", "red")) +
  labs(color="")
  #scale_y_continuous(limits=c(0,100))
fig9


# Symptoms by age and outcome
symptoms <- subset %>% 
  dplyr::group_by(age_group, outcome) %>% 
  dplyr::summarise(total=n(),
                   fever=length(no[clin_fever=="yes" & !is.na(clin_fever)]),
                   cough=length(no[clin_cough=="yes" & !is.na(clin_cough)]),
                   dyspnea=length(no[clin_dyspnea=="yes" & !is.na(clin_dyspnea)]),
                   throat=length(no[clin_throat=="yes" & !is.na(clin_throat)]),
                   headache=length(no[clin_headache=="yes" & !is.na(clin_headache)]),
                   anosmia=length(no[clin_anosmia=="yes" & !is.na(clin_anosmia)]),
                   asymptomatic=length(no[asymptomatic=="yes"]))
    
symptoms[,c("fever_mean","fever_low95","fever_up95")] <- binom.confint(symptoms$fever, symptoms$total, method="exact")[,c("mean","lower","upper")]*100
symptoms[,c("cough_mean","cough_low95","cough_up95")] <- binom.confint(symptoms$cough, symptoms$total, method="exact")[,c("mean","lower","upper")]*100
symptoms[,c("dyspnea_mean","dyspnea_low95","dyspnea_up95")] <- binom.confint(symptoms$dyspnea, symptoms$total, method="exact")[,c("mean","lower","upper")]*100
symptoms[,c("throat_mean","throat_low95","throat_up95")] <- binom.confint(symptoms$throat, symptoms$total, method="exact")[,c("mean","lower","upper")]*100
symptoms[,c("headache_mean","headache_low95","headache_up95")] <- binom.confint(symptoms$headache, symptoms$total, method="exact")[,c("mean","lower","upper")]*100
symptoms[,c("anosmia_mean","anosmia_low95","anosmia_up95")] <- binom.confint(symptoms$anosmia, symptoms$total, method="exact")[,c("mean","lower","upper")]*100
symptoms[,c("asymp_mean", "asymp_low95", "asymp_up95")] <- binom.confint(symptoms$asymptomatic, symptoms$total, method="exact")[,c("mean","lower","upper")]*100

# Fever
fig10a <- ggplot(symptoms) + ggtitle("Fever") +
  geom_point(aes(x=as.factor(age_group), y=fever_mean, colour=outcome), position=position_dodge(width=0.2)) + 
  geom_linerange(aes(x=as.factor(age_group), ymin=fever_low95, ymax=fever_up95, colour=outcome), position=position_dodge(width=0.2)) +
  theme(panel.grid.major = element_blank()) +
  ylab("reported (%)") + xlab("age (years)") + scale_y_continuous(limits=c(0,100))

fig10a

# Cough
fig10b <- ggplot(symptoms) +   ggtitle("Cough") +
  geom_point(aes(x=as.factor(age_group), y=cough_mean, colour=outcome), position=position_dodge(width=0.2)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=cough_low95, ymax=cough_up95, colour=outcome), position=position_dodge(width=0.2)) +
  ylab("reported (%)") + xlab("age (years)") + scale_y_continuous(limits=c(0,100))

fig10b

# Dyspnea
fig10c <- ggplot(symptoms) +  ggtitle("Dyspnea") +
  geom_point(aes(x=as.factor(age_group), y=dyspnea_mean, colour=outcome), position=position_dodge(width=0.2)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=dyspnea_low95, ymax=dyspnea_up95, colour=outcome), position=position_dodge(width=0.2)) +
  ylab("reported (%)") + xlab("age (years)") + scale_y_continuous(limits=c(0,100))

fig10c

# throat
fig10d <- ggplot(symptoms) +  ggtitle("Sore throat") +
  geom_point(aes(x=as.factor(age_group), y=throat_mean, colour=outcome), position=position_dodge(width=0.2)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=throat_low95, ymax=throat_up95, colour=outcome), position=position_dodge(width=0.2)) +
  ylab("reported (%)") + xlab("age (years)") + scale_y_continuous(limits=c(0,100))

fig10d

# headache
fig10e <- ggplot(symptoms) +  ggtitle("Headache") +
  geom_point(aes(x=as.factor(age_group), y=headache_mean, colour=outcome), position=position_dodge(width=0.2)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=headache_low95, ymax=headache_up95, colour=outcome), position=position_dodge(width=0.2)) +
  ylab("reported (%)") + xlab("age (years)") + scale_y_continuous(limits=c(0,100))

fig10e

# anosmia
fig10f <- ggplot(symptoms) +  ggtitle("Anosmia/Ageusia") +
  geom_point(aes(x=as.factor(age_group), y=anosmia_mean, colour=outcome), position=position_dodge(width=0.2)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=anosmia_low95, ymax=anosmia_up95, colour=outcome), position=position_dodge(width=0.2)) +
  ylab("reported (%)") + xlab("age (years)") + scale_y_continuous(limits=c(0,100))
fig10f

# proportion asymptomatic
fig10g <- ggplot(symptoms) + 
  geom_point(aes(x=as.factor(age_group), y=asymp_mean, colour=outcome), position=position_dodge(width=0.2)) +
  geom_linerange(aes(x=as.factor(age_group), ymin=asymp_low95, ymax=asymp_up95, colour=outcome), position=position_dodge(width=0.2)) +
  ylab("proportion asymptomatic (%)") + xlab("age (years)") + scale_y_continuous(limits=c(0,100))
fig10g 


save.image(paste0("linelist/ws_", date,".RData"))