# extract results of multiple mcmc fits (different priors / sample size)
fit_dirs=list.dirs("simul_output/somalia",recursive=F)[grepl("fit_",list.dirs("simul_output/somalia",recursive=F))]
l_pairwise_corrs=list()
#######
# start loop
for (k_file in 1:length(fit_dirs)){ # 
mcmc_filename=list.files(fit_dirs[k_file],pattern=".rds",full.names=T); fits_death_scaling <- readRDS(mcmc_filename)
foldertag=gsub("mcmc_","",
        paste0(paste0(paste0(c(min(fits_death_scaling[[1]]$data$date),max(fits_death_scaling[[1]]$data$date)),collapse="_"),"_"),
        paste0("seedsize_",gsub(" ","_",fits_death_scaling[[1]]$priors$seed_size),"_"),
        paste0(names(fits_death_scaling[[1]]$options[c("mcmc_burn_in","mcmc_samples")]),
        as.numeric(fits_death_scaling[[1]]$options[c("mcmc_burn_in","mcmc_samples")]),collapse = "_") ) )
mcmc_foldername=paste0("simul_output/somalia/fit_",foldertag)
if (!dir.exists(mcmc_foldername)) {dir.create(mcmc_foldername); file.copy(mcmc_filename,mcmc_foldername)}
# CDR_vals=c(baseline_daily_burials*1e4/mogadishu_popul,0.1,0.2,0.4); slope_val=linregr$coefficients["agegroup_mean"]
# fitting_date_window
for (k in 1:length(CDR_vals)) {# 
  df_posteriors=fits_death_scaling[[k]]$posterior %>% select(chain,trial,lp,all_of(fitting_params)) %>% 
    mutate(`IFR sympt. infections (%)`=1e2*sapply(ifr_logit_intercept,
          function(x) sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc) ),
         `IFR all infections (%)`=1e2*sapply(ifr_logit_intercept,function(x) 
      sum(inv.logit(x+slope_val*somalia_agegroups_IFR$agegroup_mean)*somalia_agegroups_IFR$agegroup_perc*params$pop[[1]]$y))) %>%
    mutate(CDR=round(CDR_vals[k],3))
  if (k==1){df_posteriors_comb=df_posteriors} else {df_posteriors_comb=rbind(df_posteriors_comb,df_posteriors)} }

### correlations between parameters ----------------
# scatterplots + corrs for all values of CDR
corr_scatterplot=0; if (corr_scatterplot){
png(filename=paste0(mcmc_foldername,"/param_corrs_CDR_scatterplot.png"),units="cm",width=30,height=24,res=500)
p <- ggpairs(df_posteriors_comb %>% select(all_of(fitting_params),CDR) %>% mutate(CDR=factor(CDR)),aes(color=CDR,alpha=0.4)) +
  theme_bw() + standard_theme + labs(caption=foldertag); print(p); dev.off()}
# SAVE # ggsave(paste0(mcmc_foldername,"/param_corrs_CDR_scatterplot.png"),width=30,height=24,units="cm")
# correlations only
pairwise_corrs = bind_rows(lapply(unique(df_posteriors_comb$CDR), function(x) cbind(df_posteriors_comb %>% filter(CDR==x) %>%
  select(all_of(fitting_params)) %>% corrr::correlate() %>% corrr::shave() %>% corrr::stretch() ,CDR=x))) %>%
 mutate(x=factor(x,levels=unique(x)),y=factor(y,levels=unique(y)),highlight_color=ifelse(abs(r)<0.7|is.na(r),0,1)) %>% filter(!is.na(r))
l_pairwise_corrs[[k_file]]= pairwise_corrs %>% mutate(folder_tag=gsub("2020-","",gsub("_burn.*","",foldertag)))
# plot
corr_scalar_plot=0
if (corr_scalar_plot){
ggplot(pairwise_corrs,aes(x,y,fill=r)) + geom_tile(aes(color=factor(highlight_color)),size=2,width=0.96,height=0.96) +
  facet_wrap(~CDR,labeller=labeller(CDR=label_both)) + geom_text(aes(label=round(r,2)),size=5) +
  scale_x_discrete(expand=c(0,0.02)) + scale_y_discrete(expand=c(0,0.02)) +
  scale_fill_gradient2(low="blue",mid="white",high="red",limit=c(-1,1)) + scale_color_manual(values=c(NA,"black"),guide=FALSE) +
  standard_theme + theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +xlab("") + ylab("")
# ## SAVE
ggsave(paste0(mcmc_foldername,"/param_corrs_CDR.png"),width=30,height=24,units="cm")
}
} # end of loop

# plot for all fits
ggplot(bind_rows(l_pairwise_corrs),aes(x,y,fill=r)) + geom_tile(aes(color=factor(highlight_color)),size=1,width=0.96,height=0.96) +
  facet_grid(CDR~folder_tag,labeller=labeller(CDR=label_both)) + geom_text(aes(label=round(r,2)),size=4) +
  scale_x_discrete(expand=c(0,0.02)) + scale_y_discrete(expand=c(0,0.02)) + scale_color_manual(values=c("grey","black"),guide=FALSE) +
  scale_fill_gradient2(low="blue",mid="white",high="red",limit=c(-1,1)) + standard_theme + theme_bw() + 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),axis.text.x=element_text(angle=90,vjust=0.5,hjust=0.95)) + 
  xlab("") + ylab("")
ggsave("simul_output/somalia/all_fits_param_corrs_CDR.png",width=40,height=20,units="cm")