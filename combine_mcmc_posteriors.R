posterior_CI95_comb=rbind(posterior_CI95 %>% mutate(fit_type="deaths 7x scaled"),
                          posterior_CI95_noscaling%>% mutate(fit_type="deaths not scaled"))

ggplot(posterior_CI95_comb %>% filter(!name %in% c("ifr_logit_intercept","IFR sympt. infections (%)")),
       aes(x=fit_type,group=fit_type,color=fit_type)) + scale_x_discrete(expand = expansion(1.15,0)) + 
  geom_linerange(aes(ymin=ci95_low,ymax=ci95_up),position=position_dodge2(width=0.25),alpha=0.3,size=3) +
  geom_linerange(aes(ymin=ci50_low,ymax=ci50_up),position=position_dodge2(width=0.25),alpha=0.5,size=3) +
  geom_point(aes(y=mean),pch="-",size=10,color="black") + facet_wrap(~name,scales="free") + theme_bw() + standard_theme + xlab("") + 
  ylab("mean (CI50, CI95)") + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank()) + 
  geom_text(aes(x=fit_type,y=mean,label=round(mean,2)),nudge_x=0.3,color="black") + labs(color="")
#
ggsave(paste0("simul_output/somalia_output/posteriors_mean_CIs_deaths_scaled",round(scale_factor),"fold.png"),
       width=30,height=18,units="cm")

# plot IFR estimates
i_col=2:5; 
ifr_estimates=cbind(somalia_agegroups_IFR %>% select(agegroup_mean,ifr_mean),
                    data.frame(sapply(as.numeric(array(subset(posterior_CI95,name %in% "ifr_logit_intercept")[,i_col])), function(x)
                      inv.logit(x + 0.1*somalia_agegroups_IFR$agegroup_mean)*params$pop[[1]]$y )) ) %>% rename(`estimate from literature`=ifr_mean)
colnames(ifr_estimates)[i_col+1]<-colnames(posterior_CI95)[i_col]
# ifr_estimates = ifr_estimates %>% pivot_longer(!agegroup_mean) %>% mutate(type=ifelse(grepl("estimate",name),"data","fit"))
# plot
ifr_estimates_comb = rbind(ifr_estimates %>% mutate(fit_type="deaths 7x scaled"),
                           ifr_estimates_noscaling %>% mutate(fit_type="deaths not scaled"))
ggplot(ifr_estimates_comb,aes(x=agegroup_mean)) + 
  geom_line(aes(y=`estimate from literature`*1e2,group=fit_type),size=1.05) + geom_point(aes(y=`estimate from literature`*1e2,group=fit_type),size=2) + 
  geom_line(aes(y=median*1e2,group=fit_type,color=fit_type),size=1.05) + geom_point(aes(y=median*1e2,group=fit_type,color=fit_type),size=2) + 
  geom_ribbon(aes(ymin=ci95_low*1e2,ymax=ci95_up*1e2,group=fit_type,fill=fit_type),alpha=0.2) + theme_bw() + standard_theme + 
  scale_x_continuous(breaks=2.5+(0:16)*5) + scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x,n=12),
                                                          labels=scales::trans_format("log10", scales::math_format(10^.x)) ) +
  xlab("agegroup median (year)") + ylab("IFR %") + theme(axis.text.x=element_text(vjust=0.5,size=12),axis.text.y=element_text(size=12))
# SAVE
ggsave(paste0("simul_output/somalia_output/ifr_mcmc_estimates_deaths_NOTscaled.png"),width=20,height=12,units="cm")