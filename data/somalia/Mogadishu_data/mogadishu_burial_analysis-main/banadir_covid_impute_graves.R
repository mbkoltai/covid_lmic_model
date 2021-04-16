#..........................................................................................
### + ESTIMATION OF MORTALITY DURING THE COVID-19 PANDEMIC IN BANADIR, SOMALIA (2020) + ###
#..........................................................................................

#..........................................................................................
## ----------------- R CODE TO PREPARE DATA AND FIT STATISTICAL MODELS ----------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2021)
                                          # francesco.checchi@lshtm.ac.uk 


#.........................................................................................      
### Imputing missing burial count data
    
  #...................................   
  ## Visualise correlation between graves and surface area, depending on in-filling
    # Graves vs. surface area
      plot <- ggplot(subset(obs, ! is.na(area)) ) +
        geom_point(mapping = aes(colour = infilling, x = area, y = graves), size = 2) +
        scale_colour_manual(values = brewer_pal(palette = "Dark2")(3)) +
        scale_y_continuous("number of graves", minor_breaks=NULL) +
        scale_x_continuous("surface area (square metres)", minor_breaks=NULL) +
        theme_bw() +
        theme(legend.position="top", legend.direction="horizontal") +
        labs(colour = "burial pattern:  ") +
        theme(legend.title = element_text(color="grey20", size=11),
          legend.text = element_text(color="grey20", size=11),
          axis.title.x = element_text(color="grey20", size=11), 
          axis.text.x = element_text(color = "grey20", size=11, angle=315, hjust=0, vjust=0),               
          axis.line.y = element_line(color = "grey20"),
          axis.ticks.y = element_line(color = "grey20"),
          axis.text.y = element_text(color = "grey20", size=11),
          axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
          plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
        )
      plot

    # New graves vs. new surface area
      plot <- ggplot(subset(obs, ! is.na(new_area_since)) ) +
        geom_point(mapping = aes(colour = infilling, x = new_area_since, y = new_graves_since), size = 2) +
        scale_colour_manual(values = brewer_pal(palette = "Dark2")(3)) +
        scale_y_continuous("number of new graves", minor_breaks=NULL, trans = "log2",
          breaks = c(0,50,100,250,500,1000,2000)) +
        scale_x_continuous("new surface area (square metres)", trans = "log2",
          breaks = c(0,50,100,250,500,1000,2000,5000,10000), labels = comma) +
        theme_bw() +
        theme(legend.position="top", legend.direction="horizontal") +
        labs(colour = "Burial pattern:  ") +
        theme(legend.title = element_text(color="grey20", size=11),
         legend.text = element_text(color="grey20", size=11),
         axis.title.x = element_text(color="grey20", size=11), 
         axis.text.x = element_text(color = "grey20", size=11, vjust=0),               
         axis.line.y = element_line(color = "grey20"),
         axis.ticks.y = element_line(color = "grey20"),
         axis.text.y = element_text(color = "grey20", size=11),
         axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
         plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
        )
      plot
      ggsave("out_new_graves_vs_new_area.png", width = 18, height = 10, units = "cm", dpi = "print")    

                
 #...................................   
  ## Visualise distribution of graves and area

    # Graves and area    
      # linear scale
      f_hist("graves", obs, c(NA, NA) )
      f_hist("area", obs, c(NA, NA) ) 
      # log scale
      obs[, "graves_ln"] <- log(obs[, "graves"])
      f_hist("graves_ln", obs, c(NA, NA) ) 
      obs[, "area_ln"] <- log(obs[, "area"])
      f_hist("area_ln", obs, c(NA, NA) ) 

    # New graves and new area    
      # linear scale
      f_hist("new_graves_since", obs, c(NA, NA) )
      f_hist("new_area_since", obs, c(NA, NA) ) 
      # log scale
      obs[, "new_graves_since_ln"] <- log(obs[, "new_graves_since"])
      f_hist("new_graves_since_ln", obs, c(NA, NA) ) 
      obs[, "new_area_since_ln"] <- log(obs[, "new_area_since"])
      f_hist("new_area_since_ln", obs, c(NA, NA) ) 
      
    # Does the ratio of new graves per new area change depending on epidemic period?
      obs[, "area_per_grave"] <- obs[, "new_area_since"] / obs[, "new_graves_since"]
      plot <- ggplot(obs, aes(x = date, y = area_per_grave, colour = cemetery)) + geom_point(size = 2) +
        theme_bw() +
        geom_smooth(method = "lm", aes(x = date, y = area_per_grave)) +
        scale_y_continuous("surface area (m^2) per new grave") +
        scale_x_date("", minor_breaks=NULL, date_breaks="6 months", 
          date_labels = "%b-%Y", expand = expansion(add = c(-7, 7)) ) +
        theme(legend.position = "bottom")
      print(plot)
      ggsave("out_area_per_new_grave_over_time.png", width = 18, height = 10, units = "cm", dpi = "print")    
      
      
  #...................................   
  ## Predictive model of graves as a function of surface area

    # GAMLSS generalised additive mixed model; ln(new surface area) as a growth variable

      # select data (exclude observations with zero increase in surface area)
      x1 <- c("new_graves_since", "cemetery", "new_area_since", "time_base")
      obs_fit <- obs[complete.cases(obs[, x1]), x1]
      obs_fit <- subset(obs_fit, new_area_since > 1)

      # fit model
      fit <- gamlss(new_graves_since ~ pb(log(new_area_since), df = 2) + re(random = ~1 + new_area_since | cemetery),
        sigma.formula=~pb(new_area_since), data = obs_fit, family = NBI)
      summary(fit)
      # plot(fit)

      # examine goodness of fit
      x1 <- data.frame(fitted(fit), obs_fit$new_graves_since)
      colnames(x1) <- c("fitted", "observed")
      plot <- ggplot(x1) +
        geom_point(aes(x = observed, y = fitted), size=2, colour = brewer_pal(palette = "Dark2")(2)[1]) +
        theme_bw() +
          scale_x_continuous("observed", trans = "log2", breaks = c(0,5,25,50,100,250,500,1000,2000,3000) ) +
          scale_y_continuous("fitted", trans = "log2", breaks = c(0,5,25,50,100,250,500,1000,2000,3000)) +
        geom_abline(intercept = 0, slope = 1, colour = brewer_pal(palette = "Dark2")(2)[2] ) +
        theme(axis.title = element_text(colour="grey20")) +
         ggtitle(paste("accuracy of model fit; GAMLSS model to predict ", all.vars(formula(fit))[1] , sep="") ) +
         theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") )
      print(plot)
      ggsave("out_new_graves_imputation.png", width = 18, height = 14, units = "cm", dpi = "print")

      # save final model
      filename <- "out_graves_area_pred_model.csv"
      write.csv(summary(fit), filename)
      
          
    # Generate predictions for observations for which imputation is needed
      # subset of observations with either graves, new graves or surface area non-missing
      obs_complete <- obs[! is.na(obs$graves) | ! is.na(obs$new_graves) | ! is.na(obs$area), 
        c("cemetery", "date", "graves", "new_graves", "new_graves_since", 
          "area", "new_area", "new_area_since", "impute_needed")]

      # generate new variable for new area, since now we need are since the last image, whether it has data or not
      obs_complete[, "new_area_since"] <- obs_complete[, "new_area"]
    
      # # generate likelihood profile for each values of new graves to be imputed
      # obs_complete[, paste("new_graves_pred_", c(1:n_boot), sep = "")] <- NA
      # x1 <- t(f_interval(fit, n_boot, obs_complete[obs_complete$impute_needed == "yes", all.vars(formula(fit))], TRUE))
      # obs_complete[obs_complete$impute_needed == "yes", paste("new_graves_pred_", c(1:n_boot), sep = "")] <- x1

      # point estimates
      obs_complete[, "new_graves_pred"] <- obs_complete[, "new_graves"]
      obs_complete[obs_complete$impute_needed == "yes", "new_graves_pred"] <- 
        round(predict(fit, newdata = obs_complete[obs_complete$impute_needed == "yes", all.vars(formula(fit))], type = "response", allow.new.levels = TRUE), 0)

  #...................................   
  ## Using the point estimates, construct timeline of graves for each cemetery, using imputed values
    
    # Define output
    obs_complete[, "graves_best"] <- NA
    obs_complete <- obs_complete[order(obs_complete[, "cemetery"], obs_complete[, "date"]), ]
        
    # Calculate new graves/graves for missing observations
      # for each cemetery...
        for (i in sort(unique(obs_complete$cemetery)) ) {
  
          # select data from cemetery
          x2 <- subset(obs_complete, cemetery == i)
          
          # for each of the observations...
          for (j in 2:nrow(x2) ) {
              
            # if the new graves value is missing...
            if ( is.na(x2[j, "new_graves_pred"]) ) {
              # calculate based on previous observation within the cemetery
              x3 <- x2[j, "graves"] - x2[j-1, "graves"]
              obs_complete[obs_complete$cemetery == i & obs_complete$date == x2[j, "date"], "new_graves_pred"] <- x3
            }
            
            # if the new graves value is (still) missing...
            if ( is.na(x2[j, "new_graves_pred"]) ) {
              # calculate missing previous graves observation first, then as above
              x4 <- x2[j-2, "graves"] + x2[j-1, "new_graves_pred"]
              x3 <- x2[j, "graves"] - x4
              obs_complete[obs_complete$cemetery == i & obs_complete$date == x2[j, "date"], "new_graves_pred"] <- x3
            }
              
          }
          
          # recreate graves based on imputed values, while also setting grave count to 0 at start of each cemetery time series
          obs_complete[obs_complete$cemetery == i, "graves_best"] <- c(0, cumsum(x2[2:nrow(x2), "new_graves_pred"]) )

        }

      
  # #...................................   
  # ## For a bootstrap sample, construct timeline of graves for each cemetery, using imputed values
  #   
  #   # Define reconstructed graves variable for each bootstrap replicate
  #   obs_complete[, paste("graves_pred_", c(1:n_boot), sep = "")] <- NA  
  #   
  #   # Sort
  #   obs_complete <- obs_complete[order(obs_complete[, "cemetery"], obs_complete[, "date"]), ]
  #   
  #   # Indices of values subject to imputation
  #   x1 <- which(obs_complete$impute_needed == "yes")
  # 
  #   # Prepare output
  #   out_boot <- obs_complete[, c("cemetery", "date")]
  #   out_boot[, paste("graves_boot_", c(1:n_boot), sep = "")] <- NA
  #   
  #   # For each bootstrap run..
  #   for (k in 1:n_boot) {
  #     # control statement
  #     print(paste("now on bootstrap run", k))
  #     
  #     # refresh dataset
  #     obs_boot <- obs_complete
  #     
  #     # select set of imputed values at random for this run, from their likelihood profile  
  #     for (l in x1) {
  #       obs_boot[l, "new_graves"] <- obs_boot[l, paste("new_graves_pred_", sample.int(n_boot, 1), sep = "")]
  #     }
  #     
  #     # calculate new graves/graves for missing observations, using this set of imputed values
  #       # for each cemetery...
  #       for (i in sort(unique(obs_boot$cemetery)) ) {
  # 
  #         # select data from cemetery
  #         x2 <- subset(obs_boot, cemetery == i)
  #           
  #         # for each of the observations...
  #         for (j in 2:nrow(x2) ) {
  #             
  #           # if the new graves value is missing...
  #           if ( is.na(x2[j, "new_graves"]) ) {
  #             # calculate based on previous observation within the cemetery
  #             x3 <- x2[j, "graves"] - x2[j-1, "graves"]
  #             obs_boot[obs_boot$cemetery == i & obs_boot$date == x2[j, "date"], "new_graves"] <- x3
  #           }
  #           
  #           # if the new graves value is (still) missing...
  #           if ( is.na(x2[j, "new_graves"]) ) {
  #             # calculate missing previous graves observation first, then as above
  #             x4 <- x2[j-2, "graves"] + x2[j-1, "new_graves"]
  #             x3 <- x2[j, "graves"] - x4
  #             obs_boot[obs_boot$cemetery == i & obs_boot$date == x2[j, "date"], "new_graves"] <- x3
  #           }
  #             
  #         }
  #         
  #         # recreate graves based on imputed values, while also setting grave count to 0 at start of each cemetery time series
  #         obs_boot[obs_boot$cemetery == i, "graves"] <- c(0, cumsum(x2[2:nrow(x2), "new_graves"]) )
  # 
  #       }
  # 
  #     # attribute boot output
  #     out_boot[, paste("graves_boot_", k, sep = "")] <- obs_boot[, "graves"]
  #     
  #   }
  #   
  #   
  # #...................................   
  # ## Compute 80% confidence interval of grave time series for each cemetery
  #   # Sort output
  #   out_boot[, grep("graves_boot", colnames(out_boot))] <- 
  #     t(apply(out_boot[, grep("graves_boot", colnames(out_boot))], 1, function(x) {sort(x)} ))
  #   
  #   # Select time series corresponding to the median and 80% quantiles of final grave count for each cemetery
  #     # prepare output
  #     out_boot[, c("graves_median", "graves_lci", "graves_uci")] <- NA
  #       
  #     # for each cemetery...
  #     for (i in sort(unique(out_boot$cemetery)) ) {
  #       
  #       # select final date data from cemetery
  #       x2 <- subset(out_boot, cemetery == i)
  #       x3 <- x2[nrow(x2), grep("graves_boot", colnames(x2))]
  #       
  #       # compute quantile and identify the matching column indices
  #       x4 <- findInterval(quantile(x3, c(0.5, 0.1, 0.9)), x3)
  #       out_boot[out_boot$cemetery == i, c("graves_median", "graves_lci", "graves_uci")] <- x2[, x4]
  #     }
  #   
  #     # bind to data
  #     obs_complete <- cbind(obs_complete[, c("cemetery", "date", "graves_best")],
  #       out_boot[, c("graves_median", "graves_lci", "graves_uci")] )
  #   
    
  #...................................   
  ## Merge imputation results with main database
  obs <- merge(obs, obs_complete[, c("cemetery", "date", "graves_best")], by = c("date", "cemetery"), all.x = TRUE)

 