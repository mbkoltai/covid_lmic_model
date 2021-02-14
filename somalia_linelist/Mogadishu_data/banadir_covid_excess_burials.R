#..........................................................................................
### + ESTIMATION OF MORTALITY DURING THE COVID-19 PANDEMIC IN BANADIR, SOMALIA (2020) + ###
#..........................................................................................

#..........................................................................................
## ----------------- R CODE TO PREPARE DATA AND FIT STATISTICAL MODELS ----------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2021)
                                          # francesco.checchi@lshtm.ac.uk 


#.........................................................................................                            
### Estimating excess burials across Banadir region

  #...................................   
  ## Identify period of data availability for all cemeteries
    # For cemeteries that were only established during the period, set graves value to 0 from 1 Jan 2016
    x1 <- cemeteries[cemeteries$status_start == "did not exist", "cemetery"]
    for (i in x1) {
      x2 <- min(subset(obs, cemetery == i & ! is.na(area))$date, na.rm = TRUE)
      obs[obs$cemetery == i & obs$date < x2, "graves_best"] <- 0
      obs[obs$cemetery == i & obs$date < x2, "graves_median"] <- 0
      obs[obs$cemetery == i & obs$date < x2, "graves_lci"] <- 0
      obs[obs$cemetery == i & obs$date < x2, "graves_uci"] <- 0
      
    }
    
    # Complete observations
    x1 <-  obs[complete.cases(obs[, c("graves_best")]), c("cemetery", "time_base", "date", "graves_best")]
      
    # First date on which graves are quantified in all cemeteries
    date_min <- max( aggregate(x1[, "date"], by = list(x1$cemetery), FUN = min)[, 2] )

    # Last date on which graves are quantified in all cemeteries
    date_max <- min( aggregate(x1[, "date"], by = list(x1$cemetery), FUN = max)[, 2] )

         
  #...................................   
  ## Interpolate each cemetery time series
  
    # Set up output
    obs[, c("graves_best_ipol", "graves_lci_ipol", "graves_uci_ipol")] <- NA
    
    # For each cemetery...
    for (i in sort(unique(obs$cemetery)) ) {
      
      # identify data
      x1 <- subset(obs, cemetery == i)
    
      # linearly interpolate best estimate and confidence intervals
      for (j in c("graves_best", "graves_lci", "graves_uci")) {

        # interpolate
        x2 <- approx(x1[, "time_base"], x1[, j], method = "linear", 
          xout = seq(min(x1$time_base), max(x1$time_base), by = 1), ties = "ordered", rule = 2)
        x2 <- data.frame(x2)
        colnames(x2) <- c("time_base", "y")
        x2[, "y"] <- round(x2[, "y"], 0)
        
        # plot
        x3 <- f_plot_ipol(x2, obs[obs$cemetery == i, c("cemetery", "date", "time_base", "period_covid", j)])
        
        # add to dataset
        x2 <- paste(j, "_ipol", sep = "")
        obs[obs$cemetery == i, x2] <- x3[, "predictions"]
      }

    }
      
      
  #...................................   
  ## Sum each cemetery time series into a single one for Banadir
    
    # Sum all cemetery time series
    bdr <- aggregate(obs[, c("graves_best_ipol", "graves_lci_ipol", "graves_uci_ipol")],
      by = list(obs$date), FUN = sum)
    colnames(bdr) <- c("date", "graves_best_ipol", "graves_lci_ipol", "graves_uci_ipol")

    # Restrict observations to data availability period
    bdr <- subset(bdr, date %in% c((date_min+1):date_max) )

  #...................................   
  ## Define baseline and calculate daily total and excess burials
    
    # Calculate new graves every day
    bdr <- bdr[order(bdr[, "date"]), ]
    for (i in c("graves_best_ipol", "graves_lci_ipol", "graves_uci_ipol")) {
      bdr[, paste("new_", i, sep = "")] <- c(0, diff(bdr[, i]))
    }

    # Baseline = median of time before the start of the epidemic
    baseline <- median(bdr[bdr$date < date_knot, "new_graves_best_ipol"])
    
    # Calculate excess new graves every day
    bdr[, "new_graves_excess"] <- bdr[, "new_graves_best_ipol"] - baseline
    
    # Add OCHA burials
    bdr <- merge(bdr, ocha_bdr[, c("date", "new_graves_ocha")], by = "date", all = TRUE)
    
    # Write output
    write.csv(bdr[, c("date", "new_graves_best_ipol", "new_graves_excess", "new_graves_ocha")], 
      "out_bdr_daily_burials.csv", row.names= FALSE)
          
    
  #...................................   
  ## Convert daily to weekly and monthly burials
  
    # Figure out year, month and epidemiological week
    bdr[, "year"] <- year(bdr[, "date"])
    bdr[, "month"] <- month(bdr[, "date"])
    bdr[, "week"] <- isoweek(bdr[, "date"])
      # epidemiological year (for weeks that straddle two calendar years)
      bdr[, "epi_year"] <- bdr[, "year"]
      for (i in 1:nrow(bdr)) {
        if (bdr[i, "week"] == 1 & bdr[i, "month"] == 12) {bdr[i, "epi_year"] <- bdr[i, "year"] + 1} 
        if (bdr[i, "week"] == 52 & bdr[i, "month"] == 1) {bdr[i, "epi_year"] <- bdr[i, "year"] - 1} 
      }
    
    # Weekly new graves
      # aggregate
      bdr_w <- aggregate(bdr[, grep("new_", colnames(bdr))], by = bdr[, c("epi_year", "week")], FUN = sum)
    
      # add dates back in (end of week)
      x1 <- aggregate(bdr$date, by = bdr[, c("epi_year", "week")], FUN = max)
      colnames(x1) <- c("epi_year", "week", "date")
      bdr_w <- merge(bdr_w, x1, by = c("epi_year", "week"))

      # write output
      write.csv(bdr_w[, c("epi_year", "week", "date", "new_graves_best_ipol")], 
        "out_bdr_weekly_burials.csv", row.names= FALSE)
      
    # Monthly new graves
      # aggregate
      bdr_m <- aggregate(bdr[, grep("new_", colnames(bdr))], by = bdr[, c("year", "month")], FUN = sum)
      
      # add dates back in (end of month)
      x1 <- aggregate(bdr$date, by = bdr[, c("year", "month")], FUN = max)
      colnames(x1) <- c("year", "month", "date")
      bdr_m <- merge(bdr_m, x1, by = c("year", "month"))
      
      # write output
      write.csv(bdr_m[, c("year", "month", "date", "new_graves_best_ipol")], 
        "out_bdr_monthly_burials.csv", row.names= FALSE)
      
  #...................................   
  ## Compute total excess deaths
  x1 <- sum(bdr[bdr$y == 2020, "new_graves_excess"], na.rm = TRUE)
  print(paste("total excess burials = ", x1, sep = "") )    
                
  #...................................   
  ## Plot overall time series
    
    # Plot by week
    plot <- ggplot(bdr_w) +
      geom_line(aes(x = date, y = new_graves_best_ipol), colour = palette_cb[4], size = 1.5, alpha = 0.5) +
      geom_line(aes(x = date, y = new_graves_ocha), colour = palette_cb[6], size = 1, alpha = 0.3) +
      geom_point(aes(x = date, y = new_graves_ocha), colour = palette_cb[6], size = 2, alpha = 0.8) +
      geom_line(aes(x = date, y = (baseline*7) ), linetype = "dashed", colour = brewer_pal(palette = "Reds")(9)[6],
        size = 1, alpha = 0.5) +
      annotate(geom = "rect", xmin = date_knot, xmax = max(bdr_w$date), ymin = 0, ymax = Inf,
        fill = brewer_pal(palette = "Reds")(9)[6], alpha = 0.2) +
      scale_y_continuous("new burials per week" , breaks=seq(0, 250, by = 25)) +
      scale_x_date("", minor_breaks=NULL, date_breaks="1 month", limits = c(min(bdr_w$date), max(bdr_w$date)), 
        date_labels = "%b-%Y", expand = expansion(add = c(-7, 7)) ) +
      theme_bw() +
      theme(axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_line(color = "grey20"),
        axis.text.y = element_text(color = "grey20", size=11),
        axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
        plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
      )
    
    plot
    ggsave("out_overall_weekly_trends_burials.png", width = 25, height = 15, units = "cm", dpi = "print")    
             
