#..........................................................................................
### + ESTIMATION OF MORTALITY DURING THE COVID-19 PANDEMIC IN BANADIR, SOMALIA (2020) + ###
#..........................................................................................

#..........................................................................................
## ----------------- R CODE TO PREPARE DATA AND FIT STATISTICAL MODELS ----------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2021)
                                          # francesco.checchi@lshtm.ac.uk 


#.........................................................................................                            
### Estimating excess surface area across Banadir region

  #...................................   
  ## Identify period of data availability for all cemeteries
    # For cemeteries that were only established during the period, set surface area and graves value to 0 from 1 Jan 2016
    x1 <- cemeteries[cemeteries$status_start == "did not exist", "cemetery"]
    for (i in x1) {
      x2 <- min(subset(obs, cemetery == i & ! is.na(area))$date, na.rm = TRUE)
      obs[obs$cemetery == i & obs$date < x2, "area"] <- 0
      obs[obs$cemetery == i & obs$date < x2, "graves_best"] <- 0
    }
    
    # Complete observations
    x1 <-  obs[complete.cases(obs[, c("area", "graves_best")]), c("cemetery", "time_base", "date", "area", "graves_best")]
      
    # First date on which both area and graves are quantified in all cemeteries
    date_min <- max( aggregate(x1[, "date"], by = list(x1$cemetery), FUN = min)[, 2] )

    # Last date on which both area and graves are quantified in all cemeteries
    date_max <- min( aggregate(x1[, "date"], by = list(x1$cemetery), FUN = max)[, 2] )
         
  #...................................   
  ## Interpolate each cemetery time series (both graves and area)
  
    # Set up output
    obs[, "area_ipol"] <- NA
    obs[, "graves_best_ipol"] <- NA
    
    # For each cemetery...
    for (i in sort(unique(obs$cemetery)) ) {
      
      # identify data
      x1 <- subset(obs, cemetery == i)
    
      # ignore cemeteries without any surface area
      if (all(is.na(x1$area))) {next}
      
      # interpolate area
        # interpolate
        x2 <- approx(x1[, "time_base"], x1[, "area"], method = "linear", 
          xout = seq(min(x1$time_base), max(x1$time_base), by = 1), ties = "ordered", rule = 2)
        x2 <- data.frame(x2)
        colnames(x2) <- c("time_base", "y")
        x2[, "y"] <- round(x2[, "y"], 0)
          
        # plot
        x3 <- f_plot_ipol(x2, obs[obs$cemetery == i, c("cemetery", "date", "time_base", "period_covid", "area")])
          
        # add to dataset
        obs[obs$cemetery == i, "area_ipol"] <- x3[, "predictions"]
      
      # interpolate graves
        # interpolate
        x2 <- approx(x1[, "time_base"], x1[, "graves_best"], method = "linear", 
          xout = seq(min(x1$time_base), max(x1$time_base), by = 1), ties = "ordered", rule = 2)
        x2 <- data.frame(x2)
        colnames(x2) <- c("time_base", "y")
        x2[, "y"] <- round(x2[, "y"], 0)
          
        # plot
        x3 <- f_plot_ipol(x2, obs[obs$cemetery == i, c("cemetery", "date", "time_base", "period_covid", "graves_best")])
          
        # add to dataset
        obs[obs$cemetery == i, "graves_best_ipol"] <- x3[, "predictions"]
      
    }
      

  #...................................   
  ## Calculate weights for each observation (relative grave density)

    # Restrict analysis period to the above boundaries
    obs <- subset(obs, date %in% c((date_min+1):date_max) )
        
    # Calculate new area and graves every day for each cemetery
    obs[, "new_area_ipol"] <- NA
    obs[, "new_graves_best_ipol"] <- NA
    obs <- obs[order(obs[, "cemetery"], obs[, "date"]), ]  

    # For each cemetery...
    for (i in sort(unique(obs$cemetery)) ) {
      
      # identify data
      x1 <- subset(obs, cemetery == i)
      
      # calculate new area and graves while also setting their initial value to 0    
      obs[obs$cemetery == i, "new_area_ipol"] <- c(0, diff(x1[, "area_ipol"]))
      obs[obs$cemetery == i, "new_graves_best_ipol"] <- c(0, diff(x1[, "graves_best_ipol"]))
    }  
    
    # Calculate grave density within new surface area
    obs[, "density"] <- obs[, "new_graves_best_ipol"] / obs[, "new_area_ipol"]
      
      # 'normalise' it by dividing all values by their mean
      obs[, "density_ratio"] <- obs[, "density"] / mean(obs$density, na.rm = TRUE)
      
      # set NaN values to NA (instances of new area being 0)
      obs[, "density_ratio"] <- ifelse(obs[, "new_area_ipol"] == 0, NA, obs[, "density_ratio"])
      
      # Create 'new surface area index', i.e. sum of all new area x their density ratio (weight)
      obs[, "new_area_index"] <- obs[, "new_area_ipol"] * obs[, "density_ratio"] 

  #...................................   
  ## Aggregate each cemetery time series into a single one for Banadir
    
    # Sum all cemetery time series
    bdr <- aggregate(obs[, c("new_area_ipol", "new_area_index")], by = list(obs$date), FUN = sum, na.rm = TRUE)
    colnames(bdr) <- c("date", "new_area_ipol", "new_area_index")

    # Restrict observations to data availability period
    bdr <- subset(bdr, date %in% c((date_min+1):date_max) )

  #...................................   
  ## Define baseline and calculate daily total and excess burials
    
    # Baseline = median of time before the start of the epidemic
    baseline <- median(bdr[bdr$date < date_knot, "new_area_index"])
    
    # Calculate excess new area index every day
    bdr[, "new_area_excess"] <- bdr[, "new_area_index"] - baseline
    
    # Write output
    write.csv(bdr[, c("date", "new_area_ipol", "new_area_index", "new_area_excess")], 
      "out_bdr_daily_new_area.csv", row.names= FALSE)
          
    
  #...................................   
  ## Convert daily to weekly and monthly area increase
  
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
    
    # Weekly new area
      # aggregate
      bdr_w <- aggregate(bdr[, grep("new_", colnames(bdr))], by = bdr[, c("epi_year", "week")], FUN = sum)
    
      # add dates back in (end of week)
      x1 <- aggregate(bdr$date, by = bdr[, c("epi_year", "week")], FUN = max)
      colnames(x1) <- c("epi_year", "week", "date")
      bdr_w <- merge(bdr_w, x1, by = c("epi_year", "week"))

      # write output
      write.csv(bdr_w[, c("epi_year", "week", "date", "new_area_ipol")], 
        "out_bdr_weekly_new_area.csv", row.names= FALSE)
      
    # Monthly new area
      # aggregate
      bdr_m <- aggregate(bdr[, grep("new_", colnames(bdr))], by = bdr[, c("year", "month")], FUN = sum)
      
      # add dates back in (end of month)
      x1 <- aggregate(bdr$date, by = bdr[, c("year", "month")], FUN = max)
      colnames(x1) <- c("year", "month", "date")
      bdr_m <- merge(bdr_m, x1, by = c("year", "month"))
      
      # write output
      write.csv(bdr_m[, c("year", "month", "date", "new_area_ipol")], 
        "out_bdr_monthly_new_area.csv", row.names= FALSE)
      
  #...................................   
  ## Compute total excess area index
  x1 <- sum(bdr[bdr$y == 2020, "new_area_excess"], na.rm = TRUE)
  print(paste("total excess area = ", x1, sep = "") )    
                
  #...................................   
  ## Plot overall time series
    
    # Plot by week
    plot <- ggplot(bdr_w) +
      geom_line(aes(x = date, y = new_area_index), colour = palette_cb[4], size = 1.5, alpha = 0.5) +
      geom_line(aes(x = date, y = new_area_ipol), colour = palette_cb[1], size = 1, alpha = 0.5) +
      geom_line(aes(x = date, y = (baseline*7) ), linetype = "dashed", colour = brewer_pal(palette = "Reds")(9)[6],
        size = 1, alpha = 0.5) +
      annotate(geom = "rect", xmin = date_knot, xmax = max(bdr_w$date), ymin = 0, ymax = Inf,
        fill = brewer_pal(palette = "Reds")(9)[6], alpha = 0.2) +
      scale_y_continuous("new surface area index per week (m^2)" , breaks=seq(0, 1500, by = 150)) +
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
    ggsave("out_overall_weekly_trends_area.png", width = 25, height = 15, units = "cm", dpi = "print")    
             

    
    
  #...................................  
  ## Plot trends in new graves and density index, by cemetery
    # create a date variable for the x axis
    tot[, "date"] <- dmy(paste("1", tot[, "m"], tot[, "y"], sep="/"))
          
      # create breaks for years
      year.breaks <- subset(tot, m==1)[, "date"]

      # create scaling parameters
      tot.wage.max <- max(tot$tot_wage_cereal, na.rm=TRUE)
      tot.goat.max <- max(tot$tot_goat_cereal, na.rm=TRUE)
      
    # plot
      plot <- ggplot(tot, aes(x = date) ) +
      geom_point(mapping = aes(x = date, y = tot_wage_cereal), stat = "identity", colour = "indianred3", size = 1) +
      geom_line(mapping = aes(x = date, y = tot_wage_cereal), stat = "identity", colour = "indianred3") +
      geom_point(mapping = aes(x = date, y = tot_goat_cereal*(tot.wage.max/tot.goat.max) ), colour = "dodgerblue4", size = 1) +
      geom_line(mapping = aes(x = date, y = tot_goat_cereal*(tot.wage.max/tot.goat.max) ), colour = "dodgerblue4") +
      scale_y_continuous("Terms of trade - Kcal cereal equivalent of daily wage", labels = function(x) format(x, big.mark = ",", scientific = FALSE),
        breaks = c(0, 10000, seq(20000, max(tot$tot_wage_cereal, na.rm=TRUE), by = 20000)),
        limits = c(0, 50000),
        sec.axis = sec_axis(~ (./(tot.wage.max/tot.goat.max)), labels = function(x) format(x, big.mark = ",", scientific = FALSE),
          breaks = seq(0, max(tot$tot_goat_cereal, na.rm=TRUE), by=100000), name="Terms of trade - Kcal cereal equivalent of medium goat")) +
      theme_bw() + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm") ) +
      labs(x = "\nmonth" ) +
      geom_vline(xintercept = year.breaks, color="grey50") +
      geom_hline(yintercept = 10000, color="indianred3", linetype="dashed") +
      facet_wrap(~region, ncol=3, scales = "free_y" ) +
      theme(axis.text.x=element_text(angle = 90, vjust=0.5)) +
      scale_x_date("\nmonth - year",
          expand=c(0,0) , minor_breaks=NULL, date_breaks="12 months", date_labels = "%b-%Y") +
      theme(plot.title = element_text(color="grey30"), legend.title = element_text(color="grey30"),
            axis.title.x = element_text(color="grey30"),
            axis.line.y = element_line(color = "indianred3"), axis.ticks.y = element_line(color = "indianred3"),
            axis.text.y = element_text(color = "indianred3"), axis.title.y = element_text(color="indianred3", margin = margin(r = 10) ),
            axis.line.y.right = element_line(color = "dodgerblue4"), axis.ticks.y.right = element_line(color = "dodgerblue4"),
            axis.text.y.right = element_text(color = "dodgerblue4"), axis.title.y.right = element_text(color="dodgerblue4", margin = margin(l = 10) )
           )
          
    # call plot
    print(plot)
        
    