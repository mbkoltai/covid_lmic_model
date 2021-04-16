#..........................................................................................
### + ESTIMATION OF MORTALITY DURING THE COVID-19 PANDEMIC IN BANADIR, SOMALIA (2020) + ###
#..........................................................................................

#..........................................................................................
## ----------------- R CODE TO PREPARE DATA AND FIT STATISTICAL MODELS ----------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2021)
                                          # francesco.checchi@lshtm.ac.uk 


#.........................................................................................                            
### Preparing data for analysis
    
  #...................................
  ## Time series of cemeteries and days (as date and sequential day numbers)
    # Determine start and end of analysis period (based on extent of dataset)
    date_start <- date( min(obs[, "date"], na.rm=TRUE) )
    date_end <- date( max(obs[, "date"], na.rm=TRUE) )

    # Construct time series
    dates <- seq(from = date_start, to = date_end , by = 1)
    days <- data.frame(dates, c(0:(length(dates) - 1) ) )
    colnames(days) <- c("date", "time_base")
    ts <- expand.grid(cemeteries$cemetery, dates)
    colnames(ts) <- c("cemetery", "date")
    ts <- merge(ts, days, by="date")
    ts <- ts[order(ts[, "cemetery"], ts[, "date"]), ]
    ts[, "cemetery"] <- as.character(ts[, "cemetery"])

    # Create period variables
      # for added growth model
        # backward version - baseline time
        x1 <- rep(seq(-as.integer(date_end - date_start), 0, by = 1 ),
                  times = length(unique(obs$cemetery))
                 )
        ts[, "time_base_rev"] <- x1        
    
        # forward version - epidemic time
        x1 <- rep(c(rep(0, times = (date_knot - date_start) ),
                    seq(1, as.integer(date_end - date_knot + 1), by = 1 )
                    ), 
                  times = length(unique(obs$cemetery))
                 )
        ts[, "time_covid"] <- x1
        
        # backward version - epidemic time
        x1 <- rep(c(rep(-as.integer(date_end - date_knot), times = (date_knot - date_start) ),
                    seq(-as.integer(date_end - date_knot), 0, by = 1 )
                    ), 
                  times = length(unique(obs$cemetery))
                 )
        ts[, "time_covid_rev"] <- x1
        
        # knot version - baseline time
        x1 <- rep(c(seq(-as.integer(date_knot - date_start), 0, by = 1 ),
                    seq(1, as.integer(date_end - date_knot), by = 1)
                    ),
                  times = length(unique(obs$cemetery))
                 )
        ts[, "time_base_knot"] <- x1        
        
        # knot version - epidemic time
      
      ts[, "period_covid"] <- ifelse(ts$date < date_knot, "baseline", "epidemic")
      ts[, "period_covid"] <- as.factor(ts$period_covid)
      

  #...................................    
  ## Prepare main observations dataset
    
    # Add 1 grave or 1 square metre to graves or area observations that equal 0 (so as to enable logging)
    obs[, "graves"] <- ifelse(obs$graves == 0, 1, obs$graves)  
    obs[, "area"] <- ifelse(obs$area == 0, 1, obs$area)
    obs[, "new_graves_since"] <- ifelse(obs$new_graves_since == 0, 1, obs$new_graves_since)  
    obs[, "new_area_since"] <- ifelse(obs$new_area_since == 0, 1, obs$new_area_since)
    
    # Add secondary variables
      # days since previous observation (image)
      obs <- obs[order(obs[, "cemetery"], obs[, "date"]), ]
      x1 <- c()    
      for (i in sort(unique(obs$cemetery)) ) {
        x2 <- subset(obs, cemetery == i)
        x1 <- c(x1, 0, as.integer(diff(x2$date)) )
          
      }
      obs[, "days_since"] <- x1


    # Merge dataset with time series
    obs <- merge(ts, obs, by = c("cemetery", "date"), all = TRUE)
    obs[, "date"] <- date(obs$date)

    # Factorise variables as needed
      for (i in colnames(obs)) {
        if (class(obs[, i]) == "character" ) { obs[, i] <- as.factor(obs[, i])  }
      }
    
    # Prepare additional potential predictors
      
      # days since first observation (image)
      obs <- obs[order(obs[, "cemetery"], obs[, "date"]), ]
      x1 <- c()    
      for (i in sort(unique(obs$cemetery)) ) {
        x2 <- subset(obs, cemetery == i)
        x1 <- c(x1, x2$date - min(x2$date, na.rm = TRUE) )
      }
      obs[, "days_since_start"] <- x1
      
      # starting surface area
      obs[, "area_start"] <- NA
      for (i in sort(unique(obs$cemetery)) ) {
        x2 <- subset(obs, cemetery == i)
        x3 <- min(x2$area, na.rm = TRUE)
        if (is.infinite(x3)) {x3 <- NA}
        obs[obs$cemetery == i, "area_start"] <- x3
      }

      # whether we are in the epidemic period or not
      obs[, "epidemic"] <- FALSE
      obs[obs$date >= date_knot, "epidemic"] <- TRUE
          
    # Categorise potential predictors
      # days since first image in cemetery time series
      f_hist("days_since_start", obs, c(NA, NA))
      obs[, "days_since_start_cat"] <- cut(obs[, "days_since_start"], breaks = c(-1, 1500, 2000), labels = c("< 1500", ">= 1500") )
      table(obs[, "days_since_start_cat"])
      
      # starting surface area of cemetery
      f_hist("area_start", obs, c(NA, NA))
      obs[, "area_start_cat"] <- cut(obs[, "area_start"], breaks = c(0, 2, 100000), labels = c("0", "> 0") )
      table(obs[, "area_start_cat"])

    
  #...................................
  ## Prepare population denominators and merge them into observations dataset

    # Smooth population across time series
      population[, "date"] <- ymd(paste(population$year, "-", population$month, "-15", sep = "") )
      population <- merge(population, days, by = "date", all.x = TRUE)
      population[, "time_base"] <- as.integer(population$date - date_start)
      x1 <- predict(smooth.spline(population[, "time_base"], y = population[, "pop_wp2015"], spar = 0.2),
       seq(min(days$time_base, na.rm = TRUE), max(days$time_base, na.rm = TRUE), by = 1) )
      x2 <- predict(smooth.spline(population[, "time_base"], y = population[, "pop_wp2020"], spar = 0.2),
       seq(min(days$time_base, na.rm = TRUE), max(days$time_base, na.rm = TRUE), by = 1) )
      x1 <- data.frame(x1$x, round(x1$y, digits = 0), round(x2$y, digits = 0 ))
      colnames(x1) <- c("time_base", "pop_wp2015", "pop_wp2020")

    # Plot population trends
      x2 <- reshape2::melt(x1, id = "time_base")
      colnames(x2) <- c("time_base", "series", "population")
      x2[, "estimate"] <- ifelse(x2[, "series"] == "pop_wp2020", "WorldPop (2020)", "WorldPop (2015)")
      x2 <- merge(x2, days, by = "time_base", all = TRUE)
      plot <- ggplot(x2) +
        geom_line(aes(x = date, y = population, colour = estimate, linetype = estimate), size = 1.5) +
        scale_colour_manual(values = brewer_pal(palette = "BuGn")(9)[c(5,8)]) +
        scale_linetype_manual(values = c("solid", "twodash")) +
        scale_x_date("", minor_breaks=NULL, date_breaks="4 months", date_labels = "%b-%Y" ) +
        scale_y_continuous(labels = comma) +
        theme_bw() +
        labs (colour = "Estimate: ") +
        guides(linetype = FALSE) +
        theme(legend.position="top", legend.direction="horizontal") +
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
      ggsave("out_population_over_time.png", width = 18, height = 10, units = "cm", dpi = "print")

    # Merge with observations
      obs <- merge(obs, x1[, c("time_base", "pop_wp2015", "pop_wp2020")],
        by = "time_base", all.x = TRUE)

    
  #...................................    
  ## Merge cemetery meta-data with main dataset
    obs <- merge(obs, cemeteries[, c("cemetery", "cemetery_id", "district", "ownership")], by = "cemetery", all.x = TRUE)
      # convert cemetery to factor (needed for random effect)
      obs[, "cemetery"] <- as.factor(obs[, "cemetery"])
      obs[, "cemetery_id"] <- as.factor(obs[, "cemetery_id"])
    

  #...................................    
  ## Prepare alternative datasets of mortality
    
    # OCHA data
    ocha <- subset(ocha, group == "total")
      
      # totals for Banadir region
      ocha_bdr <- aggregate(ocha$graves_ocha, by = list(ocha[, "date"]), FUN = sum)
      colnames(ocha_bdr) <- c("date", "new_graves_ocha")
      
      # add time units
      ocha[ "date"] <- date(ocha[, "date"])
      ocha_bdr[, "week"] <- isoweek(ocha_bdr[, "date"])
      ocha_bdr[, "month"] <- month(ocha_bdr[, "date"])
      ocha_bdr[, "year"] <- year(ocha_bdr[, "date"])
        # epidemiological year (for weeks that straddle two calendar years)
        ocha_bdr[, "epi_year"] <- ocha_bdr[, "year"]
        for (i in 1:nrow(ocha_bdr)) {
          if (ocha_bdr[i, "week"] == 1 & ocha_bdr[i, "month"] == 12) {ocha_bdr[i, "epi_year"] <- ocha_bdr[i, "year"] + 1} 
          if (ocha_bdr[i, "week"] == 52 & ocha_bdr[i, "month"] == 1) {ocha_bdr[i, "epi_year"] <- ocha_bdr[i, "year"] - 1} 
        }      
      ocha_bdr[, "date"] <- as.Date(ocha_bdr$date)
        
      # weekly totals
        # aggregate
        ocha_bdr_w <- aggregate(ocha_bdr[, grep("new_", colnames(ocha_bdr))], by = ocha_bdr[, c("epi_year", "week")], FUN = sum)
        colnames(ocha_bdr_w) <- c("epi_year", "week", "new_graves_ocha")
        
        # add dates back in (end of week)
        x1 <- aggregate(ocha_bdr$date, by = ocha_bdr[, c("epi_year", "week")], FUN = max)
        colnames(x1) <- c("epi_year", "week", "date")
        ocha_bdr_w <- merge(ocha_bdr_w, x1, by = c("epi_year", "week"))

      # monthly totals
        # aggregate
        ocha_bdr_m <- aggregate(ocha_bdr[, grep("new_", colnames(ocha_bdr))], by = ocha_bdr[, c("year", "month")], FUN = sum)
        colnames(ocha_bdr_m) <- c("year", "month", "new_graves_ocha")
          
        # add dates back in (end of month)
        x1 <- aggregate(ocha_bdr$date, by = ocha_bdr[, c("year", "month")], FUN = max)
        colnames(x1) <- c("year", "month", "date")
        ocha_bdr_m <- merge(ocha_bdr_m, x1, by = c("year", "month"))      
        
    # Barakaat committee data
    barakaat_committee <- subset(barakaat_committee, group == "all" & cause == "all")
    barakaat_committee[, "month"] <- month(barakaat_committee[, "date"])
    barakaat_committee[, "year"] <- year(barakaat_committee[, "date"])
    barakaat_committee_m <- aggregate(barakaat_committee$graves_bdc, 
      by = barakaat_committee[, c("year", "month")], FUN = sum)
    colnames(barakaat_committee_m) <- c("year", "month", "new_graves_bdc")
    barakaat_committee_m[, "cemetery"] <- "Barakaat 1 + 2"
        