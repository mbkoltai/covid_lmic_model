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
    }
    
    # Complete observations
    x1 <-  obs[complete.cases(obs[, c("graves_best")]), c("cemetery", "time_base", "date", "graves_best")]
      
    # First date on which graves are quantified in all cemeteries
    # date_min <- max( aggregate(x1[, "date"], by = list(x1$cemetery), FUN = min)[, 2] )
    date_min <- date(as.Date("1Jan2017", "%d%b%Y") )
      
    # Last date on which graves are quantified in all cemeteries
    date_max <- min( aggregate(x1[, "date"], by = list(x1$cemetery), FUN = max)[, 2] )

         
  #...................................   
  ## Interpolate each cemetery time series and calculate new burials per day
  
    # Set up output
    obs[, c("graves_best_ipol", "new_graves_best_ipol")] <- NA
    
    # For each cemetery...
    for (i in sort(unique(obs$cemetery)) ) {
      
      # identify data
      x1 <- subset(obs, cemetery == i)
    
      # linearly interpolate best estimate
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
        obs[obs$cemetery == i, "new_graves_best_ipol"] <- c(0, diff(x3[, "predictions"]))

    }

          
  #...................................   
  ## Sum each cemetery time series into a single one for Banadir, and calculate corresponding burial rates per 10,000 person-days
    
    # Sum all cemetery time series
    bdr <- aggregate(obs[, c("graves_best_ipol", "new_graves_best_ipol")],
      by = list(obs$date), FUN = sum)
    colnames(bdr) <- c("date", "graves_best_ipol", "new_graves_best_ipol")

    # Restrict observations to start of data availability period
    bdr <- subset(bdr, date >= date_min )
  
    # Calculate alternative burial rates (per 10,000 per day)
    bdr <- merge(bdr, unique(obs[, c("date", "pop_wp2015", "pop_wp2020")]), by = "date")
    bdr[, "br_wp2015"] <- bdr[, "new_graves_best_ipol"] * 10000 / bdr[, "pop_wp2015"]
    bdr[, "br_wp2020"] <- bdr[, "new_graves_best_ipol"] * 10000 / bdr[, "pop_wp2020"]
    
    # Add epidemic time
    bdr <- merge(bdr, unique(ts[, c("date", "time_base", "time_covid")]), by = "date", x.all = TRUE)
    
    # Smooth pre-epidemic burial rate series and extrapolate into epidemic period (as counter-factual)
    bdr[, "br_wp2015_base_s"] <- predict(smooth.spline(subset(bdr, time_covid == 0)[, c("time_base", "br_wp2015")],
      cv = TRUE), bdr$time_base )$y
    bdr[, "br_wp2020_base_s"] <- predict(smooth.spline(subset(bdr, time_covid == 0)[, c("time_base", "br_wp2020")],
      cv = TRUE), bdr$time_base )$y
    
    # Smooth epidemic burial rate series
    bdr[bdr$time_covid > 0, "br_wp2015_covid_s"] <- smooth.spline(subset(bdr, time_covid > 0)[, c("time_covid", "br_wp2015")],
      cv = TRUE)$y
    bdr[bdr$time_covid > 0, "br_wp2020_covid_s"] <- smooth.spline(subset(bdr, time_covid > 0)[, c("time_covid", "br_wp2020")],
      cv = TRUE)$y

  #...................................   
  ## Add further information and save output

    # Add OCHA burials, calculate resulting burial rates per capita and smooth the series
    bdr <- merge(bdr, ocha_bdr[, c("date", "new_graves_ocha")], by = "date", all = TRUE)
    
    bdr[, "br_wp2015_ocha"] <- bdr[, "new_graves_ocha"] * 10000 / bdr[, "pop_wp2015"]
    bdr[, "br_wp2020_ocha"] <- bdr[, "new_graves_ocha"] * 10000 / bdr[, "pop_wp2020"]
    
    bdr[! is.na(bdr$br_wp2015_ocha), "br_wp2015_ocha_s"] <- smooth.spline(subset(bdr, ! is.na(br_wp2015_ocha))[, c("time_base", "br_wp2015_ocha")],
      cv= TRUE)$y
    bdr[! is.na(bdr$br_wp2020_ocha), "br_wp2020_ocha_s"] <- smooth.spline(subset(bdr, ! is.na(br_wp2020_ocha))[, c("time_base", "br_wp2020_ocha")],
      cv= TRUE)$y

    # Add information on the number of images per day
    obs_complete[, "n_images"] <- 1
    x1 <- aggregate(obs_complete$n_images, by = list(obs_complete$date), FUN = sum)
    colnames(x1) <- c("date", "n_images")
    bdr <- merge(bdr, x1, by = "date", all.x = TRUE)
    

  #...................................   
  ## Plot cemetery time series: burials and imagery availability

    # Plot burials by day and cemetery
    plot1 <- ggplot(obs, aes(fill = cemetery, colour = cemetery, y = new_graves_best_ipol, x = date)) +
      geom_bar(position="stack", stat="identity") +
      annotate(geom = "rect", xmin = date_knot, xmax = date_max, ymin = 0, ymax = Inf,
        fill = brewer_pal(palette = "Reds")(9)[6], alpha = 0.2) +
      scale_y_continuous("new burials per day (interpolated)" , expand = c(0, 0), breaks=seq(0, 25, by = 5)) +
      scale_x_date("", minor_breaks=NULL, date_breaks="3 months", limits = c(date_min, date_max), 
        date_labels = "%b-%Y", expand = expansion(add = c(-7, 7)) ) +
      scale_fill_brewer(palette="Dark2") +
      scale_colour_brewer(palette="Dark2") +
      theme_bw() +
      guides(fill = guide_legend(nrow = 1, title = "Cemetery: ")) +
      guides(colour = FALSE) +
      theme(axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_line(color = "grey20"),
        axis.text.y = element_text(color = "grey20", size=11),
        axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
        plot.margin = unit(c(0.5,2,0.5,0.5), "cm"),
        legend.position = "top",
        legend.title = element_text(color = "grey20", size = 10),
        legend.text = element_text(color = "grey20"),
        legend.box = "horizontal"
      )
    
    plot1
    ggsave("out_daily_trends_burials.png", width = 25, height = 15, units = "cm", dpi = "print")    

    # Plot availability of imagery per day by cemetery
    x4<-obs
    obs <- merge(obs, obs_complete[, c("date", "cemetery", "n_images")], by = c("date", "cemetery"), all.x = TRUE)
    obs[, "image_today"] <- obs[, "n_images"]
    obs[is.na(obs$image_today), "image_today"] <- 0
    
    plot2 <- ggplot(subset(obs, image_today > 0), aes(x = date, colour = cemetery, fill = cemetery, y = factor(image_today))) +
      geom_point(size = 3) +
      facet_grid(cemetery ~ .) +
      scale_fill_brewer(palette="Dark2") +
      scale_colour_brewer(palette="Dark2") +
      scale_x_date("", minor_breaks=NULL, date_breaks="3 months", limits = c(date_min, date_max), 
        date_labels = "%b-%Y", expand = expansion(add = c(-7, 7)) ) +
      theme_bw() +
      theme(axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),
        strip.text.y = element_text(color = "grey20", size=10, angle=0, hjust = 0),
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = unit(c(0.5,2,0.5,0.5), "cm"),
        legend.position = "none"
      )
 
    plot2
    ggsave("out_daily_image_avail.png", width = 25, height = 5, units = "cm", dpi = "print")    
    
    # Combine the two above plots
    cowplot::plot_grid(plot1 + theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank(), axis.title.x = element_blank() ), plot2,
      ncol = 1, nrow = 2, labels = c("A", "B"), rel_heights = c(1, 0.4), align = "v", axis ="lrbt")

    ggsave("out_daily_burials_image_combined.png", width = 25, height = 25, units = "cm", dpi = "print")
        

  #...................................   
  ## Plot overall time series: burial rate

    # Prepare data
    x1 <- reshape2::melt(bdr, id.vars = "date", 
      measure.vars = grep("br_", colnames(bdr), value = TRUE) )
    x2 <- data.frame("variable" = unique(x1$variable), "estimate" = c("interpolated (high)",
      "interpolated (low)", "pre-pandemic, smoothed (high)", "pre-pandemic, smoothed (low)",
      "pandemic, smoothed (high)", "pandemic, smoothed (low)", 
      "UN OCHA (high)", "UN OCHA (low)", "UN OCHA, smoothed (high)", "UN OCHA, smoothed (low)") )
    x1 <- merge(x1, x2, by = "variable", x.all = TRUE)
    
    # Plot burial rate per 10,000 by day, by population source, including smoothing - without OCHA
    plot <- ggplot(subset(x1, ! variable %in% grep("ocha", unique(x1$variable), value = TRUE) ), 
      aes(x = date, y = value, colour = estimate, alpha = estimate, 
      size = estimate, linetype = estimate)) +
      geom_step() +
      theme_bw() +
      annotate(geom = "rect", xmin = date_knot, xmax = date_max, ymin = 0, ymax = Inf,
        fill = brewer_pal(palette = "Reds")(9)[6], alpha = 0.2) +
      scale_y_continuous("burial rate per 10,000 person-days" , limits = c(0, 0.12), 
        breaks=seq(0, 0.12, by = 0.02)) +
      scale_x_date("", minor_breaks=NULL, date_breaks="1 month", limits = c(date_min, date_max), 
        date_labels = "%b-%Y", expand = expansion(add = c(-7, 7)) ) +
      scale_colour_manual(values = c(brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8]
        ) ) +
      scale_alpha_manual(values = c(0.1, 0.1, 0.7, 0.7, 0.7, 0.7 ) ) +
      scale_size_manual(values = c(0.5, 0.5, 1.25, 1.25, 1.25, 1.25 ) ) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid", "dotted", "dotted")) +
      theme(axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_line(color = "grey20"),
        axis.text.y = element_text(color = "grey20", size=11),
        axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
        plot.margin = unit(c(0.5,2,0.5,0.5), "cm"),
        legend.position = "bottom"
      )
    
    plot
    ggsave("out_daily_trends_burial_rate.png", width = 25, height = 15, units = "cm", dpi = "print")    
   
 
    # Plot burial rate per 10,000 by day, by population source, including smoothing - with OCHA and only for 2020
    plot <- ggplot(subset(x1, date > date_knot), aes(x = date, y = value, colour = estimate, alpha = estimate, 
      size = estimate, linetype = estimate)) +
      geom_step() +
      theme_bw() +
      scale_y_continuous("burial rate per 10,000 person-days" , limits = c(0, 0.25), 
        breaks=seq(0, 0.25, by = 0.05)) +
      scale_x_date("", minor_breaks=NULL, date_breaks="1 month", limits = c(date_knot, date_max), 
        date_labels = "%b-%Y", expand = expansion(add = c(-7, 7)) ) +
      scale_colour_manual(values = c(brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "Blues")(9)[5],
                                     brewer_pal(palette = "Blues")(9)[8],
                                     brewer_pal(palette = "Blues")(9)[5],
                                     brewer_pal(palette = "Blues")(9)[8]
        ) ) +
      scale_alpha_manual(values = c(0.2, 0.2, 0.7, 0.7, 0.7, 0.7, 0.2, 0.2, 0.7, 0.7 ) ) +
      scale_size_manual(values = c(0.5, 0.5, 1.25, 1.25, 1.25, 1.25, 0.5, 0.5, 1.25, 1.25 ) ) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid", 
        "dotted", "dotted", "solid", "solid", "solid", "solid")) +
      theme(axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_line(color = "grey20"),
        axis.text.y = element_text(color = "grey20", size=11),
        axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
        plot.margin = unit(c(0.5,2,0.5,0.5), "cm"),
        legend.position = "bottom"
      )
    
    plot
    ggsave("out_daily_trends_burial_rate_2020.png", width = 30, height = 15, units = "cm", dpi = "print")    
   
  
    # Plot burial rate per 10,000 by day, by population source, including smoothing - without OCHA and only for 2020
    plot <- ggplot(subset(x1, date > date_knot & ! variable %in% grep("ocha", unique(x1$variable), value = TRUE)),
      aes(x = date, y = value, colour = estimate, alpha = estimate, 
      size = estimate, linetype = estimate)) +
      geom_step() +
      theme_bw() +
      scale_y_continuous("burial rate per 10,000 person-days" , limits = c(0, 0.12), 
        breaks=seq(0, 0.12, by = 0.02)) +
      scale_x_date("", minor_breaks=NULL, date_breaks="1 month", limits = c(date_knot, date_max), 
        date_labels = "%b-%Y", expand = expansion(add = c(-7, 7)) ) +
      scale_colour_manual(values = c(brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "BuGn")(9)[5],
                                     brewer_pal(palette = "BuGn")(9)[8],
                                     brewer_pal(palette = "Blues")(9)[5],
                                     brewer_pal(palette = "Blues")(9)[8],
                                     brewer_pal(palette = "Blues")(9)[5],
                                     brewer_pal(palette = "Blues")(9)[8]
        ) ) +
      scale_alpha_manual(values = c(0.2, 0.2, 0.7, 0.7, 0.7, 0.7, 0.2, 0.2, 0.7, 0.7 ) ) +
      scale_size_manual(values = c(0.5, 0.5, 1.25, 1.25, 1.25, 1.25, 0.5, 0.5, 1.25, 1.25 ) ) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid", 
        "dotted", "dotted", "solid", "solid", "solid", "solid")) +
      theme(axis.title.x = element_text(color="grey20", size=11), 
        axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
        axis.line.y = element_line(color = "grey20"),
        axis.ticks.y = element_line(color = "grey20"),
        axis.text.y = element_text(color = "grey20", size=11),
        axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
        plot.margin = unit(c(0.5,2,0.5,0.5), "cm"),
        legend.position = "bottom"
      )
    
    plot
    ggsave("out_daily_trends_burial_rate_2020_no_ocha.png", width = 30, height = 15, units = "cm", dpi = "print")    
   

  #...................................   
  ## Calculate excess mortality

    # Compute excess mortality multiplier for each population source (actual burial rate / counterfactual baseline)
    bdr[, "rr_wp2015"] <- bdr$br_wp2015_covid_s / bdr$br_wp2015_base_s
    bdr[, "rr_wp2020"] <- bdr$br_wp2020_covid_s / bdr$br_wp2020_base_s
        
    # Set up output - one value for each possible baseline crude death rate and population series
    out <- data.frame("population_source" = c(rep("WorldPop 2015", length(cdr_baseline)), rep("WorldPop 2020", length(cdr_baseline))) ,
      "baseline_cdr" = rep(cdr_baseline, 2) )
    out[, c("baseline_deaths", "excess_deaths", "total_deaths")] <- NA

    # For each possible set of sensitivity values...
    for (i in 1:nrow(out)) {
      
      # select data
      if (out[i, "population_source"] == "WorldPop 2015") {x1 <- bdr[, c("date", "pop_wp2015", "br_wp2015_base_s", "rr_wp2015")]}
      if (out[i, "population_source"] == "WorldPop 2020") {x1 <- bdr[, c("date", "pop_wp2020", "br_wp2020_base_s", "rr_wp2020")]}
      colnames(x1) <- c("date", "pop", "br_base", "rr")
     
      # scale baseline burial rate to baseline death rate
      x1[, "cdr_base"] <- x1[, "br_base"] * out[i, "baseline_cdr"] / x1[x1$date == date_min, "br_base"]
      
      # compute actual death rate during pandemic
      x1[, "cdr_actual"] <- x1[, "cdr_base"] * x1[, "rr"]
      
      # compute excess death rate
      x1[, "cdr_excess"] <- x1[, "cdr_actual"] - x1[, "cdr_base"]
      
      # compute death tolls
      x1 <- subset(x1, date >= date_knot & date <= date_max)
      x1[, c("toll_base",  "toll_excess", "toll_total")] <- x1[, c("cdr_base", "cdr_excess", "cdr_actual")] * x1[, "pop"] / 10000
      out[i, c("baseline_deaths", "excess_deaths", "total_deaths")] <- colSums(x1[, c("toll_base",  "toll_excess", "toll_total")], na.rm = TRUE )
      out[i, c("baseline_deaths", "excess_deaths", "total_deaths")] <- round(out[i, c("baseline_deaths", "excess_deaths", "total_deaths")]/100, digits = 0) * 100
      
    }
    
    # Save output
    write.csv(bdr, "out_bdr_daily_estimates.csv", row.names= FALSE)
    write.csv(out, "out_estimated_deaths.csv", row.names = FALSE)
    