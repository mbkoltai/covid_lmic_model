#..........................................................................................
### + ESTIMATION OF MORTALITY DURING THE COVID-19 PANDEMIC IN BANADIR, SOMALIA (2020) + ###
#..........................................................................................

#..........................................................................................
## ----------------- R CODE TO PREPARE DATA AND FIT STATISTICAL MODELS ----------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2021)
                                          # francesco.checchi@lshtm.ac.uk 


#..........................................................................................
### Preparatory steps

  #...................................      
  ## Install or load required R packages
    
    # List of required packages
    x1 <- c("ggplot2", "scales", "readxl", "data.table", "influence.ME", "lme4", "nlme", "broom.mixed", "lubridate", 
            "RColorBrewer", "lattice", "zoo", "car", "influence.ME", "MASS", "mgcv", "glmmTMB", "gamlss", "ggpubr", "gtools")
    
    # Install any packages not yet installed
    x2 <- x1 %in% row.names(installed.packages())
    if (any(x2 == FALSE)) { install.packages(x1[! x2]) }
    
    # Load all packages    
    lapply(x1, library, character.only = TRUE)
    

  #...................................      
  ## Starting steps

    # Clean up from previous code / runs
    rm(list=ls(all=TRUE) )
  
    # Set font
    windowsFonts(Arial=windowsFont("Arial"))

    # Set working directory to where this file is stored
    current_path = rstudioapi::getActiveDocumentContext()$path 
    setwd(dirname(current_path ))
    print( getwd() )
    
    # Initialise random numbers
    set.seed(123)
    
    # Colour-blind palette for graphing
    palette_cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    

#.........................................................................................
### Specify parameters
    
    # Candidate dates at which baseline period ends and epidemic epidemic period begins
    date_knot_options <- date(as.Date(c("1Jan2020", "1Feb2020", "1Mar2020", "1Apr2020"), "%d%b%Y") )
      
      # Main date for analysis
      date_knot <- date_knot_options[1]
    
    # # Number of runs for bootstrapping
    # n_boot <- 1000
    #   
    # # Number of folds for cross-validation
    # k_folds <- NA   # NA means leave-one-out cross-validation
    
    # Range of baseline crude death rate (expressed as deaths per 10,000 per day)
    cdr_baseline <- seq(0.20, 0.60, by = 0.05)

#.........................................................................................
### Bespoke functions
    
source("banadir_covid_bespoke_functions.R")
        
#.........................................................................................
### Reading in required files

  #...................................      
  ## Variable dictionary
  dict <- read_excel("banadir_covid_burial_data.xlsx", sheet = "dictionary")
    # remove tibble nonsense
    dict <- as.data.frame(dict)

  #...................................      
  ## Read in all the datasets needed
    # Which datasets
    x1 <- c("obs", "cemeteries", "population", "ocha", "barakaat_committee")
    
    # For each dataset...
    for (i in x1) {
      # read in
      assign(i, read_excel("banadir_covid_burial_data.xlsx", sheet = i) )
        # remove tibble nonsense
        assign(i, as.data.frame(get(i)) )
      # only keep needed columns
      x2 <- subset(dict, sheet == i)[, "use"]
      x2 <- which(! x2 %in% c("no") )
      x3 <- get(i)[, x2]
      assign(i, x3)
    }
    
    # Exclude any observations not to be analysed
      # for burial observations dataset
      obs <- subset(obs, include == "yes")
      obs <- obs[, ! colnames(obs) %in% "include"]
      
      # for cemetery dataset
      cemeteries <- subset(cemeteries, include == "yes")
      cemeteries <- cemeteries[, ! colnames(cemeteries) %in% "include"]
            
#.........................................................................................                            
### Preparing data for analysis
    
source("banadir_covid_prepare_data.R")
        
    
#.........................................................................................      
### Imputing missing burial count data
   
source("banadir_covid_impute_graves.R")                         
    

#.........................................................................................                            
### Estimating excess burials

source("banadir_covid_excess_burials.R")
      
            
#.........................................................................................      
### Describing trends in graves and other data characteristics
  
  #...................................   
  ## Plot trends in graves over time, by period and cemetery
    # Preparatory steps
      # data needed
      x1 <- obs[! is.na(obs$graves_best), c("date", "cemetery", "graves_best", "area", "period_covid")]

    # Draw plot
      plot <- ggplot(x1, aes(x = date, y = graves_best) ) +
              geom_line(linetype = "dashed", size = 0.7,  colour = brewer_pal(palette = "Dark2")(2)[2] ) +
              geom_point(aes(colour = period_covid), shape = 15, size = 2) +
              geom_line(aes(colour = period_covid), size = 1 ) +
              scale_colour_manual(values = brewer_pal(palette = "Dark2")(2)) +
              scale_y_continuous("number of graves") +
              theme_bw() +
              facet_wrap(~cemetery, nrow=5, scales = "free_y") +
              guides(fill = FALSE) +
              theme(legend.position="bottom", legend.direction="horizontal") +
              scale_x_date("", minor_breaks=NULL, date_breaks="4 months", date_labels = "%b-%Y" ) +
              labs(colour = "Period:  ") +
              theme(legend.title = element_text(color="grey20", size=11),
                    strip.text.x = element_text(color="grey20", size=11),
                    legend.text = element_text(color="grey20", size=11),
                    axis.title.x = element_text(color="grey20", size=11), 
                    axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
                    axis.line.y = element_line(color = "grey20"),
                    axis.ticks.y = element_line(color = "grey20"),
                    axis.text.y = element_text(color = "grey20", size=11),
                    axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
                    plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
                    )
      plot
      ggsave("out_burials_over_time.png", width = 18, height = 18, units = "cm", dpi = "print")    

   
  #...................................   
  ## Plot trends in burial rate over time, by cemetery
    # Preparatory steps
      # data needed
      x1 <- obs[is.na(obs$graves_best) == FALSE, c("date", "cemetery", "graves_best")]
      
      # generate burial rate since the previous observation
      x1 <- x1[order(x1[, "cemetery"], x1[, "date"]), ]
      x2 <- c()
      for (i in sort(unique(x1$cemetery) ) ) {
        x3 <- subset(x1, cemetery == i)
        x2 <- c(x2, NA, diff(x3$graves_best) / as.integer(diff(x3$date)) )        

      }

      x1[, "burial_rate"] <- x2
      
      # lag burial rate so as to prepare data for step graph
      x2 <- c()
      for (i in sort(unique(x1$cemetery) ) ) {
        x4 <- subset(x1, cemetery == i)
        x2 <- c(x2, x4[-1, "burial_rate"], NA)  
      }

      x1[, "burial_rate_lag"] <- x2
      x1[, "burial_rate_lag"] <- ifelse(is.na(x1[, "burial_rate_lag"]), x1[, "burial_rate"], x1[, "burial_rate_lag"])
  
      
    # Draw plot
      plot <- ggplot(x1, aes(x = date, y = burial_rate_lag) ) +
              geom_step(size = 1, colour = brewer_pal(palette = "Dark2")(2)[1] ) +
              annotate(geom = "rect", xmin = date_knot, xmax = max(x1$date), ymin = 0, ymax = Inf,
                fill = brewer_pal(palette = "Reds")(9)[6], alpha = 0.3) +
              scale_y_continuous("mean new graves per day", limits = c(0, NA) ) +
              theme_bw() +
              facet_wrap(~cemetery, ncol=2, scales = "free_y") +
              scale_x_date("", minor_breaks=NULL, date_breaks="4 months", date_labels = "%b-%Y" ) +
              theme(strip.text.x = element_text(color="grey20", size=11),
                    axis.title.x = element_text(color="grey20", size=11), 
                    axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
                    axis.line.y = element_line(color = "grey20"),
                    axis.ticks.y = element_line(color = "grey20"),
                    axis.text.y = element_text(color = "grey20", size=11),
                    axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
                    plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
                    )
      plot
 
      ggsave("out_burial_rate_over_time_wide.png", width = 26, height = 18, units = "cm", dpi = "print")    
      ggsave("out_burial_rate_over_time_long.png", width = 20, height = 20, units = "cm", dpi = "print")    
      
   
  #...................................   
  ## Plot trends in rate of surface area expansion over time, by cemetery
    # Preparatory steps
      # data needed
      x1 <- obs[is.na(obs$area) == FALSE, c("date", "cemetery", "area")]
      
      # generate area increase rate since the previous observation
      x1 <- x1[order(x1[, "cemetery"], x1[, "date"]), ]
      x2 <- c()
      for (i in sort(unique(x1$cemetery) ) ) {
        x3 <- subset(x1, cemetery == i)
        x2 <- c(x2, NA, diff(x3$area) / as.integer(diff(x3$date)) )        

      }

      x1[, "area_rate"] <- x2
      
      # lag area increase rate so as to prepare data for step graph
      x2 <- c()
      for (i in sort(unique(x1$cemetery) ) ) {
        x4 <- subset(x1, cemetery == i)
        x2 <- c(x2, x4[-1, "area_rate"], NA)  
      }

      x1[, "area_rate_lag"] <- x2
      x1[, "area_rate_lag"] <- ifelse(is.na(x1[, "area_rate_lag"]), x1[, "area_rate"], x1[, "area_rate_lag"])
  
      
    # Draw plot
      plot <- ggplot(x1, aes(x = date, y = area_rate_lag) ) +
              geom_step(size = 1, colour = brewer_pal(palette = "Dark2")(2)[1] ) +
              annotate(geom = "rect", xmin = date_knot, xmax = max(x1$date), ymin = 0, ymax = Inf,
                fill = brewer_pal(palette = "Reds")(9)[6], alpha = 0.3) +
              scale_y_continuous("mean area increase per day", limits = c(0, NA) ) +
              theme_bw() +
              facet_wrap(~cemetery, ncol=2, scales = "free_y") +
              scale_x_date("", minor_breaks=NULL, date_breaks="4 months", date_labels = "%b-%Y" ) +
              theme(strip.text.x = element_text(color="grey20", size=11),
                    axis.title.x = element_text(color="grey20", size=11), 
                    axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
                    axis.line.y = element_line(color = "grey20"),
                    axis.ticks.y = element_line(color = "grey20"),
                    axis.text.y = element_text(color = "grey20", size=11),
                    axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
                    plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
                    )
      plot
 
      ggsave("out_area_rate_over_time_wide.png", width = 26, height = 18, units = "cm", dpi = "print")    
      ggsave("out_area_rate_over_time_long.png", width = 20, height = 20, units = "cm", dpi = "print")    
      
  #...................................   
  ## Describe data for report
    # Create table
      # each cemetery
      out_descr <- data.frame(unique(obs$cemetery))
      colnames(out_descr) <- "cemetery" 
      
      # statistics for each cemetery...
      out_descr[, c("district", "status_start", "status_end", "n_obs", "start_date", "end_date", 
        "graves_cum", "area_start", "area_end")] <- NA
      for (i in 1:length(out_descr$cemetery) ) {
        # select data
        x1 <- subset(obs, cemetery == out_descr[i, "cemetery"] & image_today == 1 )
        # district
        out_descr[i, "district"] <- unique(as.character(x1$district) )
        # status at start of analysis period
        out_descr[i, "status_start"] <- cemeteries[cemeteries$cemetery == out_descr[i, "cemetery"], "status_start"]
        # status at end of analysis period
        out_descr[i, "status_end"] <- cemeteries[cemeteries$cemetery == out_descr[i, "cemetery"], "status_end"]
        # number of observations per cemetery
        out_descr[i, "n_obs"] <- nrow(x1)
        # start date
        out_descr[i, "start_date"] <- min(x1$date, na.rm = TRUE)
        # end date
        out_descr[i, "end_date"] <- max(x1$date, na.rm = TRUE)
        # ending number of graves
        out_descr[i, "graves_cum"] <- x1[which.max(x1$date), "graves_best"] 
        # starting surface area
        out_descr[i, "area_start"] <- x1[which.min(x1$date), "area"] 
        # ending surface area
        out_descr[i, "area_end"] <- x1[which.max(x1$date), "area"]             
      }
    # Save
      out_descr[, "start_date"] <- as.Date(out_descr[, "start_date"])
      out_descr[, "end_date"] <- as.Date(out_descr[, "end_date"])
      out_descr[, "district"] <- as.factor(out_descr[, "district"])
      write.csv(out_descr, "out_table_description_cemeteries.csv", row.names = FALSE)
    

#.........................................................................................                            
### Comparing estimates with OCHA and Barakaat cemetery committee data, by cemetery

  #...................................   
  ## Prepare data
     
    # Calculate new graves every day, for each cemetery
    obs <- obs[order(obs[, "cemetery"], obs[, "date"]), ]
    for (i in sort(unique(obs$cemetery)) ) {
      
      # get cemetery data
      x1 <- subset(obs, cemetery == i)
      
      # compute daily new graves
      obs[obs$cemetery == i, x2] <- c(0, diff(x1[, "new_graves_best_ipol"]))
      
    }  
    
    # Aggregate by month           
    obs[, "month"] <- month(obs[, "date"])    
    obs[, "year"] <- year(obs[, "date"])    
    obs_m <-  aggregate(obs[, "new_graves_best_ipol"], by = obs[, c("cemetery", "year", "month")], FUN = sum) 
    colnames(obs_m) <- c("cemetery", "year", "month", "new_graves_best_ipol")     

    # OCHA dataset
    ocha[, "month"] <- month(ocha[, "date"])    
    ocha[, "year"] <- year(ocha[, "date"])    
    ocha_m <-  aggregate(ocha[, "graves_ocha"], by = ocha[, c("cemetery", "year", "month")], FUN = sum) 
    colnames(ocha_m) <- c("cemetery", "year", "month", "new_graves_ocha")     
    
    # Merge OCHA into main dataset
    obs_m <- merge(obs_m, ocha_m, by = c("cemetery", "year", "month"), all = TRUE)
    obs_m[obs_m$year == 2020 & obs_m$month > 9, "new_graves_best_ipol"] <- NA
    
    # Combine Barakaat 1 and 2 into one
    x1 <- subset(obs_m, ! cemetery %in% c("Barakaat 1", "Barakaat 2"))
    x2 <- subset(obs_m, cemetery %in% c("Barakaat 1", "Barakaat 2"))
    x3 <- aggregate(x2[, c("new_graves_best_ipol", "new_graves_ocha")], by = x2[, c("year", "month")], FUN = sum, na.rm = TRUE)  
    x3[, "cemetery"] <- "Barakaat 1 + 2"
    x3 <- x3[, colnames(x3)[c(5, 1:4)]]
    obs_m <- rbind(x1, x3)
    obs_m[obs_m$cemetery == "Barakaat 1 + 2" & obs_m$new_graves_ocha == 0, "new_graves_ocha"] <- NA
    obs_m[obs_m$month == 4 & obs_m$year == 2020, "new_graves_ocha"] <- NA
    
    # Merge in Barakaaat committee data
    obs_m <- merge(obs_m, barakaat_committee_m, by = c("cemetery", "year", "month"), all = TRUE)
    
    # Comparison months
    x1 <- intersect(which(is.na(obs_m$new_graves_ocha)), which(is.na(obs_m$new_graves_bdc)))
    x2 <- obs_m[-x1, ]
    x2 <- subset(x2, ! month %in% c(10, 11))
    write.csv(x2, "out_comparison_sources.csv", row.names = FALSE)
    
#.........................................................................................
### ENDS



  