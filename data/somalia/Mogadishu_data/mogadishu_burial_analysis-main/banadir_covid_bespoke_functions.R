#..........................................................................................
### + ESTIMATION OF MORTALITY DURING THE COVID-19 PANDEMIC IN BANADIR, SOMALIA (2020) + ###
#..........................................................................................

#..........................................................................................
## ----------------- R CODE TO PREPARE DATA AND FIT STATISTICAL MODELS ----------------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2021)
                                          # francesco.checchi@lshtm.ac.uk 


#.........................................................................................
### Bespoke functions

  #...................................  
  ## Function to perform K-fold cross-validation and, if desired, plot predictions vs. observations for each fold
    f_cv <- function(f_fit, f_k_folds, f_plot) {
      
      # access dataset used for fit
      if (class(f_fit)[1] == "lm") {f_obs <- fit$data}
      if (typeof(f_fit) == "S4") {f_obs <- model.frame(fit)}
      if (class(f_fit) == "glmmTMB") {f_obs <- model.frame(fit)}
      if ("gamlss" %in% class(f_fit)) {
        f_obs <- model.frame(f_fit) ;
        colnames(f_obs) <- all.vars(formula(f_fit))
        }
     
        # if there is an offset, need to reconstitute it
        x1 <- grep("offset", colnames(f_obs), TRUE)
        if (length(x1) > 0) {
          x2 <- ! all.vars(formula(f_fit)) %in% colnames(f_obs)
          colnames(f_obs)[x1] <- all.vars(formula(f_fit))[x2]
          f_obs[, x1] <- exp(f_obs[, x1])
        }

      # select observations for which all the formula variables are non-missing
      f_obs <- f_obs[complete.cases(f_obs[, all.vars(formula(f_fit)) ] ), ] 

      # determine dependent variable
      f_dep <- all.vars(formula(f_fit))[1]
    
      # determine number of folds if f_k_folds = NA (i.e. LOOCV case)
      if (is.na(f_k_folds) == TRUE) { x1 <- nrow(f_obs) }
      if (is.na(f_k_folds) == FALSE) { x1 <- f_k_folds }

      # is this a least-squares model?
      f_ols <- FALSE
      if (class(f_fit)[1] == "lm") {f_ols <- TRUE}
      if ( typeof(f_fit) == "S4") {f_ols <- ifelse("family" %in% names(f_fit@resp), FALSE, TRUE) }
      
      # shuffle dataset
      f_obs <- f_obs[sample(nrow(f_obs), nrow(f_obs), replace=FALSE), ]
    
      # split data into K folds
        # remove a few rows so as to come up with a n row divisible by K
        f_obs <- f_obs[1:(floor(nrow(f_obs)/x1) * x1), ]
        # split
        folds <- split(f_obs, (0:(nrow(f_obs)-1) %/% (nrow(f_obs)/x1)))
     
      # fit model on all the unfolded sets and track square residuals of model fit on each fold, as well as predictions for each fold  
        # vector to hold squared residuals and other statistics
        errors <- c()
        aic <- c()
        # vector to hold observations and predictions for each fold
        observations <- c()
        predictions <- c()
        
      for (i in 1:length(folds) ) {	
        # progress statement
        print(paste("now on fold ", i, " of ", length(folds), sep = "") )
        # fit on all data but the fold
        data_now <- do.call("rbind", folds[-i])
        
        if (f_ols == FALSE) { cv.fit <- update(f_fit, formula=formula(f_fit),  family=family(f_fit)[[1]], data = data_now) }
        if (f_ols == TRUE) { cv.fit <- update(f_fit, formula=formula(f_fit),  data = data_now) }
        # calculate squared residual of model when predicting fold data, add to other errors for single test observations
        x1 <- predict(cv.fit, newdata = folds[[i]], type = "response", allow.new.levels = TRUE)
        # update output
        observations <- c(observations, folds[[i]][, f_dep])
        predictions <- c(predictions, x1)
        errors <- c(errors , (folds[[i]][, f_dep] -  x1 )^2 )
        aic <- c(aic, AIC(cv.fit))

      }
      
      # return RMSE across all folds
      print("mean RMSE across all folds:")
      print(mean(errors, na.rm = TRUE)^0.5)
      
      # return AIC across all folds
      print("mean AIC across all folds:")
      print(mean(aic, na.rm = TRUE))

      # if plot is desired...
      if (f_plot == TRUE) {
        # prepare data
        x1 <- as.data.frame(cbind(observations, predictions))
        colnames(x1) <- c("observations", "predictions")
          # if on log scale, back-transform to linear
          if (grepl("ln", f_dep) ) {x1[, c("observations", "predictions")] <- exp(x1[, c("observations", "predictions")]) }

        # plot
        plot <- ggplot(x1) +
          geom_point(aes(x = observations, y = predictions), size=2, colour = brewer_pal(palette = "Dark2")(2)[1]) + 
          theme_bw() +
          scale_x_continuous("observed") +
          scale_y_continuous("predicted") +  
          geom_abline(intercept = 0, slope = 1, colour = brewer_pal(palette = "Dark2")(2)[2] ) +
          theme(axis.title = element_text(colour="grey20")) +
         ggtitle(paste("accuracy of predictions on cross-validation; model to predict ", all.vars(formula(f_fit))[1] , sep="") ) +
         theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") )
        
        print("plot shows accuracy of predictions on cross-validation")
        print(plot)
        
        # return plot data
        invisible(x1)
        
      }
          
    }
    
       
     
  #...................................
  ## Function to examine model diagnostics of model fit
    f_diag_fit <- function(f_fit) {
      
      # is the model mixed (i.e. does it include a random effect)?  
      f_lmm <- grepl("|", deparse1(formula(f_fit)) )
      
      # prepare data to plot
      if (f_lmm == TRUE) { x1 <- data.frame(resid(f_fit), fitted(f_fit), sd(resid(f_fit)) ) }
      if (f_lmm == FALSE) { x1 <- data.frame(residuals(f_fit), fitted(f_fit), sd(residuals(f_fit)) ) }
      colnames(x1) <- c("residuals", "fitted_values", "sd_residuals")
      
      # normal distribution of residuals
      plot <- ggplot(x1) +
        geom_histogram(aes(x = residuals), colour = brewer_pal(palette = "Greens")(9)[9], 
          fill = brewer_pal(palette = "Greens")(9)[4] ) +
        theme_bw() +
        scale_x_continuous("residuals" ) +
        ggtitle(paste("distribution of residuals for ", all.vars(formula(f_fit))[1] , sep="")) +
        theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") ) +
        stat_function(fun = dnorm, args = list(mean = mean(x1$residuals), sd = sd(x1$residuals) ) )
      print(plot)
      
      # homoskedasticity
      plot <- ggplot(x1) +
        geom_point(aes(x = fitted_values, y = residuals), size = 2, colour = brewer_pal(palette = "Dark2")(2)[1] ) +
        theme_bw() +
        scale_x_continuous("fitted values" ) +
        scale_y_continuous("residuals" ) +
        geom_abline(intercept = 0, slope = 0, linetype = "dotted", colour = brewer_pal(palette = "Dark2")(2)[2] ) +
        ggtitle(paste("residuals vs. fitted values for ", all.vars(formula(f_fit))[1] , sep="")) +
        theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") )
      print(plot)
      
      # quantile-quantile plot
      plot <- ggplot(x1, aes(sample = residuals / sd_residuals)) +
        theme_bw() +        
        stat_qq(size = 2, colour = brewer_pal(palette = "Dark2")(2)[1]) + 
        stat_qq_line(colour = brewer_pal(palette = "Dark2")(2)[2]) +
        ggtitle(paste("quantile-quantile plot of residuals for ", all.vars(formula(f_fit))[1] , sep="") ) +
        theme(plot.title = element_text(color = "grey20", size = 12, face = "plain") ) +
        scale_x_continuous("theoretical quantiles" ) +
        scale_y_continuous("standardised residual quantiles")
      print(plot)      

    }  
   
    
  #...................................      
  ## Wrapper function to perform diagnostics and display fit results  
    f_do <- function(f_fit, f_diag_fit, f_cv, f_k_folds, f_plot, f_do_cv) {
      # what kind of model is it?
        # is it generalised linear?
        x1 <- TRUE
        if (class(f_fit) == "lm") {x1 <- FALSE}
        if ( typeof(f_fit) == "S4") {x1 <- ifelse("family" %in% names(f_fit@resp), TRUE, FALSE) }
        if ( class(f_fit) == "glmmTMB") {x1 <- ifelse(as.character((f_fit$modelInfo)$family)[1] == "gaussian", FALSE, TRUE) }
        # is it mixed (i.e. does it include a random effect)?  
        x2 <- grepl("|", deparse1(formula(f_fit)) )
      
      # implement and print output 
        print(paste("#.......Output of ", ifelse(x1, "generalised ", ""),
          "linear ", ifelse(x2, "mixed ", ""), "model:" , ".............", sep = "") )
        # summary
        print( summary(f_fit) )
        
        # tidy summary (fixed effects only, exponentiated if generalised linear model)
        print( tidy(f_fit, conf.int = TRUE, exponentiate = x1, effects = "fixed") )
        print("/")
        
        # model diagnostics if least-squares model
        if (x1 == FALSE) {
          print("model diagnostics:")
          f_diag_fit(f_fit)
          print("/")
        }
        
        # cross-validation if desired
        if (f_do_cv == TRUE) { 
          # cross-validation
          print("cross-validation:")
          out <- f_cv(f_fit, f_k_folds, f_plot)
          return(out)
          print("/")
        }

    }

    
  
  # #...................................
  # ## Function to fit a general linear model and display clean results
  #   f_glm <- function(f_dep, f_preds, f_offset, f_data, f_family) {
  #     # write the model formula
  #     form <- as.formula( paste(f_dep, " ~ ", paste(f_preds, collapse= " + "), 
  #     ifelse(is.na(f_offset), "", paste(" + offset(log(", f_offset, ") )", sep = "") ), sep="")  )
  # 
  #     # select observations for which all the formula variables are non-missing or non-NA
  #     f_obs <- f_data[complete.cases(f_data[, all.vars(form) ] ), ]
  # 
  #     # fit model depending on distributional assumption
  #     if (f_family == "nb") { fit <- tryCatch(glm.nb(form, data = f_obs), error=function(er) {return(FALSE)} ); 
  #       if (fit == FALSE) {fit <- glmmTMB(form, data = f_obs, family=nbinom2) } 
  #     }
  #     if (f_family != "nb") {fit <- glm(form, data = f_obs, family = f_family) }
  # 
  #     # return fit
  #       return(fit)
  #   }  
  # 

  #...................................
  ## Function to fit a mixed general linear model and display clean results (using glmmMTB package)
    f_glmm <- function(f_dep, f_preds, f_reff, f_offset, f_data, f_family) {
      # write the model formula
        form <- as.formula( paste(f_dep, " ~ ", paste(f_preds, collapse= " + "),
          ifelse(is.na(f_offset), "", paste(" + offset(log(", f_offset, ") )", sep = "") ), " + ", f_reff, sep="")  )

      # select observations for which all the formula variables are non-missing or non-NA
      f_obs <- f_data[complete.cases(f_data[, all.vars(form) ] ), ]

      # fit model
      fit <- glmmTMB(form, data = f_obs, family = f_family)

      # return fit
        return(fit)
    }


  #...................................   
  ## Function to plot histograms of model variables
  f_hist <- function(f_var, f_data, f_lims) {
      
    plot <- ggplot(f_data)
        
      # if the variable has >= 20 unique values...
        if (length(unique(na.omit(f_data[, f_var]))) >= 20) {
          plot <- plot + geom_histogram(aes(x = as.numeric(f_data[, f_var]) ), color="seagreen", fill="seagreen3" ) +
                  theme_bw() + xlab(f_var) + scale_x_continuous(expand = c(0, 0), limits = f_lims )
        }
   
      # otherwise...
        if (length(unique(na.omit(f_data[, f_var]))) < 20) {
          plot <- plot + geom_histogram(aes(x = as.numeric(f_data[, f_var]) ), stat="count", color="seagreen", fill="seagreen3") +
                  theme_bw() + xlab(f_var) + scale_x_continuous(expand = c(0, 0), limits = f_lims )
        }
          
      print(plot)
    }

  
  #...................................   
  ## Function to compute prediction confidence intervals or prediction profiles for bootstrapping, from a GAMLSS fit
    # by posterior simulation, as in https://r.789695.n4.nabble.com/Prediction-interval-with-GAM-td3460175.html 
  f_interval <- function(f_fit, f_n_boot, f_obs, f_profile) {
    
    # Extract coefficient estimates and variance-covariance matrix from fit
    beta <- na.omit( coef(f_fit) ) # only mu coefficients, i.e. for estimating linear predictor; omit NAs as random effects come out as NA
      # convert to inverse (for later operation)
      beta_inv <- ginv(matrix(beta))
    vcov_matrix <- vcov(f_fit)[1:length(beta), 1:length(beta)] # only mu coefficients
      # Cholesky decomposition
      cholesky_dec <- chol(vcov_matrix)
    
    # Simulate a large number of random beta coefficient values drawn from an assumed normal error distribution
      # around their posterior estimates
    beta_sim <- t(cholesky_dec) %*% matrix(rnorm(f_n_boot * length(beta)), length(beta), f_n_boot) + as.vector(beta) 
    
    # Predict both mu and sigma on the new data
    mu_pred <- predict(f_fit, newdata = f_obs[, all.vars(formula(f_fit))])
    sigma_pred <- predict(f_fit, newdata = f_obs[, all.vars(formula(f_fit))], what = "sigma") # for later
    
    # Produce a 'linear prediction matrix'
      # since lp_matrix %*% coef(f_fit) = mu_pred , the code below uses the inverse (i.e. matrix 'division')
      # to get the lp_matrix, which gamlss doesn't output
    lp_matrix <- mu_pred %*% beta_inv 
    
    # Compute point estimate linear predictions and back-transform them if appropriate, for all the random beta values
      # only exponential back-transform implemented here for now
    pred_sim <- lp_matrix %*% beta_sim
    if ( f_fit$mu.link == "log" ) {pred_sim <- exp(pred_sim)}
   
    # Lastly, generate random predictions by combining the point estimates with the other estimated distributional parameters
      # only  basic 'count' distributions implemented here so far
    if (family(f_fit) == "PO") 
      { rand_sim <- matrix( rPO(n = prod(dim(pred_sim)), mu = pred_sim), 
        nrow(pred_sim), ncol(pred_sim) ) }
    if (family(f_fit) == "NBI") 
      { rand_sim <- matrix( rNBI(n = prod(dim(pred_sim)), mu = pred_sim, sigma = rep(exp(sigma_pred), f_n_boot) ), 
        nrow(pred_sim), ncol(pred_sim) ) }
    if (family(f_fit) == "NBII") 
      { rand_sim <- matrix( rNBII(n = prod(dim(pred_sim)), mu = pred_sim, sigma = rep(exp(sigma_pred), f_n_boot) ), 
        nrow(pred_sim), ncol(pred_sim) ) }
    
    # If desire 95% confidence interval...
    if (f_profile == FALSE) {
      x1 <- mu_pred
      if ( f_fit$mu.link == "log") { x1 <- exp(mu_pred) }
      out <- cbind(f_obs,
                   x1,
                   t(apply(rand_sim, 1, quantile, prob = c(0.025,0.975), na.rm = TRUE) )
                   )
      out <- as.data.frame(out)
      colnames(out) <- c(colnames(f_obs), "predicted", "predicted_lci", "predicted_uci")
      out[, c("predicted", "predicted_lci", "predicted_uci")] <- round(out[, c("predicted", "predicted_lci", "predicted_uci")], digits = 0)
      return(out)
    }
    
    # ...else output entire prediction profile, sorted ascendingly
    if (f_profile == TRUE) {
      out <- apply(rand_sim, 1, sort)
      return(out)           
    }
    
  }

   
  #...................................  
  ## Function to plot and output model predictions vs. observations for each cemetery, overall and for a specific date range
    f_plot_gof <- function(f_fit, f_data, f_date_range, f_pred_smooth, f_spar) {
      
      # prepare data
        # actual dataset
        x1 <- unique(c(all.vars(formula(f_fit)), "graves", "cemetery", "period_covid", "date") )
        x2 <- f_data[complete.cases(f_data[, x1[x1 != "graves"] ]),  x1]

        # counterfactual dataset
        x4 <- grep("time_covid|time_piece2", colnames(x2), value = TRUE)
        x5 <- min(subset(f_data, period_covid == "epidemic")[, "date"], na.rm = TRUE)
        x5 <- f_data[f_data$date == x5, x4]
        x5 <- unique(x5)
        x3 <- x2
        x3[, x4] <- x5
               
      # predict
      x2[, "predictions"] <- predict(f_fit, newdata = x2[, all.vars(formula(f_fit))], type = "response")
      x2[, "counterfactuals"] <- predict(f_fit, newdata = x3[, all.vars(formula(f_fit))], type = "response")
        # back-transform if needed
        x3 <- colnames(model.frame(f_fit))[1]
        if (grepl("log", x3) | grepl("ln", x3) ) {
            x2[, "predictions"] <- exp(x2[, "predictions"]) ;
            x2[, "counterfactuals"] <- exp(x2[, "counterfactuals"])
        }
      
      # smooth predictions if desired
      if (f_pred_smooth == TRUE) {
        
        x2 <- x2[order(x2[, "cemetery"], x2[, "date"]), ]
       
        # for each cemetery...
        for (i in sort(unique(x2$cemetery))) {
          # observations from the cemetery...
          x4 <- grep("time_base", colnames(x2), value = TRUE)
          x3 <- subset(x2, cemetery == i)[, c(x4, "predictions", "counterfactuals")]

          # loess
          x2[x2$cemetery == i, "predictions"] <- predict(loess(x3[, "predictions"] ~ x3[, x4], span = f_spar) , x3[, x4] )
          x2[x2$cemetery == i, "counterfactuals"] <- predict(loess(x3[, "counterfactuals"] ~ x3[, x4], span = f_spar) , x3[, x4] )

          # # spline smooth
          # x2[x2$cemetery == i, "predictions"] <- predict(smooth.spline(x3[, x4], x3[, "predictions"], spar = f_spar) , x3[, x4] )$y
          # x2[x2$cemetery == i, "counterfactuals"] <- predict(smooth.spline(x3[, x4], x3[, "counterfactuals"], spar = f_spar) , x3[, x4] )$y
        }  
      }  
          
      # plot
        plot <- ggplot(x2) +
           geom_line(aes(x = date, y = predictions ), linetype = "dashed", size = 0.7,  
             colour = brewer_pal(palette = "Dark2")(2)[2] ) +
           geom_point(aes(x = date, y = graves, colour = period_covid), size = 1.5) +
           geom_line(aes(x = date, y = predictions, colour = period_covid), size = 0.9 ) +
           geom_line(aes(x = date, y = counterfactuals), linetype = "dashed", 
             colour = brewer_pal(palette = "Dark2")(2)[1], size = 0.7 ) +
           scale_colour_manual(values = brewer_pal(palette = "Dark2")(2)) +
           scale_y_continuous("new graves since start of analysis period") +
           theme_bw() +
           facet_wrap(~cemetery, nrow=5, scales = "free_y") +
           guides(fill = FALSE) +
           theme(legend.position="bottom", legend.direction="horizontal") +
           scale_x_date("", minor_breaks=NULL, date_breaks="4 months", date_labels = "%b-%Y" ) +
           labs(colour = "Period:  ") +
           ggtitle("within-sample predictions (lines) versus observations (dots), by cemetery" ) +
           theme(legend.title = element_text(color="grey20", size=11),
                 strip.text.x = element_text(color="grey20", size=11),
                 legend.text = element_text(color="grey20", size=11),
                 axis.title.x = element_text(color="grey20", size=11), 
                 axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
                 axis.line.y = element_line(color = "grey20"),
                 axis.ticks.y = element_line(color = "grey20"),
                 axis.text.y = element_text(color = "grey20", size=11),
                 axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
                 plot.title = element_text(color = "grey20", size = 12, face = "plain", hjust = 0.5), 
                 plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
                 )
      print(plot)
      
      # now zoom in on desired region of plot
      x3 <- seq(f_date_range[1], f_date_range[2], by = 1)
      plot <- plot %+% subset(x2, date %in% x3 ) + scale_x_date("", date_breaks="1 month", date_labels = "%b-%Y" )
      print(plot)
      
      # output predictions
      invisible(x2)
      
    }
    
    
  #...................................      
  ## Function to plot smoothed fit and data for each single cemetery
    f_plot_single <- function(f_fit, f_data) {

      # prepare data
      x1 <- f_data[, unique(c(all.vars(formula(f_fit)), "graves", "cemetery", "period_covid", "date") )]
      x2 <- as.data.frame( x1[, all.vars(formula(f_fit))] )
      colnames(x2) <- all.vars(formula(f_fit))
      x1[, "predictions"] <- predict(f_fit, newdata = x2, type = "response")
        # back-transform if needed
        x2 <- colnames(model.frame(f_fit))[1]
        if (grepl("log", x2) | grepl("ln", x2) ) {
          x1[, "predictions"] <- exp(x1[, "predictions"])
        }

      # plot
        plot <- ggplot(x1) +
           geom_line(aes(x = date, y = predictions ), linetype = "dashed", size = 0.7,  
             colour = brewer_pal(palette = "Dark2")(2)[2] ) +
           geom_point(aes(x = date, y = graves, colour = period_covid), size = 1.5) +
           geom_line(aes(x = date, y = predictions, colour = period_covid), size = 0.9 ) +
           scale_colour_manual(values = brewer_pal(palette = "Dark2")(2)) +
           scale_y_continuous("number of graves") +
           theme_bw() +
           guides(fill = FALSE) +
           theme(legend.position="bottom", legend.direction="horizontal") +
           scale_x_date("", minor_breaks=NULL, date_breaks="4 months", date_labels = "%b-%Y" ) +
           labs(colour = "Period:  ") +
           ggtitle(paste("within-sample predictions (lines) versus observations (dots)", "\n",
             "Name of cemetery: ", unique(x1$cemetery), sep = "") ) +
           theme(legend.title = element_text(color="grey20", size=11),
                 legend.text = element_text(color="grey20", size=11),
                 axis.title.x = element_text(color="grey20", size=11), 
                 axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
                 axis.line.y = element_line(color = "grey20"),
                 axis.ticks.y = element_line(color = "grey20"),
                 axis.text.y = element_text(color = "grey20", size=11),
                 axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
                 plot.title = element_text(color = "grey20", size = 12, face = "plain", hjust = 0.5), 
                 plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
                 )
      print(plot)
      
      return(x1)    
    }    
  
 
  #...................................      
  ## Function to plot linear interpolation or smooth spline and observations for each single cemetery
    f_plot_ipol <- function(f_ipol, f_data) {

      # prepare data
      f_ipol <- data.frame(f_ipol)
      colnames(f_ipol) <- c("time_base", "predictions")
      x1 <- merge(f_ipol, f_data, by = "time_base")
      colnames(x1)[grep("graves|area", colnames(x1) )] <- "observations"

      # plot
        plot <- ggplot(x1) +
           geom_line(aes(x = date, y = predictions ), linetype = "dashed", size = 0.7,  
             colour = brewer_pal(palette = "Dark2")(2)[2] ) +
           geom_point(aes(x = date, y = observations, colour = period_covid), size = 1.5) +
           geom_line(aes(x = date, y = predictions, colour = period_covid), size = 0.9 ) +
           scale_colour_manual(values = brewer_pal(palette = "Dark2")(2)) +
           scale_y_continuous("number of graves or surface area") +
           theme_bw() +
           guides(fill = FALSE) +
           theme(legend.position="bottom", legend.direction="horizontal") +
           scale_x_date("", minor_breaks=NULL, date_breaks="4 months", date_labels = "%b-%Y" ) +
           labs(colour = "Period:  ") +
           ggtitle(paste("Smoothed interpolated values (lines) versus observations (dots)", "\n",
             "Name of cemetery: ", unique(x1$cemetery), sep = "") ) +
           theme(legend.title = element_text(color="grey20", size=11),
                 legend.text = element_text(color="grey20", size=11),
                 axis.title.x = element_text(color="grey20", size=11), 
                 axis.text.x = element_text(color = "grey20", size=10, angle=315, hjust=0, vjust=0),               
                 axis.line.y = element_line(color = "grey20"),
                 axis.ticks.y = element_line(color = "grey20"),
                 axis.text.y = element_text(color = "grey20", size=11),
                 axis.title.y = element_text(color="grey20", margin = margin(r = 10), size=11 ),
                 plot.title = element_text(color = "grey20", size = 12, face = "plain", hjust = 0.5), 
                 plot.margin = unit(c(0.5,2,0.5,0.5), "cm")
                 )
      print(plot)
      
      return(x1)    
    }    
      

        
