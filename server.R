#-------------------------------------------------------------------------
#  This file is part of the TVAR forecasting application https://github.com/maximegodin/tvar_forecasting
#
#  This application is governed by the MIT license. 
#  You can  use, modify and/ or redistribute this code under the terms
#  of the MIT license:  https://opensource.org/licenses/MIT
#
#  Maxime Godin, Amina Bouchafaa, Ecole Polytechnique
#  March 29th, 2018
#-------------------------------------------------------------------------

library(shiny)
library(ggplot2)
library(gridExtra)

shinyServer(function(input, output){
  
  output$sample_id_slider <- renderUI(
    {
      input$simulate
      sliderInput("sample_id","Sample number",min=1,max=isolate(input$n_obs),ticks=FALSE,value=1)
    }
  )
  
  output$window_slider <- renderUI(
    {
      input$simulate
      sliderInput("plot_window","Subsample window",min=1,
                  max=isolate(2**input$log_Tmax),
                  ticks=FALSE,value=c(1,isolate(2**input$log_Tmax)))
    }
  )
  
  plots <- reactiveValues()
  plots$settings <- reactiveValues()
  plots$yw <- reactiveValues()
  plots$nlms <- reactiveValues()
  plots$expo_agg <- reactiveValues()
  
  # TVAR settings tab
  
  # random numbers to use in the TVAR coefficients computation
  rd_seeds <- eventReactive(input$roll_dice,
                            {matrix(runif(input$richness*input$order,min=-1,max=1),ncol=input$order)},ignoreNULL = FALSE)
  
  # function computing the TVAR coefficients
  reg_coeffs <- reactive({function(u){compute_reg_coeffs(u,rd_seeds(),isolate(input$delta))}})
  
  
  
  # plot of the TVAR coefficients
  output$coeffs_plot <- renderPlot({ 
    p <- plot_coefficients(reg_coeffs(),isolate(input$order))
    plots$settings$coeffs <- p
    p
  })
  
  # sample used in prediction
  tvar_sample <- eventReactive(input$simulate,{TVAR_gen(2^(input$log_Tmax),input$n_obs,input$sigma, input$white_noise_gen,rd_seeds(),input$delta)},ignoreNULL = FALSE)
  # simulation notifications
  observeEvent(input$simulate,showNotification("Simulation in progress",type="warning",duration = 5),ignoreNULL = FALSE)
  observeEvent(tvar_sample(),{req(tvar_sample())
    showNotification("Simulation finished",type="message")})
  
  #plot the TVAR sample
  output$sample_plot <- renderPlot({
    req(tvar_sample())
    p <- plot_sample(tvar_sample())
    plots$settings$sample <- p
    p
  })
  
  #plot the subsample
  output$subsample_plot <- renderPlot({
    req(input$plot_window)
    p <- plot_sample(tvar_sample(),input$sample_id) + facet_zoom(x = (x >=input$plot_window[1]) & (x <=input$plot_window[2]))
    plots$settings$subsample <- p
    p
  })
  
  #plot correlations
  output$acf_pacf_plot <- renderPlot({
    req(input$plot_window)
    p1 <- ggAcf(tvar_sample()[(input$plot_window[1]:input$plot_window[2])+1,input$sample_id], demean=FALSE)+
      ggtitle("")+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold"))
      
    p2 <- ggPacf(tvar_sample()[(input$plot_window[1]:input$plot_window[2])+1,input$sample_id], demean=FALSE)+
      ggtitle("")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    p <- grid.arrange(p1,p2)
    plots$settings$acf_pacf <- p
    p
  })
  
  # TVAR prediction
  
  # sliders for bandwidth              
  output$log_bandwidth_slider <- renderUI(
    {coeffs_by_M()
      isolate({sliderInput(inputId = "log_M", label = "log bandwidth", 
                           min = 0, max = input$log_Tmax+1, 
                           value = 0, step = 1, ticks=FALSE)})})
  
  
  # tapering function and plot  
  taper_pred <- reactive({function(u) {compute_taper(u,input$taper_pred,TRUE)}})
  output$taper_plot_pred <- renderPlot({plot_taper(taper_pred())},bg="transparent")

  # estimation of the coefficients for the different bandwidths
  coeffs_by_M <- eventReactive(input$predict_yw,{
    req(tvar_sample())
    showNotification("Prediction in progress",type="warning", duration=5)
    estimate_coeffs_yule_walker_by_M(tvar_sample(),
                                     isolate(2**(input$log_Tmax)),
                                     isolate(input$order),
                                     isolate(input$log_Tmax+1),
                                     input$taper_pred, TRUE)})
  
  
  # prediction for the diffrent bandwidths
  predictions_by_M <- eventReactive(coeffs_by_M(),{
    predict_TVAR_by_M(tvar_sample(), 
                      2**(input$log_Tmax),
                      input$order,
                     input$log_Tmax+1,
                      coeffs_by_M())})
  
  # notification for the end of prediction
  observeEvent(predictions_by_M(), {req(predictions_by_M())
    showNotification("Prediction finished",type="message")})
  
  # best predictions computed from the real TVAR coefficients
  best_predictions_yw <- reactive({ req(tvar_sample())
    apply(tvar_sample(),MARGIN = 2,
                                      function(x) best_predict_TVAR(x, isolate(2**(input$log_Tmax)),isolate(rd_seeds()),isolate(input$delta)))})
  
  # sliders for the bandwidth range in Romberg bias reduction
  output$log_bandwidth_slider_bandwidth_romberg <- renderUI({
    input$predict_yw
    sliderInput(inputId = "log_M_romberg", label = "log Romberg bandwidth", 
                        min = 0, max = isolate(input$log_Tmax+1), 
                        value = isolate(c(0,input$log_Tmax+1)), step = 1, ticks=FALSE)
  })
  

  # Romberg predictions for the specified bandwidth range
  yw_romberg_predictions <- reactive({req(input$log_M_romberg,predictions_by_M())
    predict_TVAR_romberg(isolate(input$n_obs),
                         predictions_by_M(),
                         input$log_M_romberg[1],
                         input$log_M_romberg[2])})
  
  # Romberg estimations for the specified bandwidth range
  yw_romberg_estimations <- reactive({req(input$log_M_romberg,coeffs_by_M())
    estim_coeffs_romberg(isolate(input$n_obs),
                         coeffs_by_M(),
                         input$log_M_romberg[1],
                         input$log_M_romberg[2],
                         isolate(input$order))})
  
  # Plots
  
  #Estimation MSE
  output$estim_mse_yw_plot <- renderPlot({req(coeffs_by_M())
    p <-
    if(input$romberg_yw){
      
      req(yw_romberg_estimations())
      mse_romberg <- mean(colSums((yw_romberg_estimations()-isolate({reg_coeffs()(1)}))**2))
      plot_estim_mse_yw(coeffs_by_M(),
                                isolate({reg_coeffs()(1)}), input$log_M, mse_romberg)
    }
    else{
      plot_estim_mse_yw(coeffs_by_M(),
                                isolate({reg_coeffs()(1)}),input$log_M)
    }
    plots$yw$estim_mse <- p
    p
  })
  
  # Prediction MSE
  output$pred_mse_yw_plot <- renderPlot({req(predictions_by_M(),best_predictions_yw(),input$log_M)
    p <- if (input$romberg_yw){
      req(yw_romberg_predictions())
      mse_romberg <- mean((yw_romberg_predictions()-best_predictions_yw())**2)
      plot_pred_mse_yw(predictions_by_M(),best_predictions_yw(),input$log_M, mse_romberg)
    }else{
      plot_pred_mse_yw(predictions_by_M(),best_predictions_yw(),input$log_M)
    }
    plots$yw$pred_mse <- p
    p
  })
  
  # Prediction diagnostic plot
  output$pred_diag_yw_plot <- renderPlot({req(predictions_by_M(),best_predictions_yw(),input$log_M)
    p <- plot_pred_diag_yw(isolate(2**(input$log_Tmax)),tvar_sample(),predictions_by_M()[input$log_M+1,],best_predictions_yw())
    plots$yw$pred_diag <- p
    p})
  
  # Romberg prediction diagnostic plot
  output$romb_pred_diag_yw_plot <- renderPlot({
    p <- plot_romb_pred_diag_yw(isolate(2**(input$log_Tmax)),tvar_sample(),predictions_by_M()[input$log_M+1,],yw_romberg_predictions(),best_predictions_yw())
    plots$yw$romb_pred_diag <- p
    p})
  
  # NLMS tab
  
  # slider for log mu
  output$nlms_log_mu_input <- renderUI({
    sliderInput("log_mu","log mu",min=input$log_mu_min,max=input$log_mu_max,ticks=FALSE, value=input$log_mu_min, step=.2,animate=animationOptions(interval=250))
  })
  
  # slider for log mu
  output$nlms_romb_log_mu_input <- renderUI({
    sliderInput("romb_log_mu","log mu (Romberg)",min=input$log_mu_min,max=input$log_mu_max,ticks=FALSE, value=input$log_mu_min, step=.2,animate=animationOptions(interval=250))
  })
  
  # slider for log mu plot window
  output$nlms_plot_window_input <- renderUI({
    sliderInput("mu_plot_lims","log mu plot window", min=input$log_mu_min,max=input$log_mu_max,ticks=FALSE, value=c(input$log_mu_min,input$log_mu_max), step = .2)
  })
  
  # NLMS estimations
  nlms_estimations <- reactive({req(input$log_mu)
    input$nlms_flush
    predict_NLMS(tvar_sample(), 10**(input$log_mu), isolate(input$order))})
  
  # NLMS romberg estimations
  nlms_romberg_estimations <- reactive({
    req(input$romb_log_mu,input$gamma)
    (predict_NLMS(tvar_sample(), 10**(input$romb_log_mu), isolate(input$order))-input$gamma*predict_NLMS(tvar_sample(), input$gamma*10**(input$romb_log_mu), isolate(input$order)))/(1-input$gamma)})
  
  # NLMS predictions
  nlms_predictions <- eventReactive(nlms_estimations(),
                                    predict_mult_TVAR(tvar_sample(), 2**input$log_Tmax, nlms_estimations()))
  # NLMS romberg predictions
  nlms_romberg_predictions <- eventReactive(nlms_romberg_estimations(),
                                            predict_mult_TVAR(tvar_sample(), 2**input$log_Tmax, nlms_romberg_estimations()))
  # NLMS best predictions
  best_predictions_nlms <- eventReactive(input$simulate,
                                         {apply(tvar_sample(),MARGIN = 2,
                                                function(x) best_predict_TVAR(x, 2**input$log_Tmax,rd_seeds(),input$delta))},ignoreNULL = FALSE)
  # data for the plots and data updates
  nlms_data <- reactiveValues()
  
  observeEvent({input$simulate
                input$nlms_flush},
               {nlms_data$mse_by_mu <- data.frame(mu=numeric(),
                                                  estim_MSE=numeric(),
                                                  pred_MSE = numeric())
                nlms_data$mse_by_mu_gamma <- data.frame(mu = numeric(),
                                                        gamma=numeric(),
                                                        estim_MSE=numeric(),
                                                        pred_MSE=numeric())}, 
               ignoreNULL = FALSE)
                 

  observeEvent(nlms_predictions(),
               { idx <- which(nlms_data$mse_by_mu[,1]==10**(input$log_mu))
                if (length(idx)>0){
                  nlms_data$mse_by_mu <- nlms_data$mse_by_mu[-idx,]
                } 
               estim_err <- nlms_estimations()-reg_coeffs()(1)
               estim_mse <- mean(colSums(estim_err**2))
               pred_mse <- mean((nlms_predictions()-best_predictions_nlms())**2)
               df <- data.frame(mu=10**(input$log_mu),estim_MSE=estim_mse,pred_MSE=pred_mse)
               nlms_data$mse_by_mu <- rbind(nlms_data$mse_by_mu,df)}
  )
  
  observeEvent(nlms_romberg_predictions(),
               {idx <- which(nlms_data$mse_by_mu_gamma[,1]==10**(input$romb_log_mu) & (nlms_data$mse_by_mu_gamma[,2]==input$gamma ))
               if (length(idx)>0){
                 nlms_data$mse_by_mu_gamma <- nlms_data$mse_by_mu_gamma[-idx,]
               } 
               estim_err <- nlms_romberg_estimations()-reg_coeffs()(1)
               estim_mse <- mean(colSums(estim_err**2))
               pred_mse <- mean((nlms_romberg_predictions()-best_predictions_nlms())**2)
               df <- data.frame(mu=10**(input$romb_log_mu),gamma=input$gamma,estim_MSE=estim_mse,pred_MSE=pred_mse)
               nlms_data$mse_by_mu_gamma <- rbind(nlms_data$mse_by_mu_gamma,df)}
               
  )
  
  # Plots
  
  # Estimation MSE
  output$estim_mse_nlms_plot <- renderPlot({
    req(input$mu_plot_lims)
    y <- nlms_data$mse_by_mu %>% filter(mu>=10**input$mu_plot_lims[1],mu<=10**input$mu_plot_lims[2]) %>% select(estim_MSE)
    ymin = min(y)
    ymax = max(y)
    p<-ggplot(nlms_data$mse_by_mu)+
      aes(x=mu,y=estim_MSE)+
      geom_point()+
      geom_line()+
      geom_point(data=nlms_data$mse_by_mu %>% filter(mu==10**(input$log_mu)),size=3)+
      scale_x_log10()+
      scale_y_log10()+
      coord_cartesian(xlim=10**input$mu_plot_lims,ylim = c(ymin,ymax))+
      xlab("mu") + ylab("Coefficients estimation MSE")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    plots$nlms$estim_mse <-p
    p
  })
  
  # Prediction MSE
  output$pred_mse_nlms_plot <- renderPlot({
    req(input$mu_plot_lims)
    y <- nlms_data$mse_by_mu %>% filter(mu>=10**input$mu_plot_lims[1],mu<=10**input$mu_plot_lims[2]) %>% select(pred_MSE)
    ymin = min(y)
    ymax = max(y)
    p <- ggplot(nlms_data$mse_by_mu)+
      geom_point(aes(x=mu,y=pred_MSE))+
      geom_line(aes(x=mu,y=pred_MSE))+
      geom_point(aes(x=mu,y=pred_MSE),data=nlms_data$mse_by_mu %>% filter(mu==10**(input$log_mu)),size=3)+
      scale_x_log10()+
      scale_y_log10()+
      coord_cartesian(xlim=10**input$mu_plot_lims,ylim = c(ymin,ymax))+
      xlab("mu") + ylab("Predictor estimation MSE")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    plots$nlms$pred_mse <-p
    p
  })
  
  # Romberg estimation MSE
  output$romb_estim_mse_nlms_plot <- renderPlot({
    req(input$mu_plot_lims)
    dat0 <- data.frame(nlms_data$mse_by_mu,gamma=rep(0,nrow(nlms_data$mse_by_mu)))
    dat_gamma <- nlms_data$mse_by_mu_gamma %>% filter(gamma == input$gamma)
    dat <- rbind(dat0,dat_gamma)
    y <- dat %>% filter(mu>=10**input$mu_plot_lims[1],mu<=10**input$mu_plot_lims[2]) %>% select(estim_MSE)
    ymin = min(y)
    ymax = max(y)
    foo <- data.frame(gamma = as.factor(c(0,input$gamma)),x = 10**c(input$log_mu,input$romb_log_mu))
    p <- ggplot(dat)+
      aes(x=mu,y=estim_MSE, color=as.factor(gamma))+
      geom_point()+
      geom_line()+
      geom_vline(aes(xintercept=x,color=gamma),data=foo)+
      scale_x_log10()+
      scale_y_log10()+
      coord_cartesian(xlim=10**input$mu_plot_lims,ylim = c(ymin,ymax))+
      xlab("mu") + ylab("Coefficients estimation MSE")+labs(color = "gamma")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            legend.title=element_text(size=12),
            legend.text=element_text(size=12))
    plots$nlms$romb_estim_mse <-p
    p
  })
  
  # Romberg prediction MSE
  output$romb_pred_mse_nlms_plot <- renderPlot({
    req(input$mu_plot_lims)
    foo <- data.frame(method = c("normal","romberg"),x = 10**c(input$log_mu,input$romb_log_mu))
    dat0 <- data.frame(nlms_data$mse_by_mu,method=rep("normal",nrow(nlms_data$mse_by_mu)))
    dat_gamma <- nlms_data$mse_by_mu_gamma %>% filter(gamma == input$gamma) %>% select(-gamma)
    dat_gamma <- data.frame(dat_gamma,method=rep("romberg",nrow(dat_gamma)))
    dat <- rbind(dat0,dat_gamma)
    y <- dat %>% filter(mu>=10**input$mu_plot_lims[1],mu<=10**input$mu_plot_lims[2]) %>% select(pred_MSE)
    ymin = min(y)
    ymax = max(y)
    p <- ggplot(dat)+
      aes(x=mu,y=pred_MSE, color=method)+
      geom_point()+
      geom_line()+
      geom_vline(aes(xintercept=x,color=method),data=foo)+
      scale_x_log10()+
      scale_y_log10()+
      coord_cartesian(xlim=10**input$mu_plot_lims,ylim = c(ymin,ymax))+
      xlab("mu") + ylab("Predictor estimation MSE")+labs(color = "method")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            legend.title=element_text(size=12),
            legend.text=element_text(size=12),
            legend.position="bottom")
    plots$nlms$romb_pred_mse <- p
    p
  })
  
  # Prediction diagnostic plot
  output$pred_diag_nlms_plot <- renderPlot({
    df <- data.frame(pred = nlms_predictions(), best_pred = best_predictions_nlms())
    
   p <- ggplot(df)+
      geom_point(aes(x=pred,y=best_pred))+
      ylab("Best predictions") + xlab("Predicted values")+
      geom_abline(slope=1,intercept=0,aes(linetype="dashed"))+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
   plots$nlms$pred_diag <-p
   p
  })
  
  # Romberg prediction diagnostic plot
  output$romb_pred_diag_nlms_plot <- renderPlot({
    df <- data.frame(normal = nlms_predictions(), romberg = nlms_romberg_predictions(), best_pred = best_predictions_nlms())
    df <- melt(df, id.vars="best_pred",variable.name="method",value.name="pred")
    p <- ggplot(df)+
      geom_point(aes(x=pred,y=best_pred,colour=method))+
      ylab("Best predictions") + xlab("Predicted values")+
      geom_abline(slope=1,intercept=0,linetype="dashed")+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            legend.title=element_text(size=12),
            legend.text=element_text(size=12),
            legend.position="bottom")
    plots$nlms$romb_pred_diag <-p
    p
  })
  
  # NLMS expo agg tab
  
  # output$expo_agg_pred_select <- renderUI({
  #   selectInput("pred_select_expo_agg","Predictions to plot",
  #               choices=c("aggregated",letters[1:N_expo_agg()]), 
  #               selected="aggregated", 
  #               multiple=TRUE, selectize = TRUE)
  # })
  # 
  # N_expo_agg <- reactive({
  #   req(tvar_sample())
  #   ceiling(log(nrow(tvar_sample())))
  # })
  # 
  # eta_expo_agg <- reactive({
  #   req(tvar_sample())
  #   Tmax <- nrow(tvar_sample())
  #   p <- ceiling(4*input$beta_expo_agg+2)
  #   ((log(ceiling(log(Tmax)))/Tmax)**(2/p))/isolate(input$sigma)**2
  # })
  # 
  # output$expo_agg_params_table <- renderTable({
  #   df <- data.frame(c("N","eta"),c(as.integer(N_expo_agg()),eta_expo_agg()))
  #   colnames(df) <- c("parameter","value")
  #   df
  # })
  # 
  # expo_agg_mu <- reactive({
  #   req(tvar_sample())
  #   beta = (((1:N_expo_agg())-1)*input$beta_expo_agg)/N_expo_agg()
  #   (nrow(tvar_sample())-1)**(-2*beta/(2*beta+1))
  # })
  # 
  # output$expo_agg_mu_table <- renderTable({
  #   req(expo_agg_mu())
  #   data.frame(predictor = letters[1:N_expo_agg()], mu = expo_agg_mu(), `log mu` = log10(expo_agg_mu()))
  # })
  # 
  # expo_agg_res <- reactive({
  #   req(tvar_sample(),expo_agg_mu())
  #   predict_exp_agg_NLMS(tvar_sample(),eta_expo_agg(),expo_agg_mu(),isolate(input$order))
  # })
  # 
  # expo_agg_predictions <- reactive({
  #         expo_agg_res()[[1]]})
  # 
  # 
  # best_predictions_expo_agg <- eventReactive(input$simulate,
  #                                            {apply(tvar_sample(),MARGIN = 2,
  #                                                   function(x) best_predict_TVAR(x, 2**input$log_Tmax,rd_seeds(),input$delta))},ignoreNULL = FALSE)
  # 
  # output$pred_mse_expo_agg <- renderPlot({
  #   req(expo_agg_predictions(),best_predictions_expo_agg())
  #   err <- t((t(expo_agg_predictions())-best_predictions_expo_agg())**2)
  #   dat <- data.frame(predictor = factor(c(letters[1:N_expo_agg()],"aggregated"),
  #                                        levels=c(letters[1:N_expo_agg()],"aggregated")),
  #                     err = err)
  #   dat <- melt(dat,id.vars = "predictor",value.name="error") %>% group_by(predictor) %>% summarise(MSE=mean(error))
  #   mse <- dat %>% filter(predictor=="aggregated") %>% select(MSE) %>% as.numeric()
  #   print(mse)
  #   ggplot(dat %>% filter(predictor != "aggregated"))+
  #     aes(x=predictor,y=MSE,group=1)+
  #     geom_point()+
  #     geom_line()+
  #     geom_hline(yintercept=mse)+
  #     xlab("Predictor")+ylab("Predictor estimation MSE")+
  #     scale_y_log10()+
  #     theme(axis.text=element_text(size=12),
  #           axis.title=element_text(size=14,face="bold"))
  # })
  # 
  # output$pred_diag_expo_agg <- renderPlot({
  #   df <- data.frame(predictor = factor(c(letters[1:N_expo_agg()],"aggregated"),
  #                                       levels=c(letters[1:N_expo_agg()],"aggregated")) ,
  #                    expo_agg_predictions())
  #   df <- melt(df,id.vars=c("predictor"),value.name="pred")
  #   df <- data.frame(df,best_pred = rep(best_predictions_expo_agg(),each=N_expo_agg()+1))
  #   df <- df %>% filter(predictor %in% input$pred_select_expo_agg)
  #   ggplot(df)+
  #     geom_point(aes(x=pred,y=best_pred, color=predictor))+
  #     ylab("Best predictions") + xlab("Predicted values")+
  #     geom_abline(slope=1,intercept=0,aes(linetype="dashed"))+
  #     theme(axis.text=element_text(size=12),
  #           axis.title=element_text(size=14,face="bold"),
  #           legend.title=element_text(size=12),
  #           legend.text=element_text(size=12),
  #           legend.position="bottom")
  # })
  
  # Expo agg customized
  
  cust_expo_agg_strategy <- reactive({
    if (input$cust_expo_agg_strat=="Quadratic loss"){
      2
    }else{
      1
    }
  })
  
  output$cust_expo_agg_pred_select <- renderUI({
    selectInput("pred_select_cust_expo_agg","Predictions to plot",
                choices=c("aggregated",letters[1:input$cust_expo_agg_N]), 
                selected="aggregated", 
                multiple=TRUE, selectize = TRUE)
  })

  
  output$cust_expo_agg_params_table <- renderTable({
    df <- data.frame(c("N","eta"),c(as.integer(input$cust_expo_agg_N),input$cust_expo_agg_eta))
    colnames(df) <- c("parameter","value")
    df
  })
  
  cust_expo_agg_mu <- reactive({
    req(input$cust_expo_agg_log_mu_max,input$cust_expo_agg_log_mu_min,input$cust_expo_agg_N)
    log_mu <- (input$cust_expo_agg_log_mu_max-input$cust_expo_agg_log_mu_min)/(input$cust_expo_agg_N-1)*c(0:(input$cust_expo_agg_N-1)) + input$cust_expo_agg_log_mu_min
    10**log_mu
  })
  
  output$cust_expo_agg_mu_table <- renderTable({
    req(cust_expo_agg_mu())
    data.frame(predictor = letters[1:input$cust_expo_agg_N], mu = cust_expo_agg_mu(), `log mu` = log10(cust_expo_agg_mu()))
  })
  
  cust_expo_agg_res <- reactive({
    req(tvar_sample(),cust_expo_agg_mu())
    predict_exp_agg_NLMS(tvar_sample(),cust_expo_agg_strategy(),input$cust_expo_agg_eta,cust_expo_agg_mu(),isolate(input$order))
  })
  
  cust_expo_agg_predictions <- reactive({
    cust_expo_agg_res()[[1]]})
  
  
  best_predictions_cust_expo_agg <- eventReactive(input$simulate,
                                             {apply(tvar_sample(),MARGIN = 2,
                                                    function(x) best_predict_TVAR(x, 2**input$log_Tmax,rd_seeds(),input$delta))},ignoreNULL = FALSE)
  
  # Plots
  
  output$pred_mse_cust_expo_agg <- renderPlot({
    req(cust_expo_agg_predictions(),best_predictions_cust_expo_agg())
    err <- t((t(cust_expo_agg_predictions())-best_predictions_cust_expo_agg())**2)
    dat <- data.frame(predictor = factor(c(letters[1:input$cust_expo_agg_N],"aggregated"),
                                         levels=c(letters[1:input$cust_expo_agg_N],"aggregated")),
                      err = err)
    dat <- melt(dat,id.vars = "predictor",value.name="error") %>% group_by(predictor) %>% summarise(MSE=mean(error))
    mse <- dat %>% filter(predictor=="aggregated") %>% select(MSE) %>% as.numeric()
    p<- ggplot(dat %>% filter(predictor != "aggregated"))+
      aes(x=predictor,y=MSE,group=1)+
      geom_point()+
      geom_line()+
      geom_hline(yintercept=mse)+
      xlab("Predictor")+ylab("Predictor estimation MSE")+
      scale_y_log10()+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    plots$expo_agg$pred_mse <- p
    p
  })
  
  output$pred_diag_cust_expo_agg <- renderPlot({
    df <- data.frame(predictor = factor(c(letters[1:input$cust_expo_agg_N],"aggregated"),
                                        levels=c(letters[1:input$cust_expo_agg_N],"aggregated")) ,
                     cust_expo_agg_predictions())
    df <- melt(df,id.vars=c("predictor"),value.name="pred")
    df <- data.frame(df,best_pred = rep(best_predictions_cust_expo_agg(),each=input$cust_expo_agg_N+1))
    df <- df %>% filter(predictor %in% input$pred_select_cust_expo_agg)
    p <- ggplot(df)+
      geom_point(aes(x=pred,y=best_pred, color=predictor))+
      ylab("Best predictions") + xlab("Predicted values")+
      geom_abline(slope=1,intercept=0,aes(linetype="dashed"))+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            legend.title=element_text(size=12),
            legend.text=element_text(size=12),
            legend.position="bottom")
    plots$expo_agg$pred_diag <- p
    p
  })
  
  output$alpha_cust_expo_agg <- renderPlot({
    df <- data.frame(predictor = factor(c(letters[1:input$cust_expo_agg_N])) ,
                     cust_expo_agg_res()[[2]])
    df <- melt(df,id.vars=c("predictor"), value.name = "alpha")
    p <- ggplot(df)+
      aes(x=predictor, y=alpha)+
      geom_boxplot()+
      geom_point()+
      coord_cartesian(ylim=c(0,1))+
      xlab("Predictor")+ylab("Predictor weight")+
      geom_hline(yintercept = 1/input$cust_expo_agg_N)+
      theme(axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"))
    plots$expo_agg$pred_weights <- p
    p
      
  })
  
  # Comparison of predictors
  
  compare<- reactiveValues()
  
  observeEvent(input$simulate,
               {               compare$data <- data.frame(algorithm = c("Local Yule-Walker",
                                                                        "Local Yule-Walker",
                                                                        "NLMS",
                                                                        "NLMS",
                                                                        "Exponential aggregation"),
                                                          method = c("Normal",
                                                                     "Romberg",
                                                                     "Normal",
                                                                     "Romberg",
                                                                     "Normal"),
                                                          estim_log_mse = rep(NA,5),
                                                          pred_log_mse = rep(NA,5))
               
               compare$pred <- data.frame(replicate(6,rep(NA,input$n_obs)))
               },
               ignoreNULL = FALSE
               )
  
  
  yw_estimations <- reactive({
    req(input$log_M)
    do.call(cbind,lapply(coeffs_by_M(), function(x) x[[input$log_M+1]]))})
  
  yw_predictions <- reactive({
    predictions_by_M()[input$log_M+1,]})
  
  observeEvent(yw_estimations(),
               {mse <- mean(colSums((yw_estimations() - reg_coeffs()(1))**2))
                compare$data[1,3] <- log10(mse)
               })
  
  observeEvent(yw_romberg_estimations(),
               {mse <- mean(colSums((yw_romberg_estimations() - reg_coeffs()(1))**2))
               compare$data[2,3] <- log10(mse)
               })
  
  observeEvent(nlms_estimations(),
               {mse <- mean(colSums((nlms_estimations() - reg_coeffs()(1))**2))
               compare$data[3,3] <- log10(mse)
               })
  
  observeEvent(nlms_romberg_estimations(),
               {mse <- mean(colSums((nlms_romberg_estimations() - reg_coeffs()(1))**2))
               compare$data[4,3] <- log10(mse)
               })
  
  observeEvent(yw_predictions(),
               {mse <- mean((yw_predictions()-best_predictions_yw())**2)
                compare$data[1,4] <- log10(mse)
                compare$pred[,1] <- yw_predictions()
               })
  
  observeEvent(yw_romberg_predictions(),
               {mse <- mean((yw_romberg_predictions()-best_predictions_yw())**2)
               compare$data[2,4] <- log10(mse)
               compare$pred[,2] <- yw_romberg_predictions()
               })
  
  observeEvent(nlms_predictions(),
               {mse <- mean((nlms_predictions()-best_predictions_nlms())**2)
               compare$data[3,4] <- log10(mse)
               compare$pred[,3] <- nlms_predictions()
               })
  
  observeEvent(nlms_romberg_predictions(),
               {mse <- mean((nlms_romberg_predictions()-best_predictions_nlms())**2)
               compare$data[4,4] <- log10(mse)
               compare$pred[,4] <- nlms_romberg_predictions()
               })
  
  observeEvent(cust_expo_agg_predictions(),
               {mse <- mean((cust_expo_agg_predictions()[input$cust_expo_agg_N+1,]-best_predictions_cust_expo_agg())**2)
               compare$data[5,4] <- log10(mse)
               compare$pred[,5] <- cust_expo_agg_predictions()[input$cust_expo_agg_N+1,]
               })
  
  output$compare_table <- renderTable(
    {df <- compare$data
    idx <- NULL
    if ("Normal" %in% input$yw_compare){
      idx <- c(idx,1)
    }
    if ("Romberg" %in% input$yw_compare){
      idx <- c(idx,2)
    }
    if ("Normal" %in% input$nlms_compare){
      idx <- c(idx,3)
    }
    if ("Romberg" %in% input$nlms_compare){
      idx <- c(idx,4)
    }
    if ("Normal" %in% input$expo_agg_compare){
      idx <- c(idx,5)
    }
    colnames(df) <- c("Predictor","Method","Estimation log MSE","Prediction log MSE")
    if(length(idx)>0){
      p <- df[idx,]
    }else{
      p <-NULL
      
    }
    compare$table <- p
    p
})
  
  output$compare_pred_diag <- renderPlot({
    pred <- compare$pred
    names <- c("yule-walker (normal)",
               "yule-walker (romberg)",
               "nlms (normal)",
               "nlms (romberg )",
               "exponential aggregation (normal)")
    idx <- NULL
    if ("Normal" %in% input$yw_compare){
      idx <- c(idx,1)
    }
    if ("Romberg" %in% input$yw_compare){
      idx <- c(idx,2)
    }
    if ("Normal" %in% input$nlms_compare){
      idx <- c(idx,3)
    }
    if ("Romberg" %in% input$nlms_compare){
      idx <- c(idx,4)
    }
    if ("Normal" %in% input$expo_agg_compare){
      idx <- c(idx,5)
    }
    idx <- setdiff(idx,which(colSums(is.na(pred))>0))
    colnames(pred) <- c("Predictor","Method","Estimation log MSE","Prediction log MSE")
    if(length(idx)>0){
      pred <- as.data.frame(cbind(best_predictions_yw(),
                    pred[,idx]))
      colnames(pred) <- c("best_pred",names[idx])
      pred <- melt(pred,id.vars=c("best_pred"),value.name="pred",variable.name="method")
      p <- ggplot(pred)+
        aes(x=pred,y=best_pred,color=method)+
        geom_point()+
        xlab("Predicted values") + ylab("Best predictors")+labs(color = "method")+
        geom_abline(slope=1,intercept=0,aes(linetype="dashed"))+
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.title=element_text(size=12),
              legend.text=element_text(size=12),
              legend.position="bottom")
    }else{
      p<-NULL
    }
    compare$pred_diag <- p
    p
  })

  # Downloads handling
  output$dl_plots_settings <- downloadHandler(filename = "TVAR_settings.zip",
                                              content = function(file) {
                                                owd <- setwd(tempdir())
                                                on.exit(setwd(owd))
                                                files = c("TVAR_coefficients.png","TVAR_sample.png","TVAR_subsample.png","ACF_PACF_subsample.png")
                                                ggsave(files[1],plots$settings$coeffs, width = 20, height = 10 * ceiling(input$order / 2), unit = "cm")
                                                ggsave(files[2],plots$settings$sample, width = 20, height = 10, unit = "cm")
                                                ggsave(files[3],plots$settings$subsample, width = 20, height = 10, unit = "cm")
                                                ggsave(files[4],plots$settings$acf_pacf, width = 20, height = 10, unit = "cm")
                                                zip(file,files)
                                              })
  
  output$dl_plots_yw <- downloadHandler(filename = "yw_plots.zip",
                                        content = function(file) {
                                          owd <- setwd(tempdir())
                                          on.exit(setwd(owd))
                                          if(input$romberg_yw){
                                            files = c("yw_romb_estim_mse.png",
                                                      "yw_romb_pred_mse.png",
                                                      "yw_romb_pred_diag.png")
                                            ggsave(files[1],plots$yw$estim_mse, width = 20, height = 10, unit = "cm")
                                            ggsave(files[2],plots$yw$pred_mse, width = 20, height = 10, unit = "cm")
                                            ggsave(files[3],plots$yw$romb_pred_diag, width = 20, height = 20, unit = "cm")
                                          }else{
                                            files = c("yw_estim_mse.png",
                                                      "yw_pred_mse.png",
                                                      "yw_pred_diag.png")
                                            ggsave(files[1],plots$yw$estim_mse, width = 20, height = 10, unit = "cm")
                                            ggsave(files[2],plots$yw$pred_mse, width = 20, height = 10, unit = "cm")
                                            ggsave(files[3],plots$yw$pred_diag, width = 20, height = 20, unit = "cm")
                                          }
                                          zip(file,files)})
  
  output$dl_plots_nlms <- downloadHandler(filename = "nlms_plots.zip",
                                          content = function(file) {
                                            owd <- setwd(tempdir())
                                            on.exit(setwd(owd))
                                            if(input$romberg_nlms){
                                              files = c("nlms_romb_estim_mse.png",
                                                        "nlms_romb_pred_mse.png",
                                                        "nlms_romb_pred_diag.png")
                                              ggsave(files[1],plots$nlms$romb_estim_mse, width = 20, height = 10, unit = "cm")
                                              ggsave(files[2],plots$nlms$romb_pred_mse, width = 20, height = 10, unit = "cm")
                                              ggsave(files[3],plots$nlms$romb_pred_diag, width = 20, height = 20, unit = "cm")
                                            }else{
                                              files = c("nlms_estim_mse.png",
                                                        "nlms_pred_mse.png",
                                                        "nlms_pred_diag.png")
                                              ggsave(files[1],plots$nlms$estim_mse, width = 20, height = 10, unit = "cm")
                                              ggsave(files[2],plots$nlms$pred_mse, width = 20, height = 10, unit = "cm")
                                              ggsave(files[3],plots$nlms$pred_diag, width = 20, height = 20, unit = "cm")
                                            }
                                            zip(file,files)})
  
  output$dl_plots_expo_agg <- downloadHandler(filename = "expo_agg_plots.zip",
                                              content = function(file) {
                                                owd <- setwd(tempdir())
                                                on.exit(setwd(owd))
                                                files = c("expo_agg_pred_mse.png",
                                                          "expo_agg_pred_diag.png",
                                                          "expo_agg_pred_weights.png")
                                                ggsave(files[1],plots$expo_agg$pred_mse, width = 20, height = 10, unit = "cm")
                                                ggsave(files[2],plots$expo_agg$pred_diag, width = 20, height = 20, unit = "cm")
                                                ggsave(files[3],plots$expo_agg$pred_weights, width = 20, height = 10, unit = "cm")
                                                zip(file,files)})
  
  output$dl_compare <- downloadHandler(filename = "compare.zip",
                                              content = function(file) {
                                                owd <- setwd(tempdir())
                                                on.exit(setwd(owd))
                                                files = c("mse_table.txt",
                                                          "pred_diag.png")
                                                ggsave(files[2],compare$pred_diag, width = 20, height = 20, unit = "cm")
                                                if (!is.null(compare$table)){
                                                 print.xtable(xtable(compare$table),file=files[1])
                                                 zip(file,files)
                                                }else{
                                                  zip(file,files[2])
                                                }
                                                
                                                })
  
})

