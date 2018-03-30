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

navbarPage("TVAR series forecasting",
           tabPanel("TVAR model",
                    sidebarLayout(
                      sidebarPanel(
                        h3("Model parameters"),
                        numericInput("order",
                                     "Order",
                                     value=1,
                                     min=1),
                        tags$hr(),
                        h4("Coefficients"),
                        numericInput("richness",
                                     "Richness",
                                     value = 1,
                                     min=1),
                        sliderInput("delta",
                                    "delta",
                                    min = 0.01, max = 1,step = 0.01, value = .99,ticks =FALSE),
                        actionButton("roll_dice","Roll the dice"),
                        tags$hr(),
                        h4("White noise"),
                        selectInput("white_noise_gen",
                                    "Distribution",
                                    c("Normal","Laplace","Rademacher")),
                        sliderInput("sigma","sigma",min=0,max=10,step=0.1,value=1,ticks =FALSE),
                        tags$hr(),
                        downloadButton("dl_plots_settings","Download plots")
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("TVAR coefficients",
                                   plotOutput("coeffs_plot")
                          ),
                          tabPanel("Simulation settings",
                                   sidebarLayout(
                                     sidebarPanel(width = 8,
                                                  sliderInput("log_Tmax",
                                                              "log length",
                                                              min = 0, max = 20, 
                                                              value = 10, step = 1, ticks=FALSE),
                                                  sliderInput("n_obs","Number of simulations",
                                                              min=100,max=500,value=500, ticks =FALSE),
                                                  actionButton("simulate","New simulation")
                                     ),
                                     mainPanel(
                                       plotOutput("sample_plot")
                                     )
                                   )
                          ),
                          tabPanel("Local stationarity",
                                   sidebarLayout(
                                     sidebarPanel(width=8,
                                                  uiOutput("sample_id_slider"),
                                                  uiOutput("window_slider")
                                     ),
                                     mainPanel(
                                       plotOutput("subsample_plot"),
                                       plotOutput("acf_pacf_plot")
                                     )
                                   )
                          )
                        )
                      )
                    )
           ),
           tabPanel("Local Yule-Walker",
                    sidebarLayout(
                      sidebarPanel(
                        h3("Configuration"),
                        selectInput("taper_pred",
                                    "Tapering function",
                                    c("Rectangular","Bartlett","Hann","Hamming","Blackman")),
                        plotOutput("taper_plot_pred",height=100,width=250),
                        tags$hr(),
                        uiOutput("log_bandwidth_slider"),
                        checkboxInput("romberg_yw","Plot Romberg's bias reduction"),
                        conditionalPanel("input.romberg_yw",
                                         uiOutput("log_bandwidth_slider_bandwidth_romberg")
                        ),
                        actionButton("predict_yw","Predict"),
                        tags$hr(),
                        downloadButton("dl_plots_yw","Download plots")
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Estimation",
                                   plotOutput("estim_mse_yw_plot")
                          ),
                          tabPanel("Prediction",
                                   plotOutput("pred_mse_yw_plot"),
                                   conditionalPanel("input.romberg_yw",
                                                    plotOutput("romb_pred_diag_yw_plot")
                                   ),
                                   conditionalPanel("!input.romberg_yw",
                                                    plotOutput("pred_diag_yw_plot")
                                   )
                          )
                        )
                      )
                    )
           ),
           tabPanel("Normalized Least Mean Squares",
                    sidebarLayout(
                      sidebarPanel(
                        h3("Configuration"),
                        splitLayout(
                          numericInput("log_mu_min","min log mu",max = 0, value = -6),
                          numericInput("log_mu_max","max log mu",max = 0, value = 0)
                        ),
                        uiOutput("nlms_log_mu_input"),
                        uiOutput("nlms_plot_window_input"),
                        checkboxInput("romberg_nlms","Plot Romberg's bias reduction"),
                        conditionalPanel("input.romberg_nlms",
                                         uiOutput("nlms_romb_log_mu_input"),
                                         sliderInput("gamma","gamma",min=0.01,max=.99,ticks=FALSE, value=0.5, step = .05,animate=animationOptions(interval=250))
                        ),
                        actionButton("nlms_flush","Flush graphics"),
                        tags$hr(),
                        downloadButton("dl_plots_nlms","Download plots")
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Estimation",
                                   conditionalPanel("!input.romberg_nlms",
                                                    plotOutput("estim_mse_nlms_plot")
                                   ),
                                   conditionalPanel("input.romberg_nlms",
                                                    renderTable("foo"),
                                                    plotOutput("romb_estim_mse_nlms_plot")
                                   )
                          ),
                          tabPanel("Prediction",
                                   conditionalPanel("!input.romberg_nlms",
                                                    splitLayout(
                                                      plotOutput("pred_mse_nlms_plot"),
                                                      plotOutput("pred_diag_nlms_plot")
                                                    )
                                                    
                                   ),
                                   conditionalPanel("input.romberg_nlms",
                                                    splitLayout(
                                                      plotOutput("romb_pred_mse_nlms_plot"),
                                                      plotOutput("romb_pred_diag_nlms_plot")
                                                    )
                                   )
                          )
                        )
                      )
                    )
           ),
           tabPanel("NLMS exponential aggregation",
                    # tabsetPanel(
                    #   tabPanel("Default parameters",
                    #            sidebarLayout(
                    #              sidebarPanel(
                    #                h3("Configuration"),
                    #                sliderInput("beta_expo_agg","Regularity parameter (beta)",min=0,max=1,value=.5, ticks=FALSE),
                    #                tags$hr(),
                    #                tableOutput("expo_agg_params_table"),
                    #                tableOutput("expo_agg_mu_table"),
                    #                tags$hr(),
                    #                uiOutput("expo_agg_pred_select")
                    #              ),
                    #              mainPanel(
                    #                plotOutput("pred_mse_expo_agg"),
                    #                plotOutput("pred_diag_expo_agg")
                    #              )
                    #            )
                    #   ),
                    #   tabPanel("Customized parameters",
                    sidebarLayout(
                      sidebarPanel(
                        h3("Configuration"),
                        selectInput("cust_expo_agg_strat","Aggregation strategy",
                                    choices = c("Quadratic loss","Gradient of quadratic loss")),
                        numericInput("cust_expo_agg_N","Number of predictors",min=0,step=1,value=7),
                        sliderInput("cust_expo_agg_eta","Temperature (eta)",
                                    min = 0,
                                    max = 1,
                                    step=.01,
                                    value = 0,
                                    ticks=FALSE),
                        splitLayout(
                          numericInput("cust_expo_agg_log_mu_min","min log mu",max = 0, value = -6),
                          numericInput("cust_expo_agg_log_mu_max","max log mu",max = 0, value = 0)
                        ),
                        tags$hr(),
                        uiOutput("cust_expo_agg_pred_select"),
                        tags$hr(),
                        downloadButton("dl_plots_expo_agg","Download plots")
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Parameters",
                                   splitLayout(
                                     tableOutput("cust_expo_agg_params_table"),
                                     tableOutput("cust_expo_agg_mu_table")
                                   )
                          ),
                          tabPanel("Aggregation",
                                   splitLayout(
                                     plotOutput("pred_mse_cust_expo_agg"),
                                     plotOutput("alpha_cust_expo_agg"))
                          ),
                          tabPanel("Prediction",
                                   plotOutput("pred_diag_cust_expo_agg")
                          )
                        )
                      )
                    )
                    #)
                    #)
                    
           ),
           tabPanel("Comparison of predictors",
                    sidebarLayout(
                      sidebarPanel(
                        h3("Predictors to compare"),
                        tags$hr(),
                        h4("Local Yule-Walker"),
                        selectInput("yw_compare","",
                                    choices=c("Normal",
                                              "Romberg"),
                                    multiple=TRUE),
                        h4("Normalized Least Mean Squares"),
                        selectInput("nlms_compare","",
                                    choices=c("Normal",
                                              "Romberg"),
                                    multiple=TRUE),
                        h4("Exponential aggregation"),
                        selectInput("expo_agg_compare","",
                                    choices=c("Normal"),
                                    multiple=TRUE),
                        tags$hr(),
                        downloadButton("dl_compare","Download comparison")
                      ),
                      mainPanel(
                        tabsetPanel(
                          tabPanel("Mean squared errors",
                                   tableOutput("compare_table")),
                          tabPanel("Diagnosis plot",
                                   plotOutput("compare_pred_diag"))
                        )
                      )
                    )
           ),
           tabPanel("About",
                    includeHTML("about.html")
                    
           )
)
