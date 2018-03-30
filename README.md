# TVAR forecasting

This is a Shiny application for time-varying autoregressive (TVAR) series forecasting. 

This application illustrates the main techniques which can be used when one wants to predict a TVAR series. 

It enables simulations of TVAR series and comparisons between the different methods by Monte-Carlo.

The different methods which are implemented are discussed in the report.pdf file. This also contains references for further readings.

## Dependencies

The software is based on R and C++. You must have R installed on your computer and a C++ compiler compatible with Rcpp.

The following R packages are required to run the application:

- dplyr
- extraDistr
- forecast
- ggforce
- ggplot2
- gridExtra
- Rcpp
- reshape2
- rmutil
- shiny
- xtable

One can easily install them by running the following command in the R console:

```r
install.packages("dplyr","extraDistr","forecast","ggforce","ggplot2","gridExtra",
                 "Rcpp","reshape2","rmutil","shiny","xtable")
```

## Running the application from Github

When the dependencies have been installed, one can simply run the application locally from the R console using:

```r
shiny::runGitHub("tvar_forecasting", "maximegodin") 
```
