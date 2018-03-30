#include <Rcpp.h>
using namespace Rcpp;

// Levison Durbin algorithm
// to generate admissible AR coefficients
// from partial autocorrelation coefficients
NumericVector levinson_durbin_acr(NumericVector acr, double delta){ 
  int p = acr.size();
  NumericVector out(p);
  NumericMatrix thetas(p, p);
  for(int j = 0; j < p; ++j){ 
      thetas(j,j) = acr[j];
  }
  
  for(int k = 1; k < p; ++k) {
      for(int j = 0; j <= k-1; ++j) {
        thetas(j,k) = thetas(j,k-1) - thetas(k,k)*thetas(k-j,k-1);
    }
  }
  for(int j = 0; j < p; ++j) {
      out[j] = std::pow(delta,(j+1)) * thetas(j,p-1);
    }
  return(out);
}

// [[Rcpp::export]]
NumericVector compute_reg_coeffs(double t, NumericMatrix seeds, double delta) {
  int F = seeds.nrow() + 1;
  double div = F*(F-1)*(2*F-1)/6;
  int order = seeds.ncol();
  NumericVector acr(order);
  for(int k = 0; k<order; ++k){
    for(int j=0; j < F-1; ++j){
      acr[k] += seeds(j,k)*std::pow((j+1),2)*std::cos((j+1)*t);
    }
    acr[k] /= div;
  }
  NumericVector coeffs = levinson_durbin_acr(acr,delta);
  return(coeffs);
}


NumericVector innovation_gen(int n, String white_noise_distrib){
  
  Rcpp::Environment package_env("package:extraDistr"); 
     
  
  if (white_noise_distrib=="Normal")
  {
    return(rnorm(n));
  } 
  else if (white_noise_distrib=="Laplace"){
    Rcpp::Function simul = package_env["rlaplace"]; 
    return(simul(n));
  }
  else if (white_noise_distrib=="Rademacher"){
    Rcpp::Function simul = package_env["rsign"]; 
    return(simul(n));
  }
  else
  {
    return(NumericVector(n));
  }
}

// [[Rcpp::export]]
NumericMatrix TVAR_gen(int Tmax, int n_obs, double sigma, String white_noise_gen, NumericMatrix seeds, double delta){
  
  int order = seeds.ncol();
  
  NumericMatrix simul(Tmax+1,n_obs);
  NumericVector coeffs(order);
  NumericVector rd(n_obs);
  
  rd = innovation_gen(n_obs,white_noise_gen);
  for(int k = 0; k < n_obs; ++k){
    simul(0,k) = rd[k];
  }
  
  for(int i = 0; i < Tmax; ++i){
    coeffs = compute_reg_coeffs((double) i/Tmax,seeds,delta);
    rd = innovation_gen(n_obs,white_noise_gen);
    for(int k = 0; k < n_obs; ++k){
      for(int l = 0; (l < order) && (i - l >=0) ; ++l){
        simul(i+1,k) += coeffs[l]*simul(i-l,k);
      }
      simul(i+1,k) += sigma * rd[k];
    }
  }
  return(simul);
}

// [[Rcpp::export]]
double best_predict_TVAR(NumericVector sample, int t, NumericMatrix seeds, double delta){
  int T_max  = sample.size() - 1;
  NumericVector coeffs = compute_reg_coeffs((double)t / T_max ,seeds,delta);
  int order = coeffs.size();
  double res = 0;
  for(int l = 0; (l < order) && (t - l >=0) ; ++l){
    res += coeffs[l]*sample[t-l];
  }
  return(res);
}

// Predict the TVAR series at time t+1 using the given regression coefficients
// [[Rcpp::export]]
double predict_TVAR(NumericVector sample, int t, NumericVector coeffs){
  int order = coeffs.size();
  double res = 0;
  for(int l = 0; (l < order) && (t - l >=0) ; ++l){
    res += coeffs[l]*sample[t-l];
  }
  return(res);
}

// [[Rcpp::export]]
NumericVector predict_mult_TVAR(NumericMatrix sample, int t, NumericMatrix coeffs){
  int n_obs = sample.ncol();
  int order = coeffs.nrow();
  NumericVector res(n_obs);
  for(int i=0; i<n_obs;++i){
    for(int l = 0; (l < order) && (t - l >=0) ; ++l){
      res[i] += coeffs(l,i)*sample(t-l,i);
    }
  }
  return(res);
}


