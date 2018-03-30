#include <Rcpp.h>
using namespace Rcpp;


const double pi = 3.14159265358979323846;

// Levinson-Durbin algorithm to solve Yule_Walker equation
// [[Rcpp::export]]
NumericVector solve_yule_walker(NumericVector gamma){
  int K = gamma.size()-1;
  NumericVector out(K);
  NumericMatrix phi(K,K);
  double acc;
  
  if (gamma[0] !=0){
    double kappa = gamma[1]/gamma[0];
    phi(0,0) = kappa;
    double sig2 = gamma[0]*(1-std::pow(kappa,2));
    if (sig2==0) return(out);
    for(int p = 1; p <= K-1; ++p){
      if (sig2==0) return(out);
      acc = 0;
      for(int k = 1; k <= p; k++){
        acc+=phi(k-1,p-1)*gamma[p+1-k];
      }
      kappa = (gamma[p+1] - acc)/sig2;
      sig2 = sig2*(1-std::pow(kappa,2));
      phi(p,p) = kappa;
      for(int m = 1; m <= p; ++m){
        phi(m-1,p) = phi(m-1,p-1) - kappa*phi(p-m,p-1);
      }
    }
    for(int i = 0; i < K; ++i){
      out[i] = phi(i,K-1);
    }
  }
  return(out);
}

double rectangular(double u){
  if (u>1 || u<0){
    return(0.);
  }
  else{
    return(1.);
  };
}

double bartlett(double u){
  if (u>1 || u<0){
    return(0.);
  }
  else if (u<.5){
    return(2*u);
  }
  else{
    return(2*(1-u));
  };
}

double hann(double u){
  if (u>1||u<0){
    return(0.);
  }
  else{
    return(0.5*(1-std::cos(2*pi*u)));
  }
}  

double hamming(double u){
  if (u<0||u>1){
    return(0.);
  }
  else{
    return(0.54-0.46*std::cos(2*pi*u));
  }
} 

double blackman(double u){
  if (u<0||u>1){
    return(0.);
  }
  else{
    return(0.42-0.5*std::cos(2*pi*u) + 0.08*cos(4*pi*u));
  }
}

// [[Rcpp::export]]
double compute_taper(double u, String taper, bool causal){
  if (causal && u >.5){
    return(0);
  }
  double out;
  if (taper=="Rectangular"){
    out = rectangular(u);
  }
  else if(taper=="Bartlett"){
    out = bartlett(u);
  }
  else if(taper=="Hann"){
    out = hann(u);
  }
  else if(taper=="Hamming"){
    out = hamming(u);
  }
  else if (taper=="Blackman"){
    out = blackman(u);
  }
  else{
    out = 0;
  }
  return(out);
}

// [[Rcpp::export]]
NumericVector estimate_gamma(NumericVector sample, int t,int order,int bandwidth,String taper, bool causal){
  NumericVector gamma(order+1);
  int Tmax = sample.size() - 1;
  int M = bandwidth;
  int kmin = std::max(1,M/2-t);
  for(int l = 0; l <= order; ++l){
    int kmax = std::min(Tmax + M/2 - l - t, M-l);
    double hM = 0;
    for(int k = kmin; k <=kmax; ++k){
      hM += std::pow(compute_taper((double) k /M,taper,causal),2);
      gamma[l] += compute_taper( (double)k /M,taper,causal) * compute_taper( (double)(k+l) /M,taper,causal) * sample[t+k-M/2] * sample[t+k+l-M/2];
    }
    if (hM>0) gamma[l] /= hM;
  }
  return(gamma);
}

// Estimation of the TVAR coefficients with local Yule Walker method
// [[Rcpp::export]]
NumericVector estimate_coeffs_yule_walker(NumericVector sample,int t,int order,int bandwidth,String taper, bool causal){
  NumericVector gamma = estimate_gamma(sample, t, order, bandwidth,taper, causal);
  return(solve_yule_walker(gamma));
}

// Predict the TVAR series at time t with local Yule-Walker method
// [[Rcpp::export]]
double predict_TVAR_yule_walker(NumericVector sample, int t, int order, int bandwidth, String taper, bool causal){
  NumericVector coeffs = estimate_coeffs_yule_walker(sample,t,order,bandwidth,taper, causal);
  double res = 0;
  for(int l = 0; (l < order) && (t - l >=0) ; ++l){
    res += coeffs[l]*sample[t-l];
  }
  return(res);
}

// Compute the Romberg acceleration weights
// [[Rcpp::export]]
NumericVector romberg_weights(int k){
  int M = 2;
  int R = k+1;
  NumericVector w(k+1);
  for(int i = 1; i<=R;++i ){
    double nom = std::pow(M,.5*(R-i)*(R-i+1));
    double p1 = 1;
    for(int j = 1; j<=i-1; ++j){
      p1 *= (1-std::pow(M,j));
    }
    double p2 = 1;
    for(int j = 1; j<=R-i; ++j){
      p2 *= (1-std::pow(M,j));
    }
    w[i-1] = std::pow(-1,R-i) * nom /(p1*p2);
  }
  return(w);
}

// Predict the TVAR series at time t with local Yule-Walker method with Romberg acceleration
// [[Rcpp::export]]
double predict_TVAR_yule_walker_romberg(NumericVector sample, int t, int order, int max_bandwidth, int min_bandwidth, String taper, bool causal){
  int k = (int) floor(log2(max_bandwidth/min_bandwidth));
  NumericVector weights = romberg_weights(k);  
  NumericVector coeffs(order);
  int bandwidth = min_bandwidth;
  for(int j =0; j<=k;j++){
    coeffs += weights[j] * estimate_coeffs_yule_walker(sample,t,order,bandwidth,taper, causal);
    bandwidth *= 2;
  }
  double res = 0;
  for(int l = 0; (l < order) && (t - l >=0) ; ++l){
    res += coeffs[l]*sample[t-l];
  }
  return(res);
}

// [[Rcpp::export]]
List estimate_coeffs_yule_walker_by_M(NumericMatrix samples, int t,int order,int max_log_bandwidth,String taper, bool causal){
  int n_obs = samples.ncol();
  List res(n_obs);
  for (int i = 0; i < n_obs; ++i){
    NumericVector sample = samples( _, i);
    List tmp(max_log_bandwidth+1);
    int M  = 1;
    for (int j = 0; j <= max_log_bandwidth; ++j){
      tmp[j] = estimate_coeffs_yule_walker(sample,t,order,M,taper, causal);
      M*=2;
    }
    res[i] = tmp;
  }
  return(res);
}
