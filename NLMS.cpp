#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix predict_NLMS(NumericMatrix sample, double mu, int order){
  int n_obs = sample.ncol();
  int n = sample.nrow();
  NumericMatrix theta(order,n_obs);
  for(int k = 0; k<n_obs; ++k){
    for(int t = 0; t<n-1; ++t){
      double num = 0;
      double denom = 0;
      for(int j = 0; j<order && t-j>=0; ++j){
        num += sample(t-j,k)*theta(j,k);
        denom += std::pow(sample(t-j,k),2);
      }
      num = mu * (sample(t+1,k)-num);
      denom = 1 + mu * denom;
      double q = num/denom;
      for(int j = 0; j<order; ++j){
        if (t-j>=0){
          theta(j,k) += q * sample(t-j,k);
        }
      }
    }
  }
  return(theta);
}

