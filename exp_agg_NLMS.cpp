//-------------------------------------------------------------------------
//  This file is part of the TVAR forecasting application https://github.com/maximegodin/tvar_forecasting
//
//  This application is governed by the MIT license.
//  You can  use, modify and/ or redistribute this code under the terms
//  of the MIT license:  https://opensource.org/licenses/MIT
//
//  Maxime Godin, Amina Bouchafaa, Ecole Polytechnique
//  March 29th, 2018
//-------------------------------------------------------------------------

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List predict_exp_agg_NLMS(NumericMatrix sample, int strategy, double eta, NumericVector mu, int order){
  int n_obs = sample.ncol();
  int n = sample.nrow();
  int N = mu.size();
  NumericMatrix alpha(N,n_obs);
  NumericMatrix res(N+1,n_obs);
  
  // for each sample
  for(int k = 0; k<n_obs; ++k){
    // initialize the weights
    for(int i=0; i<N; ++i){
      alpha(i,k) = 1./N;
    }
    
    // with no prior information we predict by 0 the next observation
    double agg_pred = 0;
    NumericMatrix theta(order,N);
    NumericVector pred(N);
    
    // for each new observation
    for(int t = 0; t<n-1; ++t){
      double s = 0;
      // for each predictor
      for(int i=0; i<N; ++i){
        // update the weights using the new observation and according to the strategy
        if (strategy==1){
          alpha(i,k)*=std::exp(-2*eta*(agg_pred-sample(t+1,k))*pred[i]);
        }
        else{
          alpha(i,k)*= std::exp(-eta*std::pow(pred[i]-sample(t+1,k),2));
        }
        s+= alpha(i,k);
      }
      
      // we now compute the prediction for the next observation
      agg_pred = 0;
      for(int i=0; i<N; ++i){
        alpha(i,k) /= s; // rescale the weights
        pred[i] =0; 
        
        // update our theta estimator
        double num = 0;
        double denom = 0;
        for(int j = 0; j<order && t-j>=0; ++j){
          num += sample(t-j,k)*theta(j,i);
          denom += std::pow(sample(t-j,k),2);
        }
        num = mu[i] * (sample(t+1,k)-num);
        denom = 1 + mu[i] * denom;
        double q = num/denom;
        for(int j = 0; j<order; ++j){
          if (t-j>=0){
            theta(j,i) += q * sample(t-j,k);
          }
        }
        
        // compute the prediction for this predictor
        for(int j = 0; j<order; ++j){
          if (t+1-j>=0){
            pred[i] += sample(t+1-j,k)*theta(j,i);
          }
        }
        // aggregate the result of the predictor
        agg_pred += alpha(i,k) * pred[i];        
      }
    }
    // store the results for this sample
    for(int i = 0; i<N; ++i){
      res(i,k) = pred[i];
    }
    res(N,k) = agg_pred;
  }
  List resu(2);
  resu[0] = res;
  resu[1] = alpha;
  return(resu);
}
