// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h> 

using namespace Rcpp;



// [[Rcpp::export()]]
Rcpp::List logistic_enet(Rcpp::NumericVector Yr, 
                                    Rcpp::NumericMatrix Xr,
                                    float lambda,
                                    float theta,
                                    int which_cov,
                                    float betaNull,
                                    Rcpp::NumericVector binitr,
                                    float delta){
  
  int n = Xr.nrow(), p = Xr.ncol()-1, i, k;
  float uj, vj, wj, sj;
  
  arma::mat X(Xr.begin(), n, p+1, false); 
  arma::colvec Y(Yr.begin(),Yr.size(), false);
  //arma::colvec gamma(gammar.begin(),gammar.size(),false);
  arma::colvec binit(binitr.begin(),binitr.size(),false);
  
  arma::colvec b1 = binit;
  arma::colvec b0 = binit;
  arma::colvec diff = arma::ones(p+1);
  
  arma::colvec b11 = binit;
  
  arma::colvec b00 = binit;
  arma::colvec diff2 = arma::ones(p+1);
  
  arma::colvec px(n);
  arma::colvec w(n);
  arma::colvec r(n);
  
  arma::colvec rj(n);
  arma::mat Xj = X;
  arma::colvec b11j(n);
  
  
  
  i = 0;
  
  while( (i < 500) & (diff.max() > delta))
  {
    b0 = b1;
    
    px = exp(X * b0) / (1 + exp(X * b0));
    w = px % (1 - px);
    r = X * b0 + ( Y - px) / w;
    
    k = 0;
    diff2 = arma::ones(p+1);
    
    while( (k < 500) & (diff2.max() > delta))
    {
      
      b00 = b11; 
      
      b11(0) = sum( w % ( r - X.cols(1,p) * b00.rows(1,p))) / sum(w);
      
      
      for(int j=1; j < (p+1) ; j++)
      {
        // if (j == which_cov) {continue;}
        Xj = X;
        Xj.shed_col(j);
        
        b11j = b11;
        b11j.shed_row(j);
        
        rj = Xj * b11j;
        uj = sum( w % ( r - rj) % X.col(j) );
        vj = theta * lambda; //* arma::as_scalar(gamma(j-1)) 
        wj = sum( w % pow(X.col(j),2)) + lambda * (1 - theta); 
        
        if (j == which_cov){
          
          if( (uj-betaNull*wj > 0) & (vj < std::abs(uj-betaNull*wj)) ) // soft-thresholding
          {
            sj = uj - vj;
            b11(j) = sj / wj;
            
          } else if( (uj-betaNull*wj < 0) & (vj < std::abs(uj-betaNull*wj)))
          {
            
            sj = uj + vj;
            b11(j) = sj / wj;
            
          } else {
            
            //sj = 0 + betaNull*wj;
            b11(j) = betaNull;
          }
          
        }else{
        
          if( (uj > 0) & (vj < std::abs(uj)) ) // soft-thresholding
          {
            sj = uj - vj;
          
          } else if( (uj < 0) & (vj < std::abs(uj)))
          {
          
            sj = uj + vj;
          
          } else {
          
            sj = 0;
          
          }
          
          b11(j) = sj / wj;
        }
        
        
        
      }
      
      
      diff2 = abs(b11 - b00);
      k++;
      
    }
    
    b1 = b11;
    
    
    diff = abs(b1 - b0);
    i++;
  }
  
  if(b1.has_nan()) 
  {				
    Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
    
  }
  
  return Rcpp::List::create(Named("b") = b1,
                            Named("lambda") = lambda,
                            Named("theta") = theta
  );
}



// [[Rcpp::export()]]
Rcpp::List logistic_enet_constraint(Rcpp::NumericVector Yr, 
                         Rcpp::NumericMatrix Xr,
                         float lambda,
                         float theta,
                         int which_cov,
                         float betaNull,
                         Rcpp::NumericVector binitr,
                         float delta){
  
  int n = Xr.nrow(), p = Xr.ncol()-1, i, k;
  float uj, vj, wj, sj;
  
  arma::mat X(Xr.begin(), n, p+1, false); 
  arma::colvec Y(Yr.begin(),Yr.size(), false);
  //arma::colvec gamma(gammar.begin(),gammar.size(),false);
  arma::colvec binit(binitr.begin(),binitr.size(),false);
  
  arma::colvec b1 = binit;
  arma::colvec b0 = binit;
  arma::colvec diff = arma::ones(p+1);
  
  arma::colvec b11 = binit;
  
  arma::colvec b00 = binit;
  arma::colvec diff2 = arma::ones(p+1);
  
  arma::colvec px(n);
  arma::colvec w(n);
  arma::colvec r(n);
  
  arma::colvec rj(n);
  arma::mat Xj = X;
  arma::colvec b11j(n);
  
  b1(which_cov) = betaNull;
  b11(which_cov) = betaNull;
  
 
  
  
  i = 0;
  
  while( (i < 500) & (diff.max() > delta))
  {
    b0 = b1;
    
    px = exp(X * b0) / (1 + exp(X * b0));
    w = px % (1 - px);
    r = X * b0 + ( Y - px) / w;
    
    k = 0;
    diff2 = arma::ones(p+1);
    
    while( (k < 500) & (diff2.max() > delta))
    {
      
      b00 = b11; 
      
      b11(0) = sum( w % ( r - X.cols(1,p) * b00.rows(1,p))) / sum(w);
      
      
      for(int j=1; j < (p+1) ; j++)
      {
        if (j == which_cov) {continue;}
        Xj = X;
        Xj.shed_col(j);
        
        b11j = b11;
        b11j.shed_row(j);
        
        rj = Xj * b11j;
        uj = sum( w % ( r - rj) % X.col(j) );
        vj = theta * lambda; //* arma::as_scalar(gamma(j-1)) 
        wj = sum( w % pow(X.col(j),2)) + lambda * (1 - theta); 
          
          if( (uj > 0) & (vj < std::abs(uj)) ) // soft-thresholding
          {
            sj = uj - vj;
            
          } else if( (uj < 0) & (vj < std::abs(uj)))
          {
            
            sj = uj + vj;
            
          } else {
            
            sj = 0;
            
          }
          
          
          b11(j) = sj / wj;
          
      }
      
      
      diff2 = abs(b11 - b00);
      k++;
      
    }
    
    b1 = b11;
    
    
    diff = abs(b1 - b0);
    i++;
  }
  
  if(b1.has_nan()) 
  {				
    Rcpp::Rcout << "warning: failure to converge due to complete or quasi-complete separation" << std::endl;
    
  }
  
  return Rcpp::List::create(Named("b") = b1,
                            Named("lambda") = lambda,
                            Named("theta") = theta
  );
}
