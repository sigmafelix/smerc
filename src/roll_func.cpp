#include <RcppArmadillo.h>
using namespace arma;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::vec roll_psum(arma::vec &X, double p = 1.){
    // https://stackoverflow.com/questions/52459711/how-to-find-cumulative-variance-or-standard-deviation-in-r
    const arma::uword n_max = X.n_elem;
    arma::vec out(n_max);
    double *x = X.begin(), *o = out.begin();

    for(arma::uword n = 1; n <= n_max; ++n, ++x, ++o){
      *o += exp(log(*x) * p);
    }

    return out;
  }


// [[Rcpp::export]]
arma::vec roll_var(arma::vec &X){
    // https://stackoverflow.com/questions/52459711/how-to-find-cumulative-variance-or-standard-deviation-in-r
    const arma::uword n_max = X.n_elem;
    double xbar = 0, M = 0; 
    arma::vec out(n_max);
    double *x = X.begin(), *o = out.begin();

    for(arma::uword n = 1; n <= n_max; ++n, ++x, ++o){
      double tmp = (*x - xbar);
      xbar += (*x - xbar) / n;
      M += tmp * (*x - xbar);
      if(n > 1L)
        *o = M / (n - 1.);
    }

    if(n_max > 0)
      out[0] = std::numeric_limits<double>::quiet_NaN();

    return out;
  }