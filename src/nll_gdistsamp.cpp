#include <RcppArmadillo.h>
#include <float.h>
#include "distprob.h"
#include "distr.h"
#include "utils.h"

#ifdef _OPENMP
  #include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double nll_gdistsamp(arma::vec beta, arma::uvec n_param, arma::vec y,
    int mixture, std::string keyfun, std::string survey,
    arma::mat X_lambda, arma::vec offset_lambda, arma::vec A, arma::mat X_phi,
    arma::vec offset_phi, arma::mat X_det, arma::vec offset_det, arma::vec db,
    arma::mat a, arma::mat u, arma::vec w, arma::vec k, arma::vec lfac_k,
    arma::vec lfac_kmyt, arma::vec kmyt, arma::uvec Kmin, int threads){

  #ifdef _OPENMP
    omp_set_num_threads(threads);
  #endif

  int M = X_lambda.n_rows;
  int T = X_phi.n_rows / M;
  int R = y.size() / M;
  unsigned J = R / T;
  int lk = k.size();
  int K = lk - 1;

  //Abundance
  const vec lambda = exp(X_lambda * beta_sub(beta, n_param, 0) + offset_lambda) % A;
  double log_alpha = beta_sub(beta, n_param, 4)(0); //length 1 vector

  //Availability
  vec phi = ones(M*T);
  if(T > 1){
    phi = inv_logit(X_phi * beta_sub(beta, n_param, 1) + offset_phi);
  }

  //Detection
  vec det_param(M*T);
  if(keyfun != "uniform"){
    det_param = exp(X_det * beta_sub(beta, n_param, 2) + offset_det);
  }
  double scale = exp(beta_sub(beta, n_param, 3)(0));

  double loglik = 0.0;

  #pragma omp parallel for reduction(+: loglik) if(threads > 1)
  for (int i=0; i<M; i++){

    int t_ind = i * T;
    int y_ind = i * T * J;
    int k_start = i * T * lk;

    vec y_sub(J);
    vec p(J);
    vec p1(lk);
    vec p3(J);
    vec p4(lk);
    double p5;

    //Some unnecessary calculations here when k < Kmin
    //These values are ignored later in calculation of site_lp
    //However hard to avoid without refactoring entirely I think
    mat mn = zeros(K+1, T);
    for(int t=0; t<T; t++){
      int y_stop = y_ind + J - 1;
      y_sub = y.subvec(y_ind, y_stop);
      uvec not_missing = find_finite(y_sub);

      if(not_missing.size() == J){

        int k_stop = k_start + lk - 1;

        p1 = lfac_kmyt.subvec(k_start, k_stop);

        p = distprob(keyfun, det_param(t_ind), scale, survey, db,
                          w, a.row(i));
        p3 = p % u.col(i) * phi(t_ind);

        p4 = kmyt.subvec(k_start, k_stop);

        p5 = 1 - sum(p3);

        mn.col(t) = lfac_k - p1 + sum(y_sub % log(p3)) + p4 * log(p5);
      }

      t_ind += 1;
      y_ind += J;
      k_start += lk;
    }

    //Note that rows of mn for k < Kmin are skipped here
    double site_lp = 0.0;
    for (int j=Kmin(i); j<(K+1); j++){
      site_lp += N_density(mixture, j, lambda(i), log_alpha) *
        exp(sum(mn.row(j)));
    }

    loglik += log(site_lp + DBL_MIN);

  }

  return -loglik;

}
