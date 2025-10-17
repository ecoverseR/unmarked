#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_colext(objective_function<Type>* obj) {

  DATA_VECTOR(y);
  DATA_MATRIX(X_psi);
  DATA_MATRIX(X_col);
  DATA_MATRIX(X_ext);
  DATA_MATRIX(X_det);
  DATA_INTEGER(M);
  DATA_INTEGER(T);
  DATA_INTEGER(J);
  DATA_MATRIX(site_sampled);
  DATA_MATRIX(nd);

  PARAMETER_VECTOR(beta_psi);
  PARAMETER_VECTOR(beta_col);
  PARAMETER_VECTOR(beta_ext);
  PARAMETER_VECTOR(beta_det);

  Type loglik = 0.0;

  vector<Type> psi = X_psi * beta_psi;
  psi = invlogit(psi);

  vector<Type> col = X_col * beta_col;
  col = invlogit(col);

  vector<Type> ext = X_ext * beta_ext;
  ext = invlogit(ext);

  vector<Type> p = X_det * beta_det;
  p = invlogit(p);

  int Tidx = 0;
  int pind = 0;

  vector<Type> Dy(2);
  matrix<Type> phi_total(2, 2);
  matrix<Type> phi(2, 2);
  vector<Type> phi_end(2);
  vector<Type> psi_vec(2);
  matrix<Type> Dy_diag(2,2);
  Dy_diag = Dy_diag.setZero();

  for (int i=0; i<M; i++){

    phi_total = phi.setIdentity();

    //All but the last time period
    for (int t=0; t<(T-1); t++){

      phi = phi.setZero();
      phi(0, 0) = 1 - ext(Tidx);
      phi(0, 1) = ext(Tidx);
      phi(1, 0) = col(Tidx);
      phi(1, 1) = 1 - col(Tidx);
      Tidx += 1;

      if(site_sampled(i, t) == 1){
        Dy = Dy.setOnes();

        if(nd(i, t) == 0){ // if at least one detection
          Dy(1) = 0;
        }

        for (int j=0; j<J; j++){
          if(R_IsNA(asDouble(y(pind)))){
            pind += 1;
            continue;
          }

          Dy(0) *= pow(p(pind), y(pind)) * pow(1-p(pind), 1-y(pind));
          pind += 1;
        }
        Dy_diag.diagonal() = Dy;
        //phi_total = phi_total * phi * Dy_diag;
        phi_total = phi_total * Dy_diag * phi;
      } else {
        phi_total = phi_total * phi;
        pind += J;
      }
    }
    
    //Final period; no psi calculation
    if(site_sampled(i, (T-1)) == 1){
      Dy = Dy.setOnes();

      if(nd(i, (T-1)) == 0){ // if at least one detection
        Dy(1) = 0;
      }

      for (int j=0; j<J; j++){
        if(R_IsNA(asDouble(y(pind)))){
          pind += 1;
          continue;
        }

        Dy(0) *= pow(p(pind), y(pind)) * pow(1-p(pind), 1-y(pind));
        pind += 1;
      }

      phi_end = phi_total * Dy;
    } else {
      Dy = Dy.setOnes();
      phi_end = phi_total * Dy;
      pind += J;
    }

    psi_vec(0) = psi(i);
    psi_vec(1) = 1 - psi(i);
    // Dot product
    vector<Type> site_lik = psi_vec.array() * phi_end.array();
    loglik += log(sum(site_lik));
  }

  return(-loglik);
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
