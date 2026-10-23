#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// name of function below **MUST** match filename
template <class Type>
Type tmb_occuCT(objective_function<Type>* obj) {

  //Input data
  DATA_MATRIX(X_state);
  DATA_MATRIX(X_det);
  DATA_MATRIX(X_int);
  DATA_IVECTOR(dep_state);
  DATA_IMATRIX(dep_det);
  DATA_IMATRIX(dep_int);
  DATA_IVECTOR(no_detects);
  DATA_IVECTOR(dep_n);
  DATA_SCALAR(interval);

  //Parameters
  PARAMETER_VECTOR(beta_state);
  PARAMETER_VECTOR(beta_det);

  int M = X_state.rows();         // number of sites
  int nD = dep_state.size();      // number of deployments

  vector<Type> psi = invlogit(X_state * beta_state);

  // lambdas for actual observations
  vector<Type> lam_det = X_det * beta_det; // log(exp()) cancels out below
  // lambdas for integral calculation
  vector<Type> lam_int = exp(X_int * beta_det);

  Type ll = 0.0;
  vector<Type> lamsum_det(M);
  lamsum_det.setZero();
  vector<Type> integral(M);
  integral.setZero();
  vector<Type> lam_det_i;
  vector<Type> lam_int_i;

  // Iterate over deployments and add values to matching site
  for (int i=0; i<nD; i++){
    if(dep_n(i) > 0){
      //Likelihood of actual observations, if there are any
      lam_det_i = lam_det.segment(dep_det(i, 0), dep_det(i, 1));
      lamsum_det(dep_state(i)) += sum(lam_det_i);
    }
    //sum of integral approximation rectangles for this deployment
    lam_int_i = lam_int.segment(dep_int(i,0), dep_int(i,1));
    integral(dep_state(i)) += sum(lam_int_i) * interval;
  }

  // Iterate over sites
  for (int i=0; i<M; i++){
    ll += log(psi(i) * exp(lamsum_det(i)) * exp(-integral(i)) + (1-psi(i)) * no_detects(i));
  }
  
  return -ll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this
